import numpy as np
from scipy.spatial import cKDTree
from collections import defaultdict

from MakroLyzer.structure_modules.structureBase import StructureAnalyzer
from MakroLyzer.math import Box
from MakroLyzer.errorOutputs.ErrorOutputs import ErrorOutputs

class OrderParameterAnalyzer(StructureAnalyzer):
    """
    Analyzer for the calculation of the order parameter S* of a given molecular graph.
    Inherits from StructureAnalyzer.
    
    S* is a measure of the the degree of orientation of molecules/units with respect to each other.
    The order parameter S* is closely related to the nematic order parameter S used in liquid crystal physics.
    S is defined as:
    
    S = 1/2 <(3 * cos²(θ) - 1)>
    
    where θ is the angle between the molecular orientation vectors and a reference direction.
    The reference direction is typically the "average orientation" of the molecules in the system,
    and <> denotes the average over all molecules.
    The expression "average orientation" is not fullycorrect, as the direction of the vectors is arbitrary.
    S is usually calculated from the Q tensor, which is defined as:
    
    Q = S[(3/2)(n⊗n) - (1/2)I]
    
    where S is the order parameter, n is the unit orientation vector (molecular vector), and I is the identity matrix.
    
    Since n is unknown, we calculate the helper tensor Q', which is defined as:
    
    Q' = (1/N) Σ[(3/2)(u⊗u) - (1/2)I]
    
    where u are the unit orientation vectors of the molecules, which are determined along the backbone of 
    a given molecule, and N is the number of molecules.
    Q and Q' share the same eigenvalues, so we can calculate S* from the largest eigenvalue of Q'.
    S and S* differ only for negative values for S.
    S ∈ [-0.5,1], while S* ∈ [0,1]. 
    
    A value of S* = 1 indicates perfect alignment of all molecules along a common direction, while
    S* = 0 indicates isotropic orientation.
    
    The order parameter S* can either be calculated for the entire graph, or the space can be divided into
    smaller cells, and S* can be calculated for each cell individually.
    The overall order parameter is then the average of the order parameters of all cells.
    """
    
    def __init__(self, BoxSize, NoCellsPerDim, MolecularVectorLength, output_handler=None, static_topology=False, backbone_cache=None):
        """
        Initialize the OrderParameterAnalyzer.

        Args:
            BoxSize (tuple of 3 floats): Length of the box in x,y and z directions. If one is given, the same box length
                                         is asssumed in all directions.
            NoCellsPerDim (tuple of 3 int): Number of cells in each direction, to calculate the order parameter for.
                                            If only one is given, the same number of cells is assumed in all directions.
            MolecularVectorLength (int): Number of atoms that are used to calculate a molecular vector from.
            output_handler (OutputHandler): Handler for writing output.
        """
        # Check if the given parameters have the correct format
        if isinstance(BoxSize, (int, float)):
            BoxSize = (BoxSize, BoxSize, BoxSize)
        elif isinstance(BoxSize, tuple) and len(BoxSize) == 3:
            pass
        else:
            raise ValueError(ErrorOutputs.WRONG_INPUT_TYPE_OP_ERROR)
        
        if isinstance(NoCellsPerDim, int):
            NoCellsPerDim = (NoCellsPerDim, NoCellsPerDim, NoCellsPerDim)
        elif isinstance(NoCellsPerDim, tuple) and len(NoCellsPerDim) == 3:
            pass
        else:
            raise ValueError(ErrorOutputs.WRONG_INPUT_TYPE_OP_ERROR)
        
        if not isinstance(MolecularVectorLength, int) or MolecularVectorLength < 2:
            raise ValueError(ErrorOutputs.WRONG_INPUT_TYPE_OP_ERROR)
        
        super().__init__(output_handler)
        self.BoxSize = BoxSize
        self.NoCellsPerDim = NoCellsPerDim
        self.MolecularVectorLength = MolecularVectorLength
        self.static_topology = static_topology
        self._cached_paths = None
        self.backbone_cache = backbone_cache
        
    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler:
            header = "Frame, Order Parameter S*"
            if self.output_handler.mode == 'streaming':
                self.output_handler.initialize_file(header)
                
    def _compute_backbone_paths(self, graph):
        """
        Compute longest backbone paths per subgraph on a reduced graph.
        Returns a list of paths (lists of node indices).
        """
        newGraph = graph.remove_1order()
        newGraph.update_degree()
        subgraphs = newGraph.get_subgraphs()

        paths = []
        for subgraph in subgraphs:
            longestPath = subgraph.find_longest_path()
            if len(longestPath) < 2:
                continue
            paths.append(longestPath)

        return paths

    def get_backbone_vectors(self, graph):
        """
        Calculate vectors with length 1 from the atoms with number of MolecularVectorLength
        consecutive atoms/nodes along the backbone.
        
        Args:
            graph (GraphManager): The graph to calculate the vectors for.
            MolecularVectorLength (int): The number of atoms/nodes that define the length of a unit.
        """
        vecAndPos = defaultdict(list)
        if self.backbone_cache is not None:
            paths = self.backbone_cache.get_paths(graph)
            for path in paths:
                pathDict = graph.get_vectors_and_positions_along_path(path, self.MolecularVectorLength)
                for midpoint, vector in pathDict.items():
                    vecAndPos[midpoint].extend(vector)
        elif self.static_topology:
            if self._cached_paths is None:
                self._cached_paths = self._compute_backbone_paths(graph)
            paths = self._cached_paths
            for path in paths:
                # Use current graph coordinates for vectors
                pathDict = graph.get_vectors_and_positions_along_path(path, self.MolecularVectorLength)
                for midpoint, vector in pathDict.items():
                    vecAndPos[midpoint].extend(vector)
        else:
            # Prepare the graph (non-static topology)
            newGraph = graph.remove_1order()
            newGraph.update_degree()
            subgraphs = newGraph.get_subgraphs()
            for subgraph in subgraphs:
                longestPath = subgraph.find_longest_path()
                if len(longestPath) < 2:
                    continue
                pathDict = subgraph.get_vectors_and_positions_along_path(longestPath, self.MolecularVectorLength)

                # append pathDict to vecAndPos
                # pathDict is a dictionary with: keys = midpoints and values = vectors
                for midpoint, vector in pathDict.items():
                    vecAndPos[midpoint].extend(vector)

        return dict(vecAndPos)
    
    def get_backbone_vectors_in_cell(self, cell, VecAndPos, midpoints, tree):
        """
        Get all vectors that are in a specific cell.

        Args:
            cell (tuple): The cell to check for vectors with format [(x_start, x_end), (y_start, y_end), (z_start, z_end)].
            vecAndPos (dict): keys = midpoints, values = vectors.
                              The midpoints are the positions of the vectors.
        Returns:
            list: list with all vectors that are in the cell.
                  If a vector has multiple positions, its returned multiple times.
        """  
        # Unpack the cell boundaries
        (x0, x1), (y0, y1), (z0, z1) = cell

        vectorsInCell = []

        # Query the tree for all midpoints within the cell
        # the midpoint of the cells is taken as the center for the query
        cellCenter = ((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2)
        # and the radius is half the size of the cell in the largest dimension
        halfCellDiagonal = np.linalg.norm([(x1 - x0) / 2, (y1 - y0) / 2, (z1 - z0) / 2])
        # This ensures that we get all midpoints that are within the cell
        # we obtain the indices of the midpoints that are potentially in the cell
        indices = tree.query_ball_point([cellCenter], r=halfCellDiagonal)

        # we need to filter the midpoints that are actually in the cell
        for idx in indices[0]:
            midpoint = midpoints[idx]
            # check if the midpoint is within the cell boundaries
            if x0 <= midpoint[0] < x1 and y0 <= midpoint[1] < y1 and z0 <= midpoint[2] < z1:
                # if yes, we add the corresponding vectors to the list
                vectorsInCell.extend(VecAndPos[tuple(midpoint)])

        return vectorsInCell
    
    def Q_prime_tensor(self, vectors):
        """
        Calculate the Q'-tensor for a given unit orientation vector and a set of vectors.

        Q' = 1/N Σ (3/2 * (u_i ⊗ u_i) - 1/2 * I)

        where u_i is the unit vector of each vector in the set, I is the identity matrix,
        and N is the number of vectors.

        Args:
            unitOrientationVector : np.ndarray
                                    The unit vector representing the orientation of the specific cell.
            vectors : np.ndarray
                      An array of vectors within the cell.
        Returns:
            np.ndarray : The Q-tensor.
        """
        if len(vectors) == 0:
            return np.full((3, 3), np.nan)

        N = len(vectors)
        Q_prime = np.zeros((3, 3))

        for vec in vectors:
            u_i = vec / np.linalg.norm(vec)
            Q_prime += (3/2) * np.outer(u_i, u_i) - (1/2) * np.eye(3)
        Q_prime /= N

        return Q_prime
                
    def compute(self, graph):
        """
        Calculate the order parameter S* from the Q' tensor for the graph within a given box.

        Args:
            graph (GraphManager): Graph to analyze.
            
        Returns: 
            Average order parameter S* (float).
        """
        # Divide the box into cells
        cells = Box.Box.divideBox(self.BoxSize, self.NoCellsPerDim)
        
        # Obtain dictionary of molecular vectors of length 1 along the backbone
        # together with their corresponding point in space.
        VecAndPos =  self.get_backbone_vectors(graph)
        
        # Build a cKDtree from the midpoints of the vectors for fast spatial queries 
        midpoints = np.asarray(list(VecAndPos.keys()), dtype=float)
        tree = cKDTree(midpoints)
        
        OrderParameters = []
        for cell in cells:
            # Get the vectors inside of the cell
            vectors = self.get_backbone_vectors_in_cell(cell, VecAndPos, midpoints, tree)
            # If there are more than two vectors inside of the cell, calculate the order parameter
            if len(vectors) < 3:
                continue
            Q_prime = self.Q_prime_tensor(vectors)
            Q_prime_eigVec,  _ = np.linalg.eig(Q_prime)
            S_star = np.max(Q_prime_eigVec)
            
            OrderParameters.append(S_star)
            
        # Average the order parameter over all cells
        if OrderParameters:
            return np.mean(OrderParameters)
        else:
            return np.nan
        
    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame.
        
        Args:
            data (list): The computed S* data.
            frame_idx (int): Current frame number.
        """
        row = (
            f"{frame_idx},"
            f"{data:.3f}"
        )
        self.output_handler.append_row(row)
        
    def finalize_output(self, header=None):
        """
        Finalize output file (write header and rows - 'collect' mode)
        """
        header = "Frame, Order Parameter S*"
        super().finalize_output(header)
