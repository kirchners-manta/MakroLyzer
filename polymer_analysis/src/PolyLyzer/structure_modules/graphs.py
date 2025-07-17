import networkx as nx
from scipy.spatial import cKDTree
import numpy as np
from IPython.display import display
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
from collections import defaultdict

from PolyLyzer import dictionaries

class GraphManager(nx.Graph):
    def __init__(self, data=None, boxSize=None, **kwargs):
        """
        Args:
            data : elements together with their coordinates.
            boxSize : optional, size of the box for periodic boundary conditions.
            **kwargs : additional keyword arguments for graph initialization.
        """
        
        super().__init__()
        
        if data is None:
            pass
        
        elif isinstance(data, GraphManager):
            # Copy constructor from another GraphManager instance
            self.add_nodes_from(data.nodes(data=True))
            self.add_edges_from(data.edges(data=True))
        elif isinstance(data, nx.Graph):
            # Initialize from a NetworkX graph object
            self.add_nodes_from(data.nodes(data=True))
            self.add_edges_from(data.edges(data=True))
        else:
            # Handle other types of data, such as initialization from raw data
            kwargs = {}
            if boxSize is not None:
                kwargs['boxSize'] = boxSize
            self.create_graph(data, **kwargs)

    
    def create_graph(self, atomData, boxSize=None):
        covalentRadii = dictionaries.dictCovalent()
        elements = atomData['atom'].values
        coords = atomData[['x','y','z']].values

        # Map each unique element to a small integer index
        # unique_elems returns a sorted array of its distinct values.
        # return_inverse=True returns an array of the same length as the original 'elements'.
        # Each entry of inverse is an integer index pointing into the unique_elements array,
        # tellingcwhich unique value the original element was.
        # EXAMPLE: elements=np.array(['H', 'C', 'O', 'H', 'C', 'N', 'O', 'O'])
        #          unique_elems=np.array(['C' 'H' 'N' 'O'])
        #          inverse=np.array([1, 0, 3, 1, 0, 2, 3, 3])
        unique_elems, inverse = np.unique(elements, return_inverse=True)
        M = len(unique_elems)

        # Build an M×M matrix of squared max distances
        maxd2 = np.zeros((M, M), dtype=float)
        for i, e1 in enumerate(unique_elems):
            for j, e2 in enumerate(unique_elems):
                d = (covalentRadii[e1] + covalentRadii[e2]) * 1.15
                maxd2[i, j] = d * d

        # Build KD-tree 
        if boxSize is not None:
            # wrap coordinates into (0, Lx], (0, Ly], (0, Lz] for periodic boundary conditions
            # since cKDTree wants coordinates in the range [0, boxSize)
            coords = np.mod(coords, boxSize)  
            # Create a cKDTree with periodic boundary conditions
            tree = cKDTree(coords, boxsize=boxSize)
            # shift coordinates to be within the periodic box
            coords = shift_coordinates(coords, boxSize)
        else:
            tree = cKDTree(coords)

        # Add all nodes 
        nodes = [
           (idx, {"index": idx, "element": elements[idx],
                  "x": coords[idx,0],
                  "y": coords[idx,1],
                  "z": coords[idx,2]})
           for idx in range(len(coords))
        ]
        self.add_nodes_from(nodes)

        # Find all candidate pairs up to the global max radius
        global_max = np.sqrt(maxd2.max())
        pairs = tree.query_pairs(global_max)

        # squared distances -> lookup in precomputed MxM matrix
        edge_list = []
        for i, j in pairs:
            ei, ej = inverse[i], inverse[j]
            # Get minimum image distance
            if boxSize is not None:
                dx, dy, dz = min_image_distance(coords[i], coords[j], boxSize)
            else:
                dx, dy, dz = coords[i] - coords[j]
            if dx*dx + dy*dy + dz*dz <= maxd2[ei, ej]:
                edge_list.append((i, j))
        # Add edges
        self.add_edges_from(edge_list)

        # Add node’s degree
        deg = dict(self.degree())
        for n, d in deg.items():
            self.nodes[n]['degree'] = d

            
    def get_coordinates(self, node):
        """
        Get the coordinates of a node in the graph.

        Args:
            node (int): The index of the node.
        
        Returns:
            np.ndarray: The coordinates of the node.
        """
        return np.array([self.nodes[node]['x'], self.nodes[node]['y'], self.nodes[node]['z']])
    
    def get_all_coordinates(self):
        """
        Get the coordinates of all nodes in the graph in a single (N, D) array.

        Returns:
            nodes (list): The list of node IDs, in the same order as rows of coords.
            coords (np.ndarray): Array of shape (N, 3) where N is the number of nodes.
        """
        nodes = list(self.nodes())                # ordering
        coords = np.stack([                       # (N,3) array
            self.get_coordinates(n) for n in nodes
        ], axis=0)
        return nodes, coords
    

    def remove_1order(self):
        """
        Remove 1-order nodes from the graph.

        Returns:
            GraphManager: A new GraphManager instance with 1-order nodes removed.
        """
        newGraph = GraphManager(self)
        for node in list(newGraph.nodes):
            if newGraph.nodes[node]['degree'] == 1:
                newGraph.remove_node(node)
        return newGraph
    
    def update_degree(self):
        """
        Update the degree of each node in the graph.

        Returns:
            GraphManager: The updated graph with new degrees.
        """
        for node in self.nodes:
            self.nodes[node]['degree'] = self.degree[node]
        return self
    
    def get_com(self):
        """
        Calculate the center of mass of the graph.

        Returns:
            np.ndarray: The coordinates of the center of mass.
        """
        total_mass = 0.0
        com = np.zeros(3)
        
        for node in self.nodes:
            element = self.nodes[node]['element']
            mass = dictionaries.dictMass().get(element, 1.0)
            total_mass += mass
            com += mass * np.array([self.nodes[node]['x'], self.nodes[node]['y'], self.nodes[node]['z']])
        if total_mass == 0:
            raise ValueError("Total mass is zero, cannot compute center of mass.")
        
        com /= total_mass
        return com           
    
    
    def surrounding(self):
        """
        Find the surrounding atoms of the graph.

        Returns:
            Modification of the graph: Adding a new attribute 'surroundingAtoms' to the graph.
        """
        for node in self.nodes:
            # Get the original atom label
            atomLabel = self.nodes[node]['element']
            
            # Get the neighbors of the atom
            neighbors = list(self.neighbors(node))
            neighborLabels = [self.nodes[neighbor]['element'] for neighbor in neighbors]
            # sort the labels alphabetically
            neighborLabels.sort()
            
            # Additional label
            addLabel = atomLabel + '_' + ''.join(neighborLabels)
            
            # Add additional label
            self.nodes[node]['surroundingAtoms'] = addLabel    
            
        return self
    
    def get_subgraphs(self):
        """
        Get the connected components of the graph.

        Returns:
            list of GraphManager: A list of GraphManager instances representing the connected components.
        """
        subgraphs = []
        for component in nx.connected_components(self):
            subgraph = GraphManager(self.subgraph(component))
            subgraphs.append(subgraph)
        return subgraphs
    
    
    def find_longest_path(self, startAtom=None):
        """
        Find the longest path in the graph using Depth First Search (DFS).

        Args:
            GraphManager: The graph to search.

        Returns:
            list of nodes: The longest path found in the graph.
        """
        def dfs(current, path, visited):
            visited.add(current)
            max_path = path

            for neighbor in self.neighbors(current):
                if neighbor == path[0] and len(path) > 1:
                    continue
                    
                if neighbor not in visited:
                    new_path = dfs(neighbor, path + [neighbor], visited)
                    if len(new_path) > len(max_path):
                        max_path = new_path

            visited.remove(current)
            return max_path

        longest = []
        # Do dfs for each node with degree == 1
        for node in self.nodes:
            if self.degree[node] != 1:
                continue
            # If the type of the start atom is provided, the longest path is
            # searched only from atoms of this type
            if startAtom is not None and self.nodes[node]['element'] != startAtom:
                continue
            # Perform DFS from the current node
            path = dfs(node, [node], set())
            if len(path) > len(longest):
                longest = path
                
        # Fallback for cycle graphs
        longest_cycle = []
        cycles = nx.cycle_basis(self, startAtom)
        if cycles:
            longest_cycle = max(cycles, key=len)
             
        if len(longest_cycle) > len(longest):
            longest = longest_cycle
                   
        return longest
    
    
    def distance(self, node1, node2):
        """
        Calculate the distance between two nodes in the graph.

        Args:
            node1 (int): The index of the first node.
            node2 (int): The index of the second node.
        Returns:
            float: The distance between the two nodes.
        """
        coords1 = np.array([self.nodes[node1]['x'], self.nodes[node1]['y'], self.nodes[node1]['z']])
        coords2 = np.array([self.nodes[node2]['x'], self.nodes[node2]['y'], self.nodes[node2]['z']])
        return np.linalg.norm(coords1 - coords2)

    def vector(self, node1, node2):
        """
        Calculate the vector between two nodes in the graph.

        Args:
            node1 (int): The index of the first node.
            node2 (int): The index of the second node.
        Returns:
            np.ndarray: The vector from node1 to node2.
        """
        coords1 = np.array([self.nodes[node1]['x'], self.nodes[node1]['y'], self.nodes[node1]['z']])
        coords2 = np.array([self.nodes[node2]['x'], self.nodes[node2]['y'], self.nodes[node2]['z']])
        return coords2 - coords1

    def dihedral(self, node1, node2, node3, node4, sign=None):
        """
        Calculate the dihedral angle between four nodes in the graph.

        Args:
            node1 (int): The index of the first node.
            node2 (int): The index of the second node.
            node3 (int): The index of the third node.
            node4 (int): The index of the fourth node.
        Returns:
            float: The dihedral angle in degrees.
        """

        # Get vectors
        b1 = self.vector(node1, node2)
        b2 = self.vector(node2, node3)
        b3 = self.vector(node3, node4)

        # Half planes (https://en.wikipedia.org/wiki/Dihedral_angle - In polymer physics)
        dihedral = np.arctan2(np.linalg.norm(b2) * np.dot(b1,(np.cross(b2,b3))), np.dot(np.cross(b1,b2),np.cross(b2,b3)))

        if sign is not None:
            dihedral = np.degrees(dihedral)
        elif sign is None:
            dihedral = np.abs(np.degrees(dihedral))

        return dihedral
    
    
    def get_hbonds(self, Acceptor, HAcceptor_dist, DonorAcceptor_dist, Angle_cut):
        H_type = "H"
        if Acceptor == H_type:
            raise ValueError("Element types must differ for H-bond detection.")

        # Collect H atoms and "other" atoms
        H_nodes, H_pos, Donor_pos, H_vectors = [], [], [], []
        type2_nodes, type2_pos = [], []
        for node, data in self.nodes(data=True):
            pos = np.array([data['x'], data['y'], data['z']])
            if data['element'] == H_type:
                # H is only considered when bound to a heavy atom
                neighbor = list(self.neighbors(node))
                if len(neighbor) == 1 and self.nodes[neighbor[0]]['element'] != "C":
                    H_nodes.append(node)
                    H_pos.append(pos)
                    Donor_pos.append(self.get_coordinates(neighbor[0]))
                    # vector of H node to covalently bound atom
                    vector = np.asarray(self.vector(node, neighbor[0]), dtype=float)
                    H_vectors.append(vector / np.linalg.norm(vector))
                elif len(neighbor) == 0:
                    raise ValueError(f"H atom is not covalently bound.")
                elif len(neighbor) > 1:
                    raise ValueError(f"H atom is covalently bound to multiple atoms.")
            elif data['element'] == Acceptor:
                type2_nodes.append(node)
                type2_pos.append(pos)

        # exit if no H's or no potential acceptors are found
        # -> no H-bonds can be formed
        if not H_nodes or not type2_nodes:
            return []

        # vertically stack the positions
        H_pos     = np.vstack(H_pos)      # shape (n_H, 3)
        type2_pos = np.vstack(type2_pos)  # shape (n_other, 3)

        # Build KD-tree on acceptor positions, because then we can query
        # for all candidates within HAcceptor_dist in O(log n) time.
        # -> we dont have to check all pairs
        tree = cKDTree(type2_pos)

        # Query for all H-bonds
        hbonds = []
        for i, nodeH in enumerate(H_nodes):
            # candidate indices within HAcceptor_dist
            # query_ball_point(x, r)
            # x : The point or points to search for neighbors of.
            # r : The radius within which to search for neighbors.
            # -> returns indices of all neighbors within radius r of x
            idxs_HA = tree.query_ball_point(H_pos[i], r=HAcceptor_dist) 
            idxs_DA = tree.query_ball_point(Donor_pos[i], r=DonorAcceptor_dist)
            # Only intersect the two sets of indices
            idxs = set(idxs_HA).intersection(idxs_DA)       
            for j in idxs:
                node2 = type2_nodes[j]
                # skip if already covalently bonded
                if self.has_edge(nodeH, node2):
                    continue
                # Convert the input to an array.
                # vector from H to Type2 atom
                HType2_vector = np.asarray(self.vector(nodeH, node2), dtype=float)
                HType2_vector /= np.linalg.norm(HType2_vector)
                cosang = np.dot(H_vectors[i], HType2_vector)
                angle = np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))
                if angle < Angle_cut or angle > (180-Angle_cut):
                    hbonds.append((nodeH, node2))

        return hbonds
    
    def get_vectors_and_positions_along_path(self, path, unitSize=None):
        """
        Get vectors between nodes in a path. The number of atoms that correspond to a 
        vector is specified by the unitSize parameter.

        Args:
            path in the form of a list of node indices.
            unitSize : size of the unit vector, i.e. how many atoms correspond to a vector.
                       If None, the unitSize is set to 1.
        Returns:
            Dict: A dictionary with keys 'vectors' and 'positions'.
                   'vectors' is a list of vectors between nodes in the path.
                   Each vector is given as a tuple of (x, y, z) coordinates, corrsponding to the
                   position of the nodes corresponding to the vector.
        """
        if unitSize is None:
            unitSize = 2
        elif unitSize < 2:
            raise ValueError("unitSize must be greater than or equal to 2.")
        elif unitSize > len(path):
            unitSize = len(path)
            
        # Adjust unitSize to be the number of atoms that correspond to a vector
        unitSize = unitSize - 1  
            
        vecToPos = defaultdict(list)
        coords = np.array([self.get_coordinates(node) for node in path])
        
        # Iterate over the path and calculate vectors
        for i in range(0, len(path) - unitSize, unitSize):
            # Get coords of the start and end nodes and calculate the vector
            start_coords, end_coords = coords[i], coords[i + unitSize]
            vector = end_coords - start_coords
            norm = np.linalg.norm(vector)
            if norm == 0:
                continue
            unitVector = vector / norm
            midpoint = tuple(start_coords + vector / 2.0)
            
            vecToPos[midpoint].append(unitVector)
            
        return dict(vecToPos)
    
    
    def find_and_tag_patterns(self, patterns, startAtom=None):
        """
        Find and tag specific patterns in the graph.

        Args:
            patterns (list of tuples): List of tuples representing the patterns to find.

        Returns:
            list of tuples: List of tuples representing the found patterns.
        """
        # Remove 1-order nodes, find subgraphs and surrounding atoms
        GraphWithout1order = self.remove_1order()
        GraphWithout1order.surrounding()
        GraphWithout1order.update_degree()
        subgraphs = GraphWithout1order.get_subgraphs()
        
        # Dictionary to store the fragment IDs
        fragmentIDs = {}
        fragmentID = 0
        
        
        # Iterate through the subgraphs and find the patterns
        for subgraph in subgraphs:
            longestPath = subgraph.find_longest_path(startAtom)

            # Add side chains to the longest path
            for node in longestPath:
                for neighbor in GraphWithout1order.neighbors(node):
                    if neighbor not in longestPath:
                        longestPath.insert(longestPath.index(node) + 1, neighbor)
                    
            # Create a list of labels for the longest path
            labels = [subgraph.nodes[node]['surroundingAtoms'] for node in longestPath]
            
            # Slide over the longest path and check for patterns
            for pattern in patterns:
                i = 0
                while i <= len(labels) - len(pattern):
                    # Check if the current window matches the pattern
                    if labels[i:i + len(pattern)] == pattern:
                        # Check if one of the nodes in the window is already tagged
                        if any(longestPath[i + j] in fragmentIDs for j in range(len(pattern))):
                            # If so, skip this window
                            i += 1
                            continue
                        
                        # Tag the nodes with the current fragment ID 
                        for j in range(len(pattern)):                       
                            fragmentIDs[longestPath[i + j]] = fragmentID
                        fragmentID += 1
                        i += len(pattern)
                    else:
                        i += 1
                        
            # Add the fragment IDs to the nodes in self
            for node in longestPath:
                if node in fragmentIDs:
                    self.nodes[node]['fragmentID'] = fragmentIDs[node]
                    self.fragents_for_neighbor(longestPath, node)
                    # Neighbors get the same fragment ID
                    #for neighbor in self.neighbors(node):
                    #    if self.degree[neighbor] == 1:
                    #        self.nodes[neighbor]['fragmentID'] = fragmentIDs[node]
                else:
                    self.nodes[node]['fragmentID'] = -1
                    
                    
    def fragents_for_neighbor(self, longestPath, node):
        """
        Assign fragment IDs to the neighbors of the given node.

        Args:
            node: The index of the node to assign fragment IDs to.

        Returns:
            None
        """
        # Get the fragment ID of the node
        fragmentID = self.nodes[node]['fragmentID']
        
        # Iterate over the neighbors of the node
        for neighbor in self.neighbors(node):
            # If the neighbor has a fragment ID of -1, or no fragment ID,
            # assign the fragment ID of the node to the neighbor
            if  'fragmentID' not in self.nodes[neighbor] and neighbor not in longestPath:
                self.nodes[neighbor]['fragmentID'] = fragmentID
                # If neighbor is not degree 1, recursively call the function
                if self.degree[neighbor] != 1:
                    self.fragents_for_neighbor(longestPath, neighbor)
                    
    def saturate(self):
        """
        Saturate the Polymers by adding atoms at the ends of the longest path.

        Returns:
            GraphManager: A new GraphManager instance with saturated Polymers.
        """
        GraphWithout1order = self.remove_1order()
        GraphWithout1order.surrounding()
        GraphWithout1order.update_degree()
        subgraphs = GraphWithout1order.get_subgraphs()
        
        # Iterate through the subgraphs and find the patterns
        for subgraph in subgraphs:
            longestPath = subgraph.find_longest_path()
            
            # find nodes in longest path with degree 1
            start_nodes = [node for node in longestPath if GraphWithout1order.nodes[node]['degree'] == 1]
            
            # Add OH to node if C atom
            for startNode in start_nodes:
                if startNode in self.nodes and self.nodes[startNode]['element'] == 'C':
                    # Add OH to the start node
                    self.add_OH_to_carboxylic_acid(startNode)
                elif startNode in self.nodes and self.nodes[startNode]['element'] == 'N':
                    # Add H to the start node
                    self.add_H_to_amide(startNode)
        # return the modified graph
        return self
            
    def add_OH_to_carboxylic_acid(self, node):
        """
        Add OH to a carboxylic acid node (assumed to be the carbon in -COOH).
        """
        import numpy as np

        coords = np.array([self.nodes[node]['x'], self.nodes[node]['y'], self.nodes[node]['z']])
        
        neighbor_C = None
        neighbor_O = None

        # Find the neighbors of the C atom
        # -> C and O
        for neighbor in self.neighbors(node):
            element = self.nodes[neighbor]['element']
            pos = np.array([self.nodes[neighbor]['x'], self.nodes[neighbor]['y'], self.nodes[neighbor]['z']])
            if element == 'C':
                neighbor_C = pos
            elif element == 'O':
                neighbor_O = pos
        if neighbor_C is None or neighbor_O is None:
            raise ValueError("Carboxylic acid carbon must have one C and one O neighbor.")
        
        # calculate the coordinates of the new O atom
        coords_mid = (neighbor_O - neighbor_C) / 2 + neighbor_C 
        length_OC_single = 1.34
        vecOH = (coords_mid - coords)/np.linalg.norm(coords_mid - coords) * length_OC_single
        coords_OH = coords - vecOH
        
        # Add the new oxygen atom
        new_index = max(self.nodes) + 1
        self.add_node(new_index, index=new_index, element='O', x=coords_OH[0], y=coords_OH[1], z=coords_OH[2])
        self.add_edge(node, new_index)
        
        # calculate the coordinates of the new H atom
        length_OH = 0.96
        vec_CC = (coords - neighbor_C)/np.linalg.norm(coords - neighbor_C)
        coords_ = coords_OH + vec_CC * length_OH
        
        # Add the new hydrogen atom
        new_index = max(self.nodes) + 1
        self.add_node(new_index, index=new_index, element='H', x=coords_[0], y=coords_[1], z=coords_[2])
        self.add_edge(new_index, new_index - 1)
        
    def add_H_to_amide(self, node):
        """
        Add H to an amide node (assumed to be the nitrogen in -NH2).        
        """
        coords = np.array([self.nodes[node]['x'], self.nodes[node]['y'], self.nodes[node]['z']])
        neighbor_C = None
        
        # Find the neighbors of the N atom
        # -> C, H
        for neighbor in self.neighbors(node):
            element = self.nodes[neighbor]['element']
            pos = np.array([self.nodes[neighbor]['x'], self.nodes[neighbor]['y'], self.nodes[neighbor]['z']])
            if element == 'C':
                # remember the C atom
                Cnode = neighbor
                neighbor_C = pos
            elif element == 'H':
                neighbor_H = neighbor
        if neighbor_C is None:
            raise ValueError("Amide nitrogen must have one C neighbor.")
        
        # Get the two H atoms coordinated to the N atom
        C_Hnodes = []
        for neighbor in self.neighbors(Cnode):
            element = self.nodes[neighbor]['element']
            pos = np.array([self.nodes[neighbor]['x'], self.nodes[neighbor]['y'], self.nodes[neighbor]['z']])
            if element == 'H':
                # remember the H atom
                C_Hnodes.append(neighbor)
        if len(C_Hnodes) != 2:
            raise ValueError("Amide carbon must have two H neighbors.")
        
        # coords of the H atoms
        coords_H1 = np.array([self.nodes[C_Hnodes[0]]['x'], self.nodes[C_Hnodes[0]]['y'], self.nodes[C_Hnodes[0]]['z']])
        coords_H2 = np.array([self.nodes[C_Hnodes[1]]['x'], self.nodes[C_Hnodes[1]]['y'], self.nodes[C_Hnodes[1]]['z']])
        
        # get C-H vectors
        vec_C_H1 = (neighbor_C - coords_H1)/np.linalg.norm(neighbor_C - coords_H1)
        vec_C_H2 = (neighbor_C - coords_H2)/np.linalg.norm(neighbor_C - coords_H2)
        
        # calculate the coordinates of the new H atoms
        length_NH = 1.02
        coords_new_H1 = coords + vec_C_H1 * length_NH
        coords_new_H2 = coords + vec_C_H2 * length_NH
        
        # Update the coordinates of the neighbor H atom
        self.nodes[neighbor_H]['x'] = coords_new_H1[0]
        self.nodes[neighbor_H]['y'] = coords_new_H1[1]
        self.nodes[neighbor_H]['z'] = coords_new_H1[2]
        # Add the new hydrogen atom
        new_index = max(self.nodes) + 1
        self.add_node(new_index, index=new_index, element='H', x=coords_new_H2[0], y=coords_new_H2[1], z=coords_new_H2[2])
        self.add_edge(node, new_index)
        
        #if np.linalg.norm(coords_new_H1 - neighbor_H) > 0.8:
        #    new_index = max(self.nodes) + 1
        #    self.add_node(new_index, index=new_index, element='H', x=coords_new_H1[0], y=coords_new_H1[1], z=coords_new_H1[2])
        #    self.add_edge(node, new_index)
        #elif np.linalg.norm(coords_new_H2 - neighbor_H) > 0.8:
        #    new_index = max(self.nodes) + 1
        #    self.add_node(new_index, index=new_index, element='H', x=coords_new_H2[0], y=coords_new_H2[1], z=coords_new_H2[2])
        #    self.add_edge(node, new_index)
        #else:
        #    raise ValueError("No new H atom can be added to the amide nitrogen.")
        
    def get_chemicalFormulas(self):
        """
        Get the chemical formulas of the polymers.

        Returns:
            dict: A dictionary with the chemical formulas of the polymers and their counts.
        """
        # Get the subgraphs of the graph
        subgraphs = self.get_subgraphs()
        formulas = {}
        
        # Iterate through the subgraphs and count the elements
        for subgraph in subgraphs:
            formula = {}
            for node in subgraph.nodes:
                element = subgraph.nodes[node]['element']
                if element not in formula:
                    formula[element] = 0
                formula[element] += 1
            
            # sort alphabetically
            formula = dict(sorted(formula.items()))
            # Convert the formula to a string
            formula_str = ''.join([f"{k}{v}" for k, v in formula.items()])
            formulas[formula_str] = formulas.get(formula_str, 0) + 1
            
        # Sort the formulas by their counts
        formulas = dict(sorted(formulas.items(), key=lambda item: item[1], reverse=True))   
        # Convert the counts to a list of tuples
        formulas = [(k, v) for k, v in formulas.items()]
        
        return formulas
    
                
    def draw_graph(self):
        """Draws the largest connected component of the graph."""
        subgraphs = self.get_subgraphs()
        largestSubgraph = max(subgraphs, key=len)
        
        # 3D visualization
        pos = {node: (self.nodes[node]['x'], self.nodes[node]['y'], self.nodes[node]['z']) for node in largestSubgraph.nodes}
        
        # Scatter plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(*zip(*pos.values()), c='b', marker='o')
        
        # Draw edges
        for edge in largestSubgraph.edges:
            x = [self.nodes[edge[0]]['x'], self.nodes[edge[1]]['x']]
            y = [self.nodes[edge[0]]['y'], self.nodes[edge[1]]['y']]
            z = [self.nodes[edge[0]]['z'], self.nodes[edge[1]]['z']]
            ax.plot(x, y, z, color='gray', alpha=0.5)
        
        # Add labels and indices
        for node, (x, y, z) in pos.items():
            label = f"{self.nodes[node]['element']}_{node}"
            label += '\n' + str(self.nodes[node]['surroundingAtoms'])
            ax.text(x, y, z, label, color='black', fontsize=6)
                
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        plt.title('3D Scatter Plot of Graph Nodes')
        plt.show()
        
    def write_xyz(self, file):
        """
        Write the graph data to an XYZ file.
        """
        with open(file, 'w') as f:
            f.write(f"{len(self.nodes)}\n")
            f.write("Graph data\n")
            for node in self.nodes:
                element = self.nodes[node]['element']
                x = self.nodes[node]['x']
                y = self.nodes[node]['y']
                z = self.nodes[node]['z']
                f.write(f"{element} {x} {y} {z}\n")


    def write_fragment_data_to_csv(self, file):
        """
        Write the fragment data to a CSV file.

        Args:
            filename (str): The name of the output CSV file.
        """
        with open(file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=';')
            writer.writerow(['Index', 'Element', 'X', 'Y', 'Z', 'FragmentID'])

            for node in self.nodes:
                data = self.nodes[node]
                writer.writerow([
                    data['index'],
                    data['element'],
                    data['x'],
                    data['y'],
                    data['z'],
                    data['fragmentID']
                ])



def min_image_distance(pos1, pos2, boxSize) -> tuple:
    """
    Calculate the minimum image distance between two positions in a periodic box.
    
    Args:
        pos1 (np.ndarray): The first position.
        pos2 (np.ndarray): The second position.
        boxSize (float or np.ndarray): The size of the box.
        
    Returns:
        float: The minimum image distance.
    """
    
    delta = pos1 - pos2
    delta -= np.round(delta / boxSize) * boxSize
    return delta[0], delta[1], delta[2]



def shift_coordinates(coords: np.ndarray, box_size):
    """
    Shift coordinates to be within the periodic box.

    Args:
        coords (np.ndarray): The coordinates to shift, shape (N, 3).
        box_size (float | np.ndarray): Box length (scalar) or per-axis lengths.

    Returns:
        np.ndarray: A copy of the coordinates with (hopefully) no particles split across the periodic boundary.
    """
    # Ensure that the particle is not split across the periodic boundary
    # For each direction (x,y,z), get the minimum and maximum coordinates
    coords   = np.asarray(coords, dtype=float)
    box_size = np.asarray(box_size, dtype=float)

    # So far only one boxSize can be given, so we repeat it for all three axes
    if box_size.ndim == 0:
        box_size = np.repeat(box_size, 3)

    shifted = coords.copy()                
    lower   = 0.05 * box_size              # 5 % threshold
    upper   = 0.95 * box_size              # 95 % threshold
    step    = 0.10 * box_size              # 10 % shift step

    # Iterate over the three Cartesian axes
    for ax in range(3):
        # Do that until the minimum is >5 % and the maximum is <95 % of the box size,
        # or we try 9 times (shift of 90%)
        for i in range(10):
            min_i, max_i = shifted[:, ax].min(), shifted[:, ax].max()
            if not (min_i < lower[ax] and max_i > upper[ax]):
                break
            # Shift the coordinates up by 10 % of the box size in this direction
            # and put them back into the box using modulo
            shifted[:, ax] = np.mod(shifted[:, ax] + step[ax], box_size[ax])

    return shifted
