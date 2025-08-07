import numpy as np 
from collections import defaultdict
from scipy.spatial import cKDTree

from MakroLyzer.math import Box
from MakroLyzer.structure_modules import graphs
import time

def get_unit_vectors(graph, unitSize):
    """
    Calculate vectors with corresponding positions along all longest paths of the graph.
    unitSize defines the length of a unit (in nodes) for the vectors.
    
    Args:
        graph (GraphManager): The graph to calculate the vectors for.
        unitSize (int): The number of nodes that define the length of a unit.
    Returns:
        dict: A dictionary with all vectors and their positions.  
    """
    
    # Prepare the graph
    newGraph = graph.remove_1order()
    newGraph.update_degree()
    
    subgraphs = newGraph.get_subgraphs()
    
    vecAndPos = defaultdict(list)
    for subgraph in subgraphs:
        longestPath = subgraph.find_longest_path()
        if len(longestPath) < 2:
            continue
        
        # Get dictionary for the longest path of the subgraph
        pathDict = subgraph.get_vectors_and_positions_along_path(longestPath, unitSize)
        
        # append pathDict to vecAndPos
        # pathDict is a dictionary with: keys = midpoints and values = vectors
        for midpoint, vector in pathDict.items():
            vecAndPos[midpoint].extend(vector)
        
    return dict(vecAndPos)

def get_unit_vectors_in_cell(cell, vecAndPos, midpoints, tree):
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
            vectorsInCell.extend(vecAndPos[tuple(midpoint)])
          
    return vectorsInCell

def get_unit_orientation(vectors):
    """
    Calculate the "molecular" orientation of the given vectors. 
    
    Args:
        vectors : np.ndarray
                An array of vectors for which the orientation should be calculated.
    Returns:
        vector : np.ndarray
                The average orientation vector of the given vectors.
    """
    
    if len(vectors) == 0:
        return np.array([np.nan, np.nan, np.nan])
    
    # Calculate the average vector
    vector = np.mean(vectors, axis=0)
    
    # Normalize the vector
    norm = np.linalg.norm(vector)
    if norm == 0:
        return np.array([np.nan, np.nan, np.nan])
    
    return vector / norm

def order_parameter(unitOrientationVector, vectors):
    """
    Calculate the order parameter for a given unit orientation vector and a set of vectors.
    
    S = 1/2 <(3 * cos²(θ) - 1)>
    
    where θ is the angle between the unit orientation vector and each vector in the set.
    <> denotes the average over all vectors.
    
    Args:
        unitOrientationVector : np.ndarray
                                The unit vector representing the orientation of the specific cell.
        vectors : np.ndarray
                  An array of vectors within the cell.
    Returns:
        float : The order parameter value.
    """
    
    if len(vectors) == 0:
        return np.nan
    
    # Calculate the cosine of the angle between the unit orientation vector and each vector
    cos_theta = np.dot(vectors, unitOrientationVector)
    
    # Calculate the order parameter
    order_param = 0.5 * np.mean(3 * cos_theta**2 - 1)
    
    return order_param

def get_order_parameter(graph, boxSize, n, unitSize):
    """
    Calculate the order parameter for the given graph within a divided box.
    
    Args:
        graph (GraphManager): The graph to calculate the order parameter for.
        boxSize (tuple): The size of the box in each dimension (x, y, z).
        n (tuple): The number of divisions along each dimension (nx, ny, nz).
        unitSize (int): The number of nodes that define the length of a unit.
    Returns:
        list: A list of order parameters for each sub-box.
    """
        
    # Divide the box into cells
    cells = Box.Box.devideBox(boxSize, n)
    
    # Get unit vectors for the graph
    vecAndPos = get_unit_vectors(graph, unitSize)
    
    # build a cKDTree from the midpoints for fast spatial queries 
    midpoints = np.asarray(list(vecAndPos.keys()), dtype=float)
    tree = cKDTree(midpoints)
    
    order_params = []
    for cell in cells:
        # Get the vectors in the cell
        vectors = get_unit_vectors_in_cell(cell, vecAndPos, midpoints, tree)
        # if there are more than 3 vectors in the cell, calculate the order parameter
        if len(vectors) < 3:
            continue    
        # Get the unit orientation vector for the cell
        unitOrientationVector = get_unit_orientation(vectors)
        # Calculate the order parameter for this cell
        order_param = order_parameter(unitOrientationVector, vectors)
        order_params.append(order_param)
        
    # average the order parameters over all cells
    if order_params:
        return np.mean(order_params)
    else:
        return np.nan
    
