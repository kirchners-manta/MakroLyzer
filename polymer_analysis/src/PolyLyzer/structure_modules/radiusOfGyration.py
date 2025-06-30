import numpy as np
import time

from PolyLyzer.structure_modules import graphs

def get_radius_of_gyration(graph):
    """
    Calculate the radius of gyration of each subgraph in the graph and of the entire graph.
    The radius of gyration is calculated as the square root of the average of the squared distances of each atom from the center of mass.
    
    R²g = (1/N) * Σ (r_i - r_cm)²

    Args:
        graph (GraphManager): The graph to calculate the radius of gyration for.
    Returns:
        tuple: A tuple containing the radius of gyration for each subgraph and the entire graph.
    """
    # Radius of gyration for the entire graph
    _, coords = graph.get_all_coordinates()
    com = graph.get_com()
    squared_distances = np.sum((coords - com)**2, axis=1)
    R_whole = np.sqrt(np.mean(squared_distances)) # 1/N * squared_distances
    
    # Radius of gyration for each subgraph
    r_subgraphs = []
    for subgraph in graph.get_subgraphs():
        _, coords_sub = subgraph.get_all_coordinates()
        com_sub = subgraph.get_com()
        squared_distances_sub = np.sum((coords_sub - com_sub)**2, axis=1)
        R_sub = np.sqrt(np.mean(squared_distances_sub))
        r_subgraphs.append(R_sub)
    
                
    return r_subgraphs, R_whole


def get_radius_of_gyration_tensor(graph):
    """"
        Calculate the radius of gyration tensor for the graph.
        Arkin, H.; Janke, W. Ground-state properties of a polymer chain in an attractive sphere. J. Phys. Chem. B 2012, 116, 10379-10386.
        
                  (∑_i (x_i - x_cm)²)  (∑_i (x_i - x_cm)(y_i - y_cm))  (∑_i (x_i - x_cm)(z_i - z_cm))
        S = 1/N * (∑_i (y_i - y_cm)(x_i - x_cm))  (∑_i (y_i - y_cm)²)  (∑_i (y_i - y_cm)(z_i - z_cm))
                  (∑_i (z_i - z_cm)(x_i - x_cm))  (∑_i (z_i - z_cm)(y_i - y_cm))  (∑_i (z_i - z_cm)²)
        
        where: 
            N is the number of atoms in the graph,
            (x_cm, y_cm, z_cm) is the center of mass of the graph,
            (x_i, y_i, z_i) are the coordinates of atom i.
            -> S is a symmetric 3x3 matrix.
            Unit of S is Å².

    Args:
        graph (GraphManager): The graph to calculate the radius of gyration tensor for.
    Returns:
        Eigenvalues, Eigenvectors: The eigenvalues and eigenvectors of the radius of gyration tensor.
    """
    
    # Get the com
    com = graph.get_com()
    
    # Initialize the tensor
    S = np.zeros((3, 3))
    
    # Calculate the tensor
    for node in graph.nodes():
        coords = graph.get_coordinates(node)
        x, y, z = coords - com
        
        S[0, 0] += x**2
        S[1, 1] += y**2
        S[2, 2] += z**2
        S[0, 1] += x * y
        S[0, 2] += x * z
        S[1, 2] += y * z
        
    S[1,0] = S[0, 1]
    S[2,0] = S[0, 2]
    S[2,1] = S[1, 2]
        
    # Normalize and Diagonalize the tensor
    S /= graph.number_of_nodes()
    eigenvalues, eigenvectors = np.linalg.eigh(S)
    
    return eigenvalues, eigenvectors


    
    