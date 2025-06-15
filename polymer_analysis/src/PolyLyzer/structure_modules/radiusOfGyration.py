import numpy as np
import csv

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
    R_whole = 1/graph.number_of_nodes() * np.sum(
        [(graph.get_coordinates(node) - graph.get_com())**2 for node in graph.nodes()]
    )
    R_whole = np.sqrt(R_whole)
    
    # Radius of gyration for each subgraph
    subgraphs = graph.get_subgraphs()
    radius_of_gyration_subgraphs = []
    
    for subgraph in subgraphs:
        R_sub = 1/subgraph.number_of_nodes() * np.sum(
            [(subgraph.get_coordinates(node) - subgraph.get_com())**2 for node in subgraph.nodes()]
        )
        R_sub = np.sqrt(R_sub)
        radius_of_gyration_subgraphs.append(R_sub)
        
    return radius_of_gyration_subgraphs, R_whole
    
    
    