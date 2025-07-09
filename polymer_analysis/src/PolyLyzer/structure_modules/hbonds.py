import numpy as np
import csv

from PolyLyzer.structure_modules import graphs

def get_Hbonds(graph, cutoffs):
    """
    Calculate the hydrogen bonds in the graph.
    The hydrogen bonds are calculated for each subgraph of the graph.
    
    Args:
        graph (GraphManager): The graph to calculate the hydrogen bonds for.
        elementTypes (list): A list of element types to consider for hydrogen bonds.
        distance (list, optional): A list of distances for each pair of element types.
        
    Returns:
        list: A list of tuples containing numbers of hydrogen bonds and their corresponding atom pairs.
    """
    
    numberOfHbonds = []
    
    for elementType, HAcceptor_dist, DonorAcceptor_dist, Angle_cut in cutoffs:        
        # Get the hydrogen bonds for the element type
        hbonds = graph.get_hbonds(elementType, HAcceptor_dist, DonorAcceptor_dist, Angle_cut)
        numberOfHbonds.append(len(hbonds))
        
    hbonds = [(elementType, HAcceptor_dist, DonorAcceptor_dist, Angle_cut, num) for (elementType, HAcceptor_dist, DonorAcceptor_dist, Angle_cut), num in zip(cutoffs, numberOfHbonds)]
    return hbonds