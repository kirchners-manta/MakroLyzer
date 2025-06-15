import numpy as np
import csv

from PolyLyzer.structure_modules import graphs

def get_Hbonds(graph, TypesDistances):
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
    
    for elementType, distance in TypesDistances:        
        # Get the hydrogen bonds for the element type
        hbonds = graph.get_hbonds(elementType, distance)
        numberOfHbonds.append(len(hbonds))
        
    hbonds = [(elementType, distance, num) for (elementType, distance), num in zip(TypesDistances, numberOfHbonds)]
    # Write the hydrogen bonds to a CSV file
    with open('hydrogen_bonds.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Element Type', 'Distance', 'Number of Hydrogen Bonds'])
        for element, distance, count in hbonds:
            writer.writerow([element, distance, count]) 
        
    return hbonds