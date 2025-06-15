import numpy as np
import csv

from PolyLyzer.structure_modules import graphs


    
def get_dihedrals(graph, sign=None):
    """
    Calculate the dihedral angles of the graph.
    The dihedral angles are calculated for each subgraph of the graph.
    The default for the sign is None, which means the dihedrals are in the range of 0 to 180 degrees.

    Args:
        GraphManager: The graph to calculate the dihedral angles for.
    """
    
    # Remove 1-order nodes, find subgraphs and surrounding atoms
    GraphWithout1order = graph.remove_1order()
    GraphWithout1order.surrounding()
    GraphWithout1order.update_degree()
    # Get the subgraphs of the graph
    subgraphs = GraphWithout1order.get_subgraphs()
    dihedrals = []
    # For each subgraph, find the longest path
    for subgraph in subgraphs:
        longestPath = subgraph.find_longest_path()
        # For each node in the longest path, find the dihedral angles
        for i in range(len(longestPath) - 3):
            node1 = longestPath[i]
            node2 = longestPath[i + 1]
            node3 = longestPath[i + 2]
            node4 = longestPath[i + 3]
            d = subgraph.dihedral(node1, node2, node3, node4, sign=sign)
            # round dihedral to integers
            d = round(d)
            dihedrals.append(d)
            
    return dihedrals
            



def get_all_dihedrals(graph, sign=None):
    """
    Get all dihedrals of the graph and write them to a .csv file.
    The dihedrals are sorted, counted and written to a .csv file with two columns: Dihedral and Count.
    They are in the range of -180 to 180 degrees if sign is True, otherwise in the range of 0 to 180 degrees.

    Args:
        graph (GraphManager): The graph to calculate the dihedral angles for.
        file (str): The name of the .csv file to write the dihedrals to.
        sign (bool, optional): If True, dihedrals are in the range of -180 to 180 degrees.
        
    Returns:
        list: A list of tuples containing the dihedral angles and their counts.
    """
    
    dihedrals = get_dihedrals(graph, sign=sign)
    
    # round dihedrals to integers
    dihedrals = [round(d) for d in dihedrals]

    # Group dihedrals 
    Dihedrals = dict(sorted({x: dihedrals.count(x) for x in set(dihedrals)}.items(), key=lambda item: item[0]))
    # Add 0 count for missing dihedrals between -180/0 and 180
    if sign is not None:
        for i in range(-180, 181):
            if i not in Dihedrals:
                Dihedrals[i] = 0
    elif sign is None:
        for i in range(0, 181):
            if i not in Dihedrals:
                Dihedrals[i] = 0
    # Convert the counts to a list of tuples
    dihedrals = [(k, v) for k, v in Dihedrals.items()]
    # Sort the dihedrals by size
    dihedrals = sorted(dihedrals, key=lambda x: x[0])
            
    return dihedrals

def get_CisTrans(graph):
    """
    Get the cis and trans counts of the graph and write them to a .csv file.
    Args:
        graph (GraphManager): The graph to calculate the cis and trans counts for.
        file (str): The name of the .csv file to write the cis and trans counts to.
        
    Returns:
        list: A list containing the counts of cis and trans.
    """
    
    dihedrals = get_dihedrals(graph, sign=None)
    
    # 0 to 90 degrees are cis, 90 to 180 degrees are trans
    cis = sum(1 for d in dihedrals if 0 <= d <=90)
    trans = sum(1 for d in dihedrals if 90 < d <=180)
    
    cisTrans = [('Cis', cis), ('Trans', trans)]
    
    return cisTrans