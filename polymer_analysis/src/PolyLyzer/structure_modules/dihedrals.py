import numpy as np
import csv

from PolyLyzer.structure_modules import graphs


    
def get_dihedrals(graph, file, sign=None):
    """
    Calculate the dihedral angles of the graph.
    The dihedral angles are calculated for each subgraph of the graph.
    The dihedral angles are sorted and counted.
    The dihedral angles are written to a .csv file.

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
    # Print dihedrals to .cvs file
    with open(file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        writer.writerow(['Dihedral', 'Count'])
        for d in dihedrals:
            writer.writerow([d[0], d[1]])
            
    return dihedrals