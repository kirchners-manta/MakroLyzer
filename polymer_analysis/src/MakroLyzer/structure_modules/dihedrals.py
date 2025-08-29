import numpy as np
import csv

from MakroLyzer.structure_modules import graphs


    
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
        sub_list =  []
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
            sub_list.append(d)
        dihedrals.append(sub_list)

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
    
    dihedral_list = get_dihedrals(graph, sign=sign)

    # round dihedrals to integers
    dihedrals = [round(d) for sublist in dihedral_list for d in sublist]

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
            
    return dihedrals, dihedral_list

def get_CisTrans(graph):
    """
    Get the cis and trans counts of the graph and write them to a .csv file.
    Args:
        graph (GraphManager): The graph to calculate the cis and trans counts for.
        file (str): The name of the .csv file to write the cis and trans counts to.
        
    Returns:
        list: A list containing the counts of cis and trans.
    """

    dihedral_list = get_dihedrals(graph, sign=None)

    # 0 to 90 degrees are cis, 90 to 180 degrees are trans
    cis = sum(1 for sublist in dihedral_list for d in sublist if 0 <= d <=90)
    trans = sum(1 for sublist in dihedral_list for d in sublist if 90 < d <=180)

    cisTrans = [('Cis', cis), ('Trans', trans)]
    
    return cisTrans

def get_Ramachandran(graph):
    """
    Get the Ramachandran plot data of the graph.
    The Polymer needs to be an Amino Acid or an alternative Amino Acid. 

    Args:
        graph (GraphManager): The graph to calculate the Ramachandran plot data for.

    Returns:
        list: A list containing the phi and psi angles for the Ramachandran plot.
    """  
    
    # Prepare a ramachandran matrix of ints and initialize to 0
    ramachandran = [[0 for _ in range(360)] for _ in range(360)]
    
    # Get subgraphs
    subgraphs = graph.get_subgraphs()

    for subgraph in subgraphs:

        # Get the backbone of the AA which is always like this: 
        # N-C-C-N-...
        # N-Calpha-CarbonylC-N-...
        backbone = subgraph.AminoAcidBackbone()

        # Slide over the backbone and calculate phi/psi angles
        # Phi is defined as the dihedral angle C(i)-N(i+1)-Ca(i+1)-C(i+1)
        # Psi is defined as the dihedral angle N(i)-Ca(i)-C(i+1)-N(i+1)

        # We start with the first dihedral at the third atom in the backbone which is the
        # first CarbonylC = Cprime
        i = 2
        while i < len(backbone) - 4:
            Cprime1 = backbone[i]
            N1 = backbone[i+1]
            Calpha = backbone[i+2]
            Cprime2 = backbone[i+3]
            N2 = backbone[i+4]

            # phi C(i)-N(i+1)-Ca(i+1)-C(i+1)
            phi = graph.dihedral(Cprime1, N1, Calpha, Cprime2, sign=True)
            phi = round(phi) + 180
            # psi N(i)-Ca(i)-C(i+1)-N(i+1)
            psi = graph.dihedral(N1, Calpha, Cprime2, N2, sign=True)
            psi = round(psi) + 180

            # Increment the ramachandran matrix
            ramachandran[phi][psi] += 1

            # move to the next Calpha
            i += 3 

    return ramachandran