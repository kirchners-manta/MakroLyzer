import numpy as np

from MakroLyzer.structure_modules import graphs

def count_subgraphs(graph):
    """
    Count the number of subgraphs in a graph object.
    
    Args:
        GraphManager: The graph to count subgraphs for.
    Returns:
        int: The number of subgraphs in the graph.
    """
    subgraphs = graph.get_subgraphs()
    return len(subgraphs)


def count_rings(graph):
    """
    Count the number of rings and strands in a graph.

    Returns:
        [ring_count, strand_count]
    """
    # Graph without 1-order nodes
    G_wo1 = graph.remove_1order()

    ring_count = 0
    strand_count = 0
    for subgraph in graph.get_subgraphs():
        backbone = subgraph.find_longest_path()

        # restrict the filtered graph to the subgraph's nodes
        sg_wo1 = G_wo1.subgraph(subgraph.nodes()).copy()

        # consider only backbone nodes that survived in the filtered graph
        bb_in = [n for n in backbone if sg_wo1.has_node(n)]

        # decide strand vs ring
        if any(sg_wo1.degree(n) == 1 for n in bb_in):
            strand_count += 1
        else:
            ring_count += 1

    return [ring_count, strand_count]
