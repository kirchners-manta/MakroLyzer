import numpy as np

from PolyLyzer.structure_modules import graphs

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