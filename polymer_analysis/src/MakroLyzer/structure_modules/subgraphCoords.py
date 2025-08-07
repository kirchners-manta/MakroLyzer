import numpy as np

from MakroLyzer.structure_modules import graphs

def get_subgraph_coords(graph):
    """
    Get the coordinates of a subgraph within a graph.
    
    Args:
        graph (GraphManager): The graph containing the subgraph.
        
    Returns:
        numpy array: An array of coordinates for each node in the subgraph.
        The coordinates are in the format [(element_name, (x, y, z)), ...
    """
    
    subgraphs = graph.get_subgraphs()
    
    subgraph_coords = []
    for subgraph in subgraphs:
        coords_list = []
        for node in subgraph.nodes:
            name = subgraph.nodes[node]['element']
            coords = subgraph.get_coordinates(node)
            coords_list.append((name, *coords))
        subgraph_coords.append(coords_list)            
            
    return subgraph_coords
            
        