import numpy as np

from PolyLyzer.structure_modules import graphs


def end_to_end_dist(graph):
    """
    Calculate the end-to-end distance of a graph object.
    
    Args:
        GraphManager: The graph to calculate the distance for.
    Returns:
        list of floats: The end-to-end distances of all subgraphs.
    """
    subgraphs = graph.get_subgraphs()
    prepared = []
    for subgraph in subgraphs:
        sub = subgraph.remove_1order()
        sub.update_degree()
        prepared.append(sub)
        
    distances = []
    for subgraph in prepared:
        longestPath = subgraph.find_longest_path()
        
        if len(longestPath) < 2:
            # Subgraph is probably a ring
            distances.append(np.nan)
            continue

        startNode = longestPath[0]
        endNode = longestPath[-1]
        distance = subgraph.distance(startNode, endNode)
        distances.append(distance)
    
    return np.array(distances)


def end_to_end_dist_ensemble(graph):
    """
    Calculate the ensemble average of the end-to-end distance (root of ⟨R²⟩) 
    for all subgraphs representing polymer chains.
    
    Args:
        GraphManager: The graph to calculate the distance for.
    Returns:
        float: Ensemble average of the end-to-end distance.
    """
    subgraphs = graph.get_subgraphs()
    subgraphsWithout1order = []
    for subgraph in subgraphs:
        sub = subgraph.remove_1order()
        sub.update_degree()
        subgraphsWithout1order.append(sub)
    subgraphScalarProduct = []
    for subgraph in subgraphsWithout1order:
        longestPath = subgraph.find_longest_path()
        if len(longestPath) < 2:
            continue
        vectors = []
        for i in range(len(longestPath) - 1):
            node1 = longestPath[i]
            node2 = longestPath[i + 1]
            vector = subgraph.vector(node1, node2)
            vectors.append(vector)
        vectors = np.array(vectors)
        # Calculate the scalar product of the vectors
        R2 = np.sum(vectors @ vectors.T)  
        subgraphScalarProduct.append(R2)
    if not subgraphScalarProduct:
        return 0.0  
    ensembleAverageR2 = np.mean(subgraphScalarProduct)
    return np.sqrt(ensembleAverageR2), np.std(subgraphScalarProduct)