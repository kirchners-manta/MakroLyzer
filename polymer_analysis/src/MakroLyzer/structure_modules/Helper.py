import numpy as np

from MakroLyzer.errorOutputs.ErrorOutputs import ErrorOutputs


def get_Gtensor_eigVal_eigVec(graph):
    """
    Calculate the radius of gyration tensor for the graph and return its eigenvalues and eigenvectors.

    Arkin, H.; Janke, W. Ground-state properties of a polymer chain in an attractive sphere. J. Phys. Chem. B 2012, 116, 10379-10386.

              (∑_i (x_i - x_cm)^2)  (∑_i (x_i - x_cm)(y_i - y_cm))  (∑_i (x_i - x_cm)(z_i - z_cm))
    S = 1/N * (∑_i (y_i - y_cm)(x_i - x_cm))  (∑_i (y_i - y_cm)^2)  (∑_i (y_i - y_cm)(z_i - z_cm))
              (∑_i (z_i - z_cm)(x_i - x_cm))  (∑_i (z_i - z_cm)(y_i - y_cm))  (∑_i (z_i - z_cm)^2)

    Args:
        graph (GraphManager): The graph to calculate the radius of gyration tensor for.

    Returns:
        tuple: (eigenvalues, eigenvectors) of the radius of gyration tensor.
    """

    com = graph.get_com()
    S = np.zeros((3, 3))

    for node in graph.nodes():
        x, y, z = graph.get_coordinates(node) - com
        S[0, 0] += x ** 2
        S[1, 1] += y ** 2
        S[2, 2] += z ** 2
        S[0, 1] += x * y
        S[0, 2] += x * z
        S[1, 2] += y * z

    S[1, 0] = S[0, 1]
    S[2, 0] = S[0, 2]
    S[2, 1] = S[1, 2]

    S /= graph.number_of_nodes()
    eigenvalues, eigenvectors = np.linalg.eigh(S)
    
    # Conditions check
    if np.any(eigenvalues < 0):
        raise ValueError(ErrorOutputs.NEGATIVE_G_EIGENVALUE_WARNING)
    elif np.all(eigenvalues == 0):
        raise ValueError(ErrorOutputs.ZERO_G_EIGENVALUE_WARNING)
    
    return eigenvalues, eigenvectors