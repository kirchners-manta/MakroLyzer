import numpy as np

from PolyLyzer.structure_modules.radiusOfGyration import get_radius_of_gyration_tensor

def asphericity_parameter(eigenvalues):
    """
    Calculate the asphericity parameter from the eigenvalues of the radius of gyration tensor.
    Arkin, H.; Janke, W. Ground-state properties of a polymer chain in an attractive sphere. J. Phys. Chem. B 2012, 116, 10379-10386.
    
    The asphericity parameter is defined as:
        b = λ1 - 1/2 * (λ2 + λ3)
        
    where λ1 >= λ2 >= λ3 are the eigenvalues of the radius of gyration tensor.
    
    Args:
        eigenvalues (np.ndarray): The eigenvalues of the radius of gyration tensor.
    Returns:
        float: The asphericity parameter.
    """
    # Sort eigenvalues in descending order
    eigenvalues = np.sort(eigenvalues)[::-1]  
    
    # Calculate the asphericity parameter
    b = eigenvalues[0] - 0.5 * (eigenvalues[1] + eigenvalues[2])
    
    return b


def get_asphericity_parameter(graph):
    # Get the radius of gyration tensor
    eigenvalues, eigenvectors = get_radius_of_gyration_tensor(graph)
    
    # Calculate the asphericity parameter
    b = asphericity_parameter(eigenvalues)
    
    return b