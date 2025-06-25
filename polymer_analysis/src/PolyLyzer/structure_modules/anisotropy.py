import numpy as np

from PolyLyzer.structure_modules.radiusOfGyration import get_radius_of_gyration_tensor

def anisotropy_factor(eigenvalues):
    """
    Calculate the anisotropy factor from the eigenvalues of the radius of gyration tensor.
    
    The anisotropy factor is defined as:
        K = 1 - 3*(λ1λ2 + λ2λ3 + λ3λ1) / (λ1 + λ2 + λ3)² 
        
    where λ1 >= λ2 >= λ3 are the eigenvalues of the radius of gyration tensor.
    
    Args:
        eigenvalues (np.ndarray): The eigenvalues of the radius of gyration tensor.
    Returns:
        float: The anisotropy factor.
    """
    # Sort eigenvalues in descending order
    eigenvalues = np.sort(eigenvalues)[::-1]
    
    # Calculate the anisotropy factor
    K = 1 - 3 * (eigenvalues[0] * eigenvalues[1] + eigenvalues[1] * eigenvalues[2] + eigenvalues[2] * eigenvalues[0]) / (eigenvalues.sum() ** 2)
    
    return K

def get_anisotropy_factor(graph):
    # Get the radius of gyration tensor
    eigenvalues, eigenvectors = get_radius_of_gyration_tensor(graph)
    
    # Calculate the anisotropy factor
    K = anisotropy_factor(eigenvalues)
    
    return K