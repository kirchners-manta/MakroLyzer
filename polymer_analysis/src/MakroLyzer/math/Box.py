import numpy as np
from itertools import product

class Box:
    
    @staticmethod
    def devideBox(boxSize, n):
        """
        Divide the box into n^3 equal parts.
        
        Args:
            boxSize : tuple of floats
                The size of the box in each dimension (x, y, z).

            n : tuple of ints
                The number of divisions along each dimension (nx, ny, nz).
        Returns:
            list
                A list of tuples representing the start and end coordinates of each sub-box
                in the format [(x_start, x_end), (y_start, y_end), (z_start, z_end)].
        """
        
        # If only one dimension is provided, assume it is a cubic box
        if isinstance(boxSize, (int, float)):
            boxSize = (boxSize, boxSize, boxSize)
            
        if isinstance(n, (int, float)):
            n = (n, n, n)
            
        # Calculate the size of each sub-box
        x_size = boxSize[0] / n[0]
        y_size = boxSize[1] / n[1]
        z_size = boxSize[2] / n[2]
        
        # precompute the start and end coordinates for each dimension
        x_dots = [(i * x_size, (i+1) * x_size) for i in range(n[0])]
        y_dots = [(i * y_size, (i+1) * y_size) for i in range(n[1])]
        z_dots = [(i * z_size, (i+1) * z_size) for i in range(n[2])]
        
        # Create the sub-boxes by combining the start and end coordinates
        sub_boxes = list(product(x_dots, y_dots, z_dots))
        
        return sub_boxes
        
        
        