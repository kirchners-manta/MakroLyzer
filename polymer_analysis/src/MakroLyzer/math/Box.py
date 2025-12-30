import numpy as np
from itertools import product

class Box:
    
    @staticmethod
    def divideBox(BoxSize, NoCellsPerDim):
        """
        Divide the box into NoCellsPerDim_1 * NoCellsPerDim_2 * NoCellsPerDim_3 equal parts.
        
        Args:
            boxSize : tuple of floats
                The size of the box in each dimension (x, y, z).

            NoCellsPerDim : tuple of ints
                            The number of cells along each dimension.
        Returns:
            list
                A list of tuples representing the start and end coordinates of each sub-box
                in the format [(x_start, x_end), (y_start, y_end), (z_start, z_end)].
        """
        
        # If only one dimension is provided, assume it is a cubic box
        if isinstance(BoxSize, (int, float)):
            BoxSize = (BoxSize, BoxSize, BoxSize)
            
        if isinstance(NoCellsPerDim, (int, float)):
            NoCellsPerDim = (NoCellsPerDim, NoCellsPerDim, NoCellsPerDim)
            
        # Calculate the size of each sub-box
        x_size = BoxSize[0] / NoCellsPerDim[0]
        y_size = BoxSize[1] / NoCellsPerDim[1]
        z_size = BoxSize[2] / NoCellsPerDim[2]
        
        # precompute the start and end coordinates for each dimension
        x_dots = [(i * x_size, (i+1) * x_size) for i in range(NoCellsPerDim[0])]
        y_dots = [(i * y_size, (i+1) * y_size) for i in range(NoCellsPerDim[1])]
        z_dots = [(i * z_size, (i+1) * z_size) for i in range(NoCellsPerDim[2])]
        
        # Create the sub-boxes by combining the start and end coordinates
        sub_boxes = list(product(x_dots, y_dots, z_dots))
        
        return sub_boxes
        
        
        