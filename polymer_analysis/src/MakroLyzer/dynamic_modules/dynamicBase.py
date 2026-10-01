"""
This is the base class for all dynamic analysis modules.
Each dynamic analyzer handles computation and its own output formatting.

-> All dynamic modules inherit from this class. 
   Thus, they implement the methods defined here.
"""

# abc - Abstract Base Classes
# Infrastructure for defining abstract base classes in Python
from abc import ABC, abstractmethod
# os - Miscellaneous operating system interfaces
# Provides a way of using operating system dependent functionality
import os


class OutputHandler:
    """
    Base class for handling the output of a dynamic analyzer.
    """
    
    def __init__(self, file_path, mode='collect'):
        """
        Initialize Output Handler.
        _initialized: Flag to check if output file is already created (for 'streaming' mode)

        Args:
            file_path (str): Path where to save the output file.
            mode (str, optional): Mode of writing output. 'collect' (accumulate and write later) or 'streaming' (write immediately). Defaults to 'collect'.
        """
        self.file_path = file_path
        self.mode = mode
        self.accumulated_rows = []
        self.accumulated_matrices = {}
        self._initialized = False
        
    def write_csv(self, header, rows):
        """
        Write output to csv.
        We need this if we treat chunks of the trajecotry.

        Args:
            header (str): Header row for the CSV file.
            rows (list): List of rows to write.
        """
        # If the file exists, append a number to the filename to avoid overwriting
        base, ext = os.path.splitext(self.file_path)
        counter = 1
        while os.path.exists(self.file_path):
            self.file_path = f"{base}_{counter}{ext}"
            counter += 1
            
        # Write data
        with open(self.file_path, 'w') as f:
            f.write(header + '\n')
            for row in rows:
                f.write(row + '\n')
                
    def write_raw(self, data):
        """
        Write raw data to output file.

        Args:
            data (str): Raw data string to write.
        """
        # If the file exists, append a number to the filename to avoid overwriting
        base, ext = os.path.splitext(self.file_path)
        counter = 1
        while os.path.exists(self.file_path):
            self.file_path = f"{base}_{counter}{ext}"
            counter += 1
            
        with open(self.file_path, 'w') as f:
            f.write(data)
                
    def initialize_file(self, header):
        """
        Initialize file (for 'streaming' mode)
        """
        if not self._initialized:
            # If the file exists, append a number to the filename to avoid overwriting
            base, ext = os.path.splitext(self.file_path)
            counter = 1
            while os.path.exists(self.file_path):
                self.file_path = f"{base}_{counter}{ext}"
                counter += 1
            
            with open(self.file_path, 'w') as f:
                f.write(header + '\n')
            self._initialized = True
            
    def append_row(self, row):
        """
        Add row of current frame to outputfile ('streaming' mode)
        Or add row to the accumulated rows ('collect' mode)

        Args:
            row (list): Result row to append to the file.
        """
        if self.mode == 'streaming':
            with open(self.file_path, 'a') as f:
                f.write(row + '\n')
        else:
            self.accumulated_rows.append(row)
            
    def append_matrix(self, matrix, frame_idx):
        """
        Create an outputfile for the current frame ('streaming' mode)
        Or add matrix to the accumulated matrices ('collect' mode)

        Args:
            matrix (list of lists): Matrix to append.
            frame_idx (int): Current frame number.
        """
        if self.mode == 'streaming':
            # If the file exists, append a number to the filename to avoid overwriting
            base, ext = os.path.splitext(self.file_path)
            counter = 1
            file_path = f"{base}_frame{frame_idx}{ext}"
            while os.path.exists(file_path):
                file_path = f"{base}_frame{frame_idx}_{counter}{ext}"
                counter += 1
                
            with open(file_path, 'w') as f:
                for row in matrix:
                    f.write(','.join(map(str, row)) + '\n')
        else:
            self.accumulated_matrices[frame_idx] = matrix
            
    def finalize(self, header=None):
        """
        Write accumulated rows to file ('collect' mode)

        Args:
            header (str, optional): Header of the file. Defaults to None.
        """
        if self.mode == 'collect' and self.accumulated_rows:
            # If the file exists, append a number to the filename to avoid overwriting
            base, ext = os.path.splitext(self.file_path)
            counter = 1
            while os.path.exists(self.file_path):
                self.file_path = f"{base}_{counter}{ext}"
                counter += 1
                
            # Write data
            with open(self.file_path, 'w') as f:
                if header:
                    f.write(header + '\n')
                for row in self.accumulated_rows:
                    f.write(row + '\n')
                    
    def finalize_matrices(self):
        """
        Write accumulated matrices to files ('collect' mode)
        Each frame gets its own file named with the frame index.
        """
        if self.mode == 'collect' and self.accumulated_matrices:
            for frame_idx, matrix in self.accumulated_matrices.items():
                # If the file exists, append a number to the filename to avoid overwriting
                base, ext = os.path.splitext(self.file_path)
                counter = 1
                file_path = f"{base}_frame{frame_idx}{ext}"
                while os.path.exists(file_path):
                    file_path = f"{base}_frame{frame_idx}_{counter}{ext}"
                    counter += 1
                    
                with open(file_path, 'w') as f:
                    for row in matrix:
                        f.write(','.join(map(str, row)) + '\n')
                    
class DynamicAnalyzer(ABC):
    """
    Base interface for per-frame and accumulated dynamic analyses.

    Main calls initialize_output(), run() for each analyzed frame, then
    finalize() and finalize_output(), in that order. Concrete analyzers own
    their timing parameters, reference coordinates, and history buffers.
    """

    def __init__(self, output_handler=None):
        """
        Initialize dynamic analyzer.

        Args:
            output_handler (OutputHandler): Handler for writing output. 
        """
        self.output_handler = output_handler
        self.frame_number = 0  # Number of analyzed frames processed successfully.

    @abstractmethod
    def compute(self, graph):
        """
        Perform some crazy calculations to obtain the results. 
        (Optionally return a result)

        Args:
            graph : The molecular graph to analyze.
            
        Returns; 
            Analysis results (format depends on the analyzer.)
            Per-frame analyzers such as RMSD return the current value. Accumulating
            analyzers such as MSD may update internal state and return None.
            Copy retained coordinates because main may update the graph in place.
        """
        pass

    @abstractmethod
    def render_output(self, data, frame_idx):
        """
        Format and write the analysis output for a single frame.

        Args:
            data : Results from the compute function.
            frame_idx (int): Current frame number.
            
        Returns:
            Formatted output.
        """
        pass

    def initialize_output(self):
        """
        Initialize output file with header.
        Called once before first frame.
        """
        pass

    def finalize_output(self, header=None):
        """
        Finalize output after all frames are processed ('collect' mode)
        """
        if self.output_handler and self.output_handler.mode == 'collect':
            self.output_handler.finalize(header)
            
    def run(self, graph, frame_idx):
        """
        Execute analysis for a single frame: compute + output

        Args:
            graph : The molecular graph to analyze
            frame_idx (int): Current frame number
            
        Returns: 
            Analysis Result.
        """
        result = self.compute(graph)
        if self.output_handler is not None and result is not None:
            self.render_output(result, frame_idx)
        self.frame_number += 1
        return result

    def finalize(self):
        """
        Finish accumulated calculations; per-frame analyzers need no action.
        """
        return None
