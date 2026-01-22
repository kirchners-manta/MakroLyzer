from MakroLyzer.modify_modules.structureModifierBase import StructureModifier

class FunctionalizePEModifier(StructureModifier):
    """
    A class for functionalizing polyethylene (PE) features. 
    1. Replacing CH2 by CO
    2. Replacing C2H4 by COO 
    3. Replacing C2H4 by CONH
    A number between 0 and 100 indicating the percentage of functionalization.
    Inherits fron StructureModifier.
    """
    
    def __init__(self, mode, func_type, percentage=None, neighbor_exclusion=None, distance=None, output_handler=None):
        """
        Initializes the FunctionalizePE class with the given percentage and functionalization type.

        Args:
            percentage (float): The percentage of functionalization (0-100).
            func_type (str): The type of functionalization ('CO', 'COH', 'ester' or 'amide').
            neighbor_exclusion (int): Minimum C-C hop distance to avoid selecting adjacent nodes.
        """
        super().__init__(output_handler)
        if func_type not in ['CO', 'COH', 'ester', 'amide']:
            raise ValueError("func_type must be 'CO', 'COH', 'ester' or 'amide'")
        if mode not in ['random', 'periodic']:
            raise ValueError("mode must be 'random' or 'periodic'.")
        if mode == 'random':
            if percentage is None:
                raise ValueError("percentage is required for random mode.")
            if not (0 <= percentage <= 100):
                raise ValueError("percentage must be between 0 and 100")
        else:
            if distance is None:
                raise ValueError("distance is required for periodic mode.")
        
        self.mode = mode
        self.func_type = func_type
        self.percentage = percentage
        self.neighbor_exclusion = neighbor_exclusion
        self.distance = distance
        
    def compute(self, graph):
        """
        Functionalizes the polyethylene features in the given graph.
        """
        graph.functionalizePE(
            self.mode,
            self.func_type,
            percentage=self.percentage,
            neighbor_exclusion=self.neighbor_exclusion,
            distance=self.distance,
        )
        return graph

    def render_output(self, graph, frame_idx):
        """
        Write graph xyz to file.

        Args:
            graph (GraphManager): The modified graph.
            frame_idx (int): Current frame number.
        """
        NoNodes = graph.number_of_nodes()
        header = f"{NoNodes}\nelement  x         y        z"
        self.output_handler.write_output(header, graph) 
