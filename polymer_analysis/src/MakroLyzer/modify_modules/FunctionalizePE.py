from MakroLyzer.modify_modules.structureModifierBase import StructureModifier

class FunctionalizePEModifier(StructureModifier):
    """
    A class for functionalizing polyethylene (PE) features. 
    1. Replacing CH2 by CO
    2. Replacing CH2 by NH
    A number between 0 and 100 indicating the percentage of functionalization.
    Inherits fron StructureModifier.
    """
    
    def __init__(self, percentage, func_type, output_handler=None):
        """
        Initializes the FunctionalizePE class with the given percentage and functionalization type.

        Args:
            percentage (int): The percentage of functionalization (0-100).
            func_type (str): The type of functionalization ('CO' or 'NH').
        """
        super().__init__(output_handler)
        if func_type not in ['CO', 'NH']:
            raise ValueError("func_type must be 'CO'")
        if not (0 <= percentage <= 100):
            raise ValueError("percentage must be between 0 and 100")
        
        self.percentage = percentage
        self.func_type = func_type
        
    def compute(self, graph):
        """
        Functionalizes the polyethylene features in the given graph.
        """
        graph.functionalizePE(self.percentage, self.func_type)
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
