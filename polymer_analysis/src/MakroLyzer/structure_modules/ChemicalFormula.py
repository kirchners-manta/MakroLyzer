from MakroLyzer.structure_modules.structureBase import StructureAnalyzer

class ChemicalFormulaAnalyzer(StructureAnalyzer):
    """
    Analyzer for the chemical formulas of the molecules/subgraphs in the given molecular graph.
    Inherits from StructureAnalyzer.
    """
    
    def __init__(self, output_handler=None):
        """
        Initialize the ChemicalFormulasAnalyzer.

        Args:
            output_handler (OutputHandler): Handler for writing output.
        """
        super().__init__(output_handler)
        
    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler:
            header = "Frame, Chemical Formula, Count"
            if self.output_handler.mode == 'streaming':
                self.output_handler.initialize_file(header)
                
    def compute(self, graph):
        """
        Get the Chemical Formulas and their counts for the given graph.

        Args:
            graph (GraphManager): Graph to analyze.
            
        Returns:
            dict: Chemical Formula, Count
        """
        # Get the subgraphs of the graph
        subgraphs = graph.get_subgraphs()
        formulas = {}
        
        # Iterate through the subgraphs and count the elements
        for subgraph in subgraphs:
            formula = {}
            for node in subgraph.nodes:
                element = subgraph.nodes[node]['element']
                if element not in formula:
                    formula[element] = 0
                formula[element] += 1
            
            # sort alphabetically
            formula = dict(sorted(formula.items()))
            # Convert the formula to a string
            formula_str = ''.join([f"{k}{v}" for k, v in formula.items()])
            formulas[formula_str] = formulas.get(formula_str, 0) + 1
            
        # Sort the formulas by their counts
        formulas = dict(sorted(formulas.items(), key=lambda item: item[1], reverse=True))   
        # Convert the counts to a list of tuples
        formulas = [(k, v) for k, v in formulas.items()]
        
        return formulas
    
    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame.
        
        Args:
            data (dict): Chemical Formulas and their counts for the given graph.
            frame_idx (int): Current frame number.
        """
        if self.output_handler:
            for ChemFormula, count in data:
                row = f"{frame_idx},{ChemFormula},{count}"
                self.output_handler.append_row(row)
                
    def finalize_output(self, header=None):
        """
        Finalize output file (write header and rows - 'collect' mode)
        """
        header = "Frame, Chemical Formula, Count"
        super().finalize_output(header)