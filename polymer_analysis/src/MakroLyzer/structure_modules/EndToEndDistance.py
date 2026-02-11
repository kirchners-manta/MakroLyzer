import numpy as np

from MakroLyzer.structure_modules.structureBase import StructureAnalyzer

class EndToEndDistanceAnalyzer(StructureAnalyzer):
    """
    Analyzer to calculate the end-to-end distances of all connected subgraphs in a 
    given graph.
    Inherits from StructureAnalyzer. 
    """
    
    def __init__(self, output_handler=None, static_topology=False, backbone_cache=None):
        """
        Initialize EndToEndDistanceAnalyzer.
        
        Args:
            output_handler (OutputHandler): Handler for writing output.
        """
        super().__init__(output_handler)
        self.static_topology = static_topology
        self._cached_endpoints = None
        self.backbone_cache = backbone_cache

    def _compute_endpoints(self, graph):
        """
        Compute end-to-end endpoints per subgraph on a reduced graph.
        Returns a list of (startNode, endNode).
        """
        subgraphs = graph.get_subgraphs()
        prepared = []
        for subgraph in subgraphs:
            sub = subgraph.remove_1order()
            sub.update_degree()
            prepared.append(sub)

        endpoints = []
        for subgraph in prepared:
            backbone = subgraph.find_longest_path()
            if len(backbone) < 2:
                continue
            endpoints.append((backbone[0], backbone[-1]))

        return endpoints
        
    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler:
            header = "Frame, End-to-End Distance (subgraphs) / Å"
            if self.output_handler.mode == 'streaming':
                self.output_handler.initialize_file(header)
                
    def compute(self, graph):
        """
        Calculate the end-to-end distance for all subgraphs of the given graph.

        Args:
            graph (GraphManager): Graph to analyze.
        """
        
        distances = []
        if self.backbone_cache is not None:
            endpoints = self.backbone_cache.get_endpoints(graph)
            for startNode, endNode in endpoints:
                distance = graph.distance(startNode, endNode)
                distances.append(distance)
        elif self.static_topology:
            if self._cached_endpoints is None:
                self._cached_endpoints = self._compute_endpoints(graph)

            for startNode, endNode in self._cached_endpoints:
                distance = graph.distance(startNode, endNode)
                distances.append(distance)
        else:
            subgraphs = graph.get_subgraphs()
            prepared = []
            for subgraph in subgraphs:
                sub = subgraph.remove_1order()
                sub.update_degree()
                prepared.append(sub)

            for subgraph in prepared:
                backbone = subgraph.find_longest_path()
                if len(backbone) < 2:
                    continue
                else:
                    startNode = backbone[0]
                    endNode = backbone[-1]
                    distance = subgraph.distance(startNode, endNode)
                    distances.append(distance)
                
        return np.array(distances)
                
    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame.
        
        Args:
            data (list): The computed b data.
            frame_idx (int): Current frame number.
        """
        if self.output_handler:
            row = f"{frame_idx}," + ','.join(f"{dist:.3f}" if not np.isnan(dist) else "nan" for dist in data)
        
        self.output_handler.append_row(row)
        
    def finalize_output(self):
        """
        Finalize output file (write header and rows - 'collect' mode)
        """
        header = "Frame, End-to-End Distance (subgraphs) / Å"
        super().finalize_output(header)
        
