import pandas as pd
import numpy as np
from collections import Counter

from MakroLyzer.structure_modules.structureBase import StructureAnalyzer, OutputHandler
from MakroLyzer.structure_modules import dihedrals as dihedral_funcs


class DihedralsAnalyzer(StructureAnalyzer):
    """
    Analyzer for the calculation of dihedral angles and cis-trans classification.
    Inherits from StructureAnalyzer.
    
    This analyzer computes dihedrals once and uses them to generate three output files:
    1. Dihedral angle counts (dihedrals.csv) - only in 'collect' mode
    2. Dihedral angles per subgraph (dihedrals_list.csv) 
    3. Cis-trans counts (cisTrans.csv)
    
    By caching the dihedral calculations, the cis-trans classification reuses the same 
    dihedral data without recalculating.
    """
    
    def __init__(self, dihedral_output_handler=None, dihedral_list_output_handler=None, 
                 cistrans_output_handler=None, dihedral_range='abs'):
        """
        Initialize the DihedralsAnalyzer.

        Args:
            dihedral_output_handler (OutputHandler): Handler for dihedral counts output.
            dihedral_list_output_handler (OutputHandler): Handler for dihedral list output.
            cistrans_output_handler (OutputHandler): Handler for cis-trans counts output.
            dihedral_range (str): 'abs' for 0-180 degrees, 'nonabs' for -180-180 degrees.
            -> sign parameter is True for 'nonabs', None for 'abs'.
        """
        # Initialize with the primary output handler (dihedral counts)
        super().__init__(dihedral_output_handler)
        
        # Store all three output handlers
        self.dihedral_handler = dihedral_output_handler
        self.dihedral_list_handler = dihedral_list_output_handler
        self.cistrans_handler = cistrans_output_handler
        
        # Set sign parameter based on dihedral_range
        self.sign = True if dihedral_range == 'nonabs' else None
        
        # Cache for compiled dihedral counts across frames
        # Necessary since OutputHandler can only append rows, not columns
        self.all_dihedral_data = []
        
    def initialize_output(self):
        """
        Initialize output files with headers. ('streaming' mode)
        """
        if self.dihedral_list_handler and self.dihedral_list_handler.mode == 'streaming':
            header = "Frame,Subgraph,Dihedral Angles / °"
            self.dihedral_list_handler.initialize_file(header)
            
        if self.cistrans_handler and self.cistrans_handler.mode == 'streaming':
            header = "Frame,Cis count,Trans count"
            self.cistrans_handler.initialize_file(header)
    
    def compute(self, graph):
        """
        Calculate dihedral angles and cis-trans counts for the given graph.
        Computes both in one call to avoid redundant calculations.

        Args:
            graph (GraphManager): Graph to analyze.
            
        Returns:
            dict: Dictionary containing:
                - 'dihedrals': List of (angle, count) tuples
                - 'dihedral_list': List of dihedral angles per subgraph
                - 'cistrans': List of (label, count) tuples [('Cis', count), ('Trans', count)]
        """
        # Get all dihedrals and their list representation
        dihedrals, dihedral_list = self._get_all_dihedrals(graph)
        
        # Calculate cis-trans from the dihedral_list (reusing the same calculation)
        cistrans = self._get_cistrans_from_dihedrals(dihedral_list)
        
        return {
            'dihedrals': dihedrals,
            'dihedral_list': dihedral_list,
            'cistrans': cistrans
        }
        
    def _get_all_dihedrals(self, graph):
        """
        Get all dihedrals and counts of the graph in a single pass.
        Combines the logic of calculating dihedral angles and counting them efficiently.
        The dihedrals are in the range of -180 to 180 degrees if sign is True, 
        otherwise in the range of 0 to 180 degrees.

        Args:
            graph (GraphManager): The graph to calculate the dihedral angles for.

        Returns:
            tuple: (dihedrals_list, dihedral_counts)
                - dihedrals_list: List of dihedral angles per subgraph
                - dihedral_counts: List of (angle, count) tuples sorted by angle
        """
        
        # Remove 1-order nodes, find subgraphs and surrounding atoms
        GraphWithout1order = graph.remove_1order()
        GraphWithout1order.surrounding()
        GraphWithout1order.update_degree()
        
        # Get the subgraphs of the graph
        subgraphs = GraphWithout1order.get_subgraphs()
        dihedral_list = []
        angle_counter = Counter()
        
        # For each subgraph, calculate dihedrals and count them in one pass
        for subgraph in subgraphs:
            sub_list = []
            longestPath = subgraph.find_longest_path()
            
            # For each node in the longest path, calculate dihedral angles
            for i in range(len(longestPath) - 3):
                node1 = longestPath[i]
                node2 = longestPath[i + 1]
                node3 = longestPath[i + 2]
                node4 = longestPath[i + 3]
                
                d = round(subgraph.dihedral(node1, node2, node3, node4, self.sign))
                # 180 and -180 are equivalent, count them together
                if d == -180:
                    d = 180
                    
                sub_list.append(d)
                angle_counter[d] += 1
            
            dihedral_list.append(sub_list)
        
        # Create result with all angles in range filled with 0 counts
        if self.sign is not None:
            # For signed angles: -179 to 180,
            angle_range = range(-179, 181)
        else:
            # For unsigned angles: 0 to 180
            angle_range = range(0, 181)
            
        dihedral_counts = [(angle, angle_counter.get(angle, 0)) for angle in angle_range]
            
        return dihedral_counts, dihedral_list        
    
    def _get_cistrans_from_dihedrals(self, dihedral_list):
        """
        Calculate cis and trans counts from dihedral list.
        This avoids recalculating dihedrals.
        
        Args:
            dihedral_list: List of dihedral angles per subgraph
            
        Returns:
            list: [('Cis', cis_count), ('Trans', trans_count)]
        """
        # 0 to 90 degrees are cis, 90 to 180 degrees are trans
        cis = sum(1 for sublist in dihedral_list for d in sublist if 0 <= np.abs(d) <= 90)
        trans = sum(1 for sublist in dihedral_list for d in sublist if 90 < np.abs(d) <= 180)
        
        return [('Cis', cis), ('Trans', trans)]
    
    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame to all three output files.
        
        Args:
            data (dict): Dictionary with 'dihedrals', 'dihedral_list', and 'cistrans' keys.
            frame_idx (int): Current frame number.
        """
        # Cache for collect mode
        if self.dihedral_handler and self.dihedral_handler.mode == 'collect':
            self.all_dihedral_data.append(data['dihedrals'])
        
        # Write/Save dihedral angles list to file (streaming or collect)
        if self.dihedral_list_handler:
            for sub_idx, angles in enumerate(data['dihedral_list']):
                angle_str = ','.join(str(a) for a in angles)
                row = f"{frame_idx},{sub_idx},{angle_str}"
                self.dihedral_list_handler.append_row(row)
        
        # Write/Save cis-trans counts (streaming)
        if self.cistrans_handler:
            ct_dict = dict(data['cistrans'])
            cis_count = ct_dict.get('Cis', 0)
            trans_count = ct_dict.get('Trans', 0)
            row = f"{frame_idx},{cis_count},{trans_count}"
            self.cistrans_handler.append_row(row)
    
    def finalize_output(self):
        """
        Finalize all output files (write header and rows - 'collect' mode)
        """
        # Finalize dihedral angles counts file
        if self.dihedral_handler and self.dihedral_handler.mode == 'collect':
            # For collect mode, we need to merge dihedral data across frames
            self._finalize_dihedral_counts()
            
        # Finalize dihedral list file
        if self.dihedral_list_handler and self.dihedral_list_handler.mode == 'collect':
            header = "Frame,Subgraph,Dihedral Angles / °"
            self.dihedral_list_handler.finalize(header)
            
        # Finalize cis-trans counts file
        if self.cistrans_handler and self.cistrans_handler.mode == 'collect':
            header = "Frame,Cis count,Trans count"
            self.cistrans_handler.finalize(header)
    
    def _finalize_dihedral_counts(self):
        """
        Finalize dihedral angle counts by merging across frames.
        Creates a single output with Angle column once, followed by Frame columns.
        Uses pandas to create a properly formatted output.
        """
        # Build a DataFrame for each frame from collected dihedral data
        dfs = []
        
        # Process dihedral data for each frame
        for frame_idx, frame_data in enumerate(self.all_dihedral_data):
            # frame_data is a list of (angle, count) tuples
            if frame_data:
                df = pd.DataFrame(frame_data, columns=["Angle", f"Frame {frame_idx}"])
                dfs.append(df)
        
        if dfs:
            # Merge all DataFrames on the "Angle" column
            df_merged = dfs[0]
            for df in dfs[1:]:
                df_merged = pd.merge(df_merged, df, on="Angle", how="outer")
            
            # Sort by angle and fill NaN with 0, ensuring Angle column stays first
            df_merged = df_merged.sort_values("Angle").fillna(0)
            # Convert all frame columns to int, keep Angle as is
            for col in df_merged.columns:
                if col != "Angle":
                    df_merged[col] = df_merged[col].astype(int)
                    
            # write data with write raw from outoput handler
            csv_data = df_merged.to_csv(index=False)
            self.dihedral_handler.write_raw(csv_data)
    