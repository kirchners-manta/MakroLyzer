class EstimateFrames:

    @staticmethod
    def estimateFramesGROMACS(trajectory_path: str):
        """Return the MDAnalysis frame count for GRO, PDB, TRR, or XTC.

        Counting coordinates does not require a topology. GRO contributes one
        frame; PDB contributes one per model. TRR/XTC readers may build a frame
        offset index on first access, but do not load all coordinates into memory.
        """
        from MDAnalysis.coordinates.core import reader

        with reader(str(trajectory_path)) as trajectory:
            return trajectory.n_frames
    
    @staticmethod
    def estimateFramesXYZ(xyz_path: str):
        with open(xyz_path, 'r') as f:
            first_line = f.readline()
            try: 
                n_atoms = int(first_line.strip())
            except ValueError:
                raise ValueError("Could not read number of atoms from file.")
            # Get number of lines
            total_lines = sum(1 for _ in f)
        lines_per_frame = n_atoms + 2
        return round(total_lines/lines_per_frame)
    
    @staticmethod
    def estimateFramesLMP(lmp_path: str):
        with open(lmp_path, 'r') as f:
            for _ in range(3):
                next(f)
            n_atoms = f.readline()
            try:
                n_atoms = int(n_atoms.strip())
            except ValueError:
                raise ValueError("Could not read number of atoms from file.")
            for _ in range(5):
                next(f)
            total_lines = sum(1 for _ in f)
        lines_per_frame = n_atoms + 9
        return round(total_lines / lines_per_frame)
