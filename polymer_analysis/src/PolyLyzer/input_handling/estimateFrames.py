class EstimateFrames:
    
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