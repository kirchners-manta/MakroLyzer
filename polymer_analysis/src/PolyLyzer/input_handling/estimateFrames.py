def estimateFrames(xyz_path: str):
    with open(xyz_path) as f:
        first_line = f.readline()
        try: 
            n_atoms = int(first_line.strip())
        except ValueError:
            raise ValueError("Could not read number of atoms from file.")
        # Get number of lines
        total_lines = sum(1 for _ in f)
    lines_per_frame = n_atoms + 2
    return round(total_lines/lines_per_frame)
        