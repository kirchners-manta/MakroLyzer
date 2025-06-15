import pandas as pd

def readXYZ(xyz_path: str):
    """
    Generator that yields one frame (as DataFrame) at a time from an extended XYZ trajectory.
    """
    with open(xyz_path, 'r') as file:
        while True:
            num_atoms_line = file.readline()
            if not num_atoms_line:
                break  
            try:
                num_atoms = int(num_atoms_line.strip())
            except ValueError:
                raise ValueError(f"Expected number of atoms, got: {num_atoms_line}")
            # skip comment line
            file.readline() 
            
            data = []
            for _ in range(num_atoms):
                line = file.readline()
                if not line:
                    break  
                parts = line.strip().split()
                if len(parts) == 4:
                    atom, x, y, z = parts
                    data.append([atom, float(x), float(y), float(z)])
                else:
                    raise ValueError(f"Invalid XYZ line: {line.strip()}")
                    
            if len(data) != num_atoms:
                break  # incomplete frame
            df = pd.DataFrame(data, columns=["atom", "x", "y", "z"])
            df["index"] = df.index
            yield df
