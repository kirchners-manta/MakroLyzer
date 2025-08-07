import pandas as pd

def readLMP(lmp_path: str):
    """
    Generator that yields one frame (as DataFrame) at a time from an extended LMP trajectory.
    """
    with open(lmp_path, 'r') as file:
        while True:
            try:
                for _ in range(3):  # skip the first 3 header lines
                    next(file)
            except StopIteration:
                break   
            num_atoms_line = file.readline()
            if not num_atoms_line:
                break  
            try:
                num_atoms = int(num_atoms_line.strip())
            except ValueError:
                raise ValueError(f"Expected number of atoms, got: {num_atoms_line}")

            try:
                for _ in range(4):  # skip the first 4 header lines
                    next(file)
            except StopIteration:
                break

            header = file.readline().strip().split()
            
            # Check if "element", "Element", or "type" is in the header
            try:
                atom_type_pos = header.index("element") - 2
            except ValueError:
                try:
                    atom_type_pos = header.index("Element") - 2
                except ValueError:
                    atom_type_pos = header.index("type") - 2
                    
            x_keywords = ["xu", "x", "ix"]
            y_keywords = ["yu", "y", "iy"]
            z_keywords = ["zu", "z", "iz"]
            
            def findIndex(header, keywords):
                for keyword in keywords:
                    if keyword in header:
                        return header.index(keyword) - 2
                raise ValueError("No position keyword found")
            
            x = findIndex(header, x_keywords)
            y = findIndex(header, y_keywords)
            z = findIndex(header, z_keywords)
            
            atom_mol_pos = header.index("mol") - 2 if "mol" in header else None
            atom_charge_pos = header.index("q") - 2 if "q" in header else None
            atom_id_pos = header.index("id") - 2 if "id" in header else None
            
            positions = {
                "id": atom_id_pos,
                "atom": atom_type_pos,
                "x": x,
                "y": y,
                "z": z,
                "Molecule": atom_mol_pos,
                "Charge": atom_charge_pos,
            }
            
            data = []
            for _ in range(num_atoms):
                line = file.readline()
                if not line:
                    break  
                parts = line.strip().split()
                if len(parts) < 3:
                    raise ValueError(f"Invalid LMP line: {line.strip()}")
                
                atom_data = {
                    "atom": parts[positions["atom"]],
                    "x": float(parts[positions["x"]]),
                    "y": float(parts[positions["y"]]),
                    "z": float(parts[positions["z"]]),
                }
                
                if positions["id"] is not None:
                    atom_data["id"] = parts[positions["id"]]
                if positions["Molecule"] is not None:
                    atom_data["Molecule"] = parts[positions["Molecule"]]
                if positions["Charge"] is not None:
                    atom_data["Charge"] = float(parts[positions["Charge"]])
                
                data.append(atom_data)
            if len(data) != num_atoms:
                break
            df = pd.DataFrame(data)
            df["index"] = df.index
            df = df[["index", "id", "atom", "x", "y", "z", "Molecule", "Charge"]]
            
            yield df