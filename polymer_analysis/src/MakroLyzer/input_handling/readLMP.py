import pandas as pd
import numpy as np

def readLMP(lmp_path: str, unwrap: bool = False, box_size: list = None):
    """
    Yield LAMMPS frames, optionally removing periodic jumps from x/y/z.

    xu/yu/zu are already unwrapped and always pass through unchanged.
    Temporal unwrapping requires a supplied fixed orthorhombic box and motion
    smaller than half a box length per axis between saved frames. Consume all
    frames before applying an analysis stride. The first frame is preserved.
    Atom IDs establish consistent ordering; without IDs, input order must be stable.
    """
    wrapped_previous = None
    unwrapped_previous = None
    elements_previous = None
    ids_previous = None
    mode_previous = None
    box = None
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
                    
            # Select a complete coordinate triplet. Image counters ix/iy/iz
            # are not coordinates and must never be used as position fallbacks.
            if all(name in header for name in ("xu", "yu", "zu")):
                coordinate_names = ("xu", "yu", "zu")
                already_unwrapped = True
            elif all(name in header for name in ("x", "y", "z")):
                coordinate_names = ("x", "y", "z")
                already_unwrapped = False
            else:
                raise ValueError("Expected complete x/y/z or xu/yu/zu coordinate columns.")
            x, y, z = (header.index(name) - 2 for name in coordinate_names)

            if unwrap:
                if mode_previous is not None and already_unwrapped != mode_previous:
                    raise ValueError("Coordinate convention changed between frames.")
                mode_previous = already_unwrapped
                # Native unwrapped coordinates need no box or minimum-image correction.
                if not already_unwrapped and box is None:
                    if box_size is None:
                        raise ValueError("Unwrapping x/y/z requires box_size.")
                    box = np.asarray(box_size, dtype=float)
                    if box.shape not in ((), (3,)) or not np.all(np.isfinite(box) & (box > 0)):
                        raise ValueError("box_size must be a positive scalar or three positive lengths.")

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
            if unwrap:
                has_ids = "id" in df.columns
                if has_ids:
                    if df["id"].duplicated().any():
                        raise ValueError("Duplicate atom IDs in LAMMPS frame.")
                    # Dump rows may be reordered between frames. Sort numerically
                    # so each array row continues to represent the same atom.
                    df = df.sort_values("id", key=lambda values: values.astype(int)).reset_index(drop=True)
                ids_current = tuple(df["id"]) if has_ids else None
                elements_current = df["atom"].to_numpy(copy=True)
                if elements_previous is not None:
                    if len(elements_current) != len(elements_previous):
                        raise ValueError("Atom count changed between frames during unwrapping.")
                    if ids_current != ids_previous:
                        raise ValueError("Atom IDs changed between frames during unwrapping.")
                    if not np.array_equal(elements_current, elements_previous):
                        raise ValueError("Element ordering changed between frames during unwrapping.")

                if not already_unwrapped:
                    wrapped_current = df[["x", "y", "z"]].to_numpy(copy=True)
                    if wrapped_previous is None:
                        unwrapped_current = wrapped_current.copy()
                    else:
                        displacement = wrapped_current - wrapped_previous
                        # Correct boundary jumps independently along each box axis.
                        displacement -= box * np.round(displacement / box)
                        # Retain all earlier crossings by accumulating on unwrapped positions.
                        unwrapped_current = unwrapped_previous + displacement
                    df[["x", "y", "z"]] = unwrapped_current
                    wrapped_previous = wrapped_current
                    # Caller edits to yielded frames must not alter retained state.
                    unwrapped_previous = unwrapped_current.copy()
                elements_previous = elements_current
                ids_previous = ids_current

            df["index"] = df.index
            df = df.reindex(columns=["index", "id", "atom", "x", "y", "z", "Molecule", "Charge"])
            
            yield df
