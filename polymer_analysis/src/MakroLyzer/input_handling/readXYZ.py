import pandas as pd
import numpy as np  

def readXYZ(xyz_path: str, unwrap: bool = False, box_size: list = None):
    """
    Yield independent XYZ frames, optionally removing periodic jumps over time.

    Unwrapping requires a fixed orthorhombic box (one length or three lengths),
    unchanged atom ordering, and motion smaller than half a box length per axis
    between saved frames. The first frame is preserved, not made whole.
    Consume every frame before applying any analysis stride.
    """
    wrapped_previous = None
    unwrapped_previous = None
    elements_previous = None
    if unwrap:
        if box_size is None:
            raise ValueError("Unwrapping requires box_size.")
        box = np.asarray(box_size, dtype=float)
        if box.shape not in ((), (3,)) or not np.all(np.isfinite(box) & (box > 0)):
            raise ValueError("box_size must be a positive scalar or three positive lengths.")

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
            if unwrap:
                # Each row represents the same atom across successive frames.
                wrapped_current = df[["x", "y", "z"]].to_numpy(copy=True)
                elements_current = df["atom"].to_numpy(copy=True)
                if wrapped_previous is None:
                    unwrapped_current = wrapped_current.copy()
                else:
                    if wrapped_current.shape != wrapped_previous.shape:
                        raise ValueError("Atom count changed between frames during unwrapping.")
                    # XYZ has no IDs: swaps of identical elements cannot be detected.
                    if not np.array_equal(elements_current, elements_previous):
                        raise ValueError("Element ordering changed between frames during unwrapping.")
                    displacement = wrapped_current - wrapped_previous
                    # Remove periodic jumps using the nearest image per axis.
                    displacement -= box * np.round(displacement / box)
                    # Accumulate on unwrapped positions to retain earlier crossings.
                    unwrapped_current = unwrapped_previous + displacement

                df[["x", "y", "z"]] = unwrapped_current
                # Retain independent arrays so caller edits cannot change our state.
                wrapped_previous = wrapped_current
                unwrapped_previous = unwrapped_current.copy()
                elements_previous = elements_current
            yield df
