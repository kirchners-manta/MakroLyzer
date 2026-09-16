"""Read GROMACS coordinates into the same frame format as the XYZ reader."""

from functools import partial
import warnings

import pandas as pd


def readGROMACS(topology_path: str, trajectory_path: str = None):
    """Yield independent frames with element labels and coordinates in Ångström.

    GRO/PDB files can be read alone; TRR/XTC files require a topology with
    matching atom order. Missing elements are guessed from atom names by
    MDAnalysis. Bonds and periodic boxes are still handled by MakroLyzer.
    """
    import MDAnalysis as mda
    from MDAnalysis.exceptions import NoDataError

    paths = (str(topology_path),) if trajectory_path is None else (
        str(topology_path), str(trajectory_path)
    )
    universe = mda.Universe(*paths, convert_units=True)
    try:
        try:
            elements = universe.atoms.elements
            needs_guess = any(not str(element).strip() for element in elements)
        except NoDataError:
            needs_guess = True
        if needs_guess:
            warnings.warn(
                "Missing element information: MDAnalysis is guessing elements from "
                "atom names. Check these assignments, especially for nonstandard "
                "names or coarse-grained systems.",
                UserWarning,
                stacklevel=2,
            )
            universe.guess_TopologyAttrs(to_guess=['elements'])
        elements = [str(element).strip().capitalize() for element in universe.atoms.elements]
        if any(not element for element in elements):
            raise ValueError("Could not determine all elements. Supply a topology with element information.")

        for timestep in universe.trajectory:
            # Readers reuse their coordinate buffers; retained frames must own their data.
            frame = pd.DataFrame(timestep.positions.copy(), columns=['x', 'y', 'z'])
            frame.insert(0, 'atom', elements)
            frame['index'] = universe.atoms.indices
            yield frame
    finally:
        universe.trajectory.close()


def readGROMACS_trr(tpr_path: str, trr_path: str):
    """Yield frames from a topology and TRR trajectory."""
    yield from readGROMACS(tpr_path, trr_path)


def readGROMACS_gro(gro_path: str):
    """Yield the structure in a GRO file (MDAnalysis reads its first frame)."""
    yield from readGROMACS(gro_path)


def readGROMACS_pdb(pdb_path: str):
    """Yield each model in a PDB file."""
    yield from readGROMACS(pdb_path)


def input_reader(args):
    """Return the input path and reader for either processing pipeline."""
    for key in ('groFile', 'pdbFile'):
        if args.get(key):
            return args[key], readGROMACS
    for key in ('trrFile', 'xtcFile'):
        if args.get(key):
            topology = args.get('topologyFile')
            if not topology:
                raise ValueError('TRR/XTC input requires -tpr/--topology.')
            return args[key], partial(readGROMACS, topology)
    raise ValueError('No supported input file was provided.')
