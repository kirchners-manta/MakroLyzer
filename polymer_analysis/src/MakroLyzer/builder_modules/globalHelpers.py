"""Shared helpers for recipe-based structure builders."""

import json
from pathlib import Path


def read_recipe(recipe_file):
    """Read a UTF-8 JSON recipe and return its dictionary.

    Accepts a string or path-like file path. The JSON root must be an object;
    individual builders validate their own recipe fields.

    Example recipe::

        {"monomer": ["G", "G"], "linker": ["4-O-beta"]}

    File access errors propagate to the caller. Invalid JSON or a non-object
    root raises ValueError with the recipe path for context.
    """
    path = Path(recipe_file)
    with path.open(encoding='utf-8') as file:
        try:
            recipe = json.load(file)
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"Invalid JSON in recipe '{path}' at line {exc.lineno}, "
                f"column {exc.colno}: {exc.msg}"
            ) from exc

    if not isinstance(recipe, dict):
        raise ValueError(f"Recipe '{path}' must contain a JSON object.")
    return recipe


def optimize_structure(graph, max_iterations=1500):
    """Return a UFF-relaxed graph copy, retaining IDs, attributes and bonds.

    Requires RDKit and explicit hydrogens. Small reproducible coordinate noise
    breaks exactly linear initial angles. Store convergence and energies in
    graph.graph['uff_optimization']; energies are in kcal/mol. If the iteration
    limit is reached, warn and return the latest geometry.
    """
    import warnings
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem

    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer.")
    nodes = list(graph)
    if not nodes:
        raise ValueError("Cannot optimize an empty graph.")
    # Graph IDs may contain gaps after atom deletion. RDKit instead numbers
    # atoms consecutively from 0; this mapping connects the two representations.
    # The same node order is used later to copy coordinates back correctly.
    indices = {node: index for index, node in enumerate(nodes)}
    # RWMol is an editable RDKit molecule: first populate atoms, then bonds.
    molecule = Chem.RWMol()
    positions = []
    for node in nodes:
        data = graph.nodes[node]
        atom = Chem.Atom(data['element'])
        atom.SetNoImplicit(True)  # Do not invent hydrogens absent from the graph.
        # Use a stored formal charge, or neutral atoms when none is specified.
        atom.SetFormalCharge(int(data.get('formal_charge', 0)))
        molecule.AddAtom(atom)
        positions.append([data[key] for key in ('x', 'y', 'z')])
    # Transfer explicit connectivity; never infer bonds from atom distances.
    # left/right are the two graph node IDs at the ends of each bond.
    bond_types = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    for left, right, data in graph.edges(data=True):
        order = data.get('bond_order')
        if order not in bond_types:
            raise ValueError(f"Missing or unsupported bond_order {order!r} on {(left, right)}.")
        molecule.AddBond(indices[left], indices[right], bond_types[order])
    # Convert to a regular molecule and let RDKit assign chemical properties
    # such as valence, aromaticity and hybridization, needed for force-field typing.
    molecule = molecule.GetMol()
    Chem.SanitizeMol(molecule)
    # This helper is intended for closed-shell builder products. An incomplete
    # valence can be interpreted as a radical, so reject it before optimization.
    if any(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()):
        raise ValueError("Expected a closed-shell structure; check atom valences.")
    # UFF needs parameters for every atom type, not just an element name.
    if not AllChem.UFFHasAllMoleculeParams(molecule):
        raise ValueError("UFF parameters are unavailable for one or more atoms.")

    positions = np.asarray(positions, dtype=float)
    if not np.isfinite(positions).all():
        raise ValueError("Optimization requires finite coordinates.")
    # Add tiny displacements (standard deviation 0.001 angstrom) so exactly
    # linear angles can bend. Seed 0 makes this perturbation reproducible.
    positions += np.random.default_rng(0).normal(0, 0.001, positions.shape)
    # A conformer stores one set of 3D coordinates for the molecule's atoms.
    # Constructed geometry is used directly, rather than generating a new one.
    conformer = Chem.Conformer(len(nodes))
    conformer.Set3D(True)
    for index, position in enumerate(positions):
        conformer.SetAtomPosition(index, tuple(map(float, position)))
    molecule.AddConformer(conformer)
    # Include interactions between disconnected fragments if any are present.
    # All atoms remain free to move; no positional constraints are imposed.
    force_field = AllChem.UFFGetMoleculeForceField(molecule, ignoreInterfragInteractions=False)
    if force_field is None:
        raise ValueError("Could not construct the UFF force field.")
    # Initial energy refers to the slightly perturbed starting coordinates.
    # Minimize updates the conformer in place: 0 means converged, 1 means
    # the iteration limit was reached and more iterations may be required.
    initial_energy = float(force_field.CalcEnergy())
    status = force_field.Minimize(maxIts=max_iterations)
    final_energy = float(force_field.CalcEnergy())
    optimized = molecule.GetConformer().GetPositions()
    if not np.isfinite(optimized).all() or not np.isfinite([initial_energy, final_energy]).all():
        raise RuntimeError("UFF produced nonfinite coordinates or energies.")

    # Keep the original graph untouched. Transfer only x/y/z; retain its atom
    # labels, IDs, bond orders and all other existing attributes.
    result = graph.copy()
    for node, position in zip(nodes, optimized):
        result.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))
    # Store diagnostics on the graph itself, separate from individual atoms.
    # These describe this optimization call, not the entire construction history.
    result.graph['uff_optimization'] = {
        'converged': status == 0, 'status': status,
        'initial_energy': initial_energy, 'final_energy': final_energy,
        'max_iterations': max_iterations,
    }
    # A finite geometry may still be useful without convergence. Return it with
    # an explicit warning instead of silently treating the step as successful.
    if status != 0:
        warnings.warn(
            f"UFF did not converge within {max_iterations} iterations; returning latest geometry.",
            RuntimeWarning, stacklevel=2,
        )
    return result
