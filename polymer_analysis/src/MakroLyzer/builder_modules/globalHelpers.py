"""Shared helpers for recipe-based structure builders."""

import json
from pathlib import Path


def read_recipe(recipe_file):
    """Read a UTF-8 JSON recipe and return its dictionary.

    Accepts a string or path-like file path. The JSON root must be an object;
    individual builders validate their own recipe fields.

    Example recipe::

        {"monomer": ["G", "G"], "linker": ["beta-O-4"]}

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


def minimize_energy():
    raise NotImplementedError("Energy minimization is not implemented yet.")
