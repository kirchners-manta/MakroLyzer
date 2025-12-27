"""This module contains functions for generating error outputs for MakroLyzer analyses."""

class ErrorOutputs:
    NEGATIVE_G_EIGENVALUE_WARNING = (
        "Warning: Negative eigenvalue encountered in the radius of gyration tensor calculation. "
        "You might want to check your input structure."
    )
    ZERO_G_EIGENVALUE_WARNING = (
        "Warning: All gyration tensor eigenvalues are zero. "
        "You might want to check your input structure."
    )