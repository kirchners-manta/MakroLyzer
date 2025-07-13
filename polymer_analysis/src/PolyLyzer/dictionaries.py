"""Here you can find all dictionaries used in the code."""


# Dictionaries for element masses, van der Waals radii, and covalent radii
def dictMass() -> dict:
    elem_masses = {
        "H": 1.008,
        "Li": 6.941,
        "Na": 22.990,
        "K": 39.098,
        "Mg": 24.305,
        "Ca": 40.078,
        "Zn": 65.38,
        "C": 12.011,
        "N": 14.007,
        "O": 15.999,
        "B": 10.811,
        "F": 18.998,
        "P": 30.974,
        "S": 32.065,
        "Cl": 35.453,
        "Br": 79.904,
        "I": 126.904,
        "Ag": 107.868,
        "Au": 196.967,
        "D": 0.00,
        "X": 0.00,
    }
    return elem_masses

def dictVdW() -> dict:
    elem_vdW = {
        "H": 1.20,
        "Li": 1.81,
        "Na": 2.27,
        "K": 2.75,
        "Mg": 1.73,
        "Ca": 2.31,
        "Zn": 1.39,
        "C": 1.70,
        "N": 1.55,
        "O": 1.52,
        "B": 1.65,
        "F": 1.47,
        "P": 1.80,
        "S": 1.80,
        "Cl": 1.75,
        "Br": 1.85,
        "I": 1.98,
        "Ag": 1.72,
        "Au": 1.66,
        "D": 0.00,
        "X": 0.00,
    }
    return elem_vdW

def dictCovalent() -> dict:
    elem_covalent = {
        "H": 0.31,
        "Li": 1.28,
        "Na": 1.66,
        "K": 2.03,
        "Mg": 1.41,
        "Ca": 1.71,
        "Zn": 1.18,
        "C": 0.79,#0.76
        "N": 0.71,
        "O": 0.66,
        "B": 0.84,
        "F": 0.57,
        "P": 1.07,
        "S": 1.05,
        "Cl": 1.02,
        "Br": 1.20,
        "I": 1.39,
        "Ag": 1.45,
        "Au": 1.36,
        "D": 0.00,
        "X": 0.00,
    }
    return elem_covalent