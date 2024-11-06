import numpy as np
from .Rotation import RotationMatrix
from .dimer import find_interface


def get_center_of_mass(atoms):
    total_mass = 0.0
    com = np.zeros(3)
    for atom in atoms:
        atom_mass = atom.element.mass
        total_mass += atom_mass
        com += atom_mass * atom.position
    return com / total_mass
