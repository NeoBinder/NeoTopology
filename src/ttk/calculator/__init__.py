import numpy as np

from .dimer import find_interface


def get_center_of_mass(atoms):
    total_mass = 0.0
    com = np.zeros(3)
    for atom in atoms:
        atom_mass = atom.element.mass
        total_mass += atom_mass
        com += atom_mass * atom.position
    if total_mass == 0:
        return np.full(3, np.nan)  # 返回 [nan, nan, nan]
    return com / total_mass
