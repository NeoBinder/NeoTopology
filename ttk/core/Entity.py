import numbers
from .utils import atom_types_dict, filter_atom_types, to_unitless
import ttk
import numpy as np
from ttk import unit


class Entity(object):

    def __init__(self):
        self.atoms = []
        pass

    def __hash__(self):
        return id(self)

    @property
    def positions(self):
        return np.vstack([atom.position for atom in self.atoms])

    @positions.setter
    def positions(self, pos):
        atoms = self.atoms
        assert isinstance(pos, np.ndarray) or isinstance(pos, ttk.unit.Quantity)
        assert pos.shape[0] == len(atoms)
        assert pos.shape[1] == 3

        pos_value = np.copy(pos)
        if isinstance(pos_value, ttk.unit.Quantity):
            for idx, atom in enumerate(atoms):
                atom.position = pos_value[idx]
        else:
            print("no unit detected, use default unit (nm)")
            for idx, atom in enumerate(atoms):
                atom.position = pos_value[idx] * ttk.unit.nanometer

    def get_positions(
        self, heavy_atoms=False, unit=ttk.unit.nanometer, magnitude=False
    ):
        if heavy_atoms:
            atoms = self.heavy_atoms
        else:
            atoms = self.atoms

        res = np.vstack([atom.position for atom in atoms])
        res = res.to(unit)
        if magnitude:
            res = res.magnitude
        return res

    @property
    def heavy_atoms(self):
        return [atom for atom in self.atoms if atom.element.symbol != "H"]

    # maybe create select classes for any selections
    def select_atoms(self, idx_list):
        atoms_map = {atom.name: atom for atom in self.atoms}
        atoms = [atoms_map[idx] for idx in idx_list]
        return atoms

    def select_by_dist(self, coords, threshold, atom_type="default"):
        threshold, unit_defined = to_unitless(threshold)
        if not unit_defined:
            print("threshold unit not defined, use nanometer")
        coords, unit_defined = to_unitless(coords)
        if not unit_defined:
            print("coords unit not defined, use nanometer")

        selected_atoms = []

        if len(coords.shape) == 1:
            dim = 1
        elif len(coords.shape) == 2:
            dim = 2

        for atom in self.atoms:
            dist = np.linalg.norm(
                atom.position.to(ttk.unit.nanometer).magnitude - coords, axis=dim - 1
            )
            if dim == 2 and (dist < threshold).any():
                selected_atoms.append(atom)
            elif dim == 1 and dist < threshold:
                selected_atoms.append(atom)

        selected_atoms = filter_atom_types(selected_atoms, atom_type)

        return selected_atoms

    def is_close(self, entity, threshold):
        if isinstance(threshold, numbers.Number):
            threshold *= ttk.unit.nanometer
        target_positions = entity.get_positions()
        source_positions = self.get_positions()
        if target_positions.shape[1] != source_positions.shape[1]:
            raise ValueError(
                "Both input sets of observations must have the same number of features."
            )
        for i, xa in enumerate(source_positions):
            # Loop over all samples in set XB
            for j, xb in enumerate(target_positions):
                dist = np.sqrt(np.sum((xa - xb) ** 2))
                if dist < threshold:
                    return True
        return False
