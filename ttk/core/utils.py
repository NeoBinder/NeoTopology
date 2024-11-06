import numbers

atom_types_dict = {
    "heavy_atoms": "heavy_atoms",
    "heavy": "heavy_atoms",
    "heavy_atom": "heavy_atoms",
    "default": "heavy_atoms",
    "all_atoms": None,
}


def filter_atom_types(atoms, atom_type):
    if atom_type in atom_types_dict:
        atom_type = atom_types_dict[atom_type]

    if atom_type == "heavy_atoms":
        atoms = [atom for atom in atoms if atom.element.symbol != "H"]
    return atoms


def to_unitless(coords):
    if isinstance(coords, numbers.Number):
        unit_defined = False
    else:
        coords = coords.to("nanometer").magnitude
        unit_defined = True
    return coords, unit_defined
