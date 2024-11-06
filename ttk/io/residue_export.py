from rdkit import Chem
from rdkit.Chem import AllChem

import ttk


def res2mol(res, removeHs=False, template=None):
    mol = Chem.RWMol()
    atom2molidx = {}
    bonds = set()
    conf = Chem.rdchem.Conformer(len(res.atoms))
    for atom in res.atoms:
        mol_atom = Chem.Atom(atom.symbol)
        idx = mol.AddAtom(mol_atom)
        atom2molidx[atom] = idx
        for bond in atom.bonds:
            bonds.add(bond)
        conf.SetAtomPosition(idx, atom.get_position(ttk.unit.angstrom, magnitude=True))

    for bond in bonds:
        if bond.atom1 in atom2molidx and bond.atom2 in atom2molidx:
            mol.AddBond(
                atom2molidx[bond.atom1], atom2molidx[bond.atom2], bond.type.to_rdkit()
            )
        else:
            if bond.atom1 not in atom2molidx:
                source = atom2molidx[bond.atom2]
                atom = bond.atom1
                symbol = atom.symbol
            else:
                source = atom2molidx[bond.atom1]
                atom = bond.atom2
                symbol = atom.symbol

            mol.GetAtomWithIdx(source).SetProp("conn", symbol)
            if symbol not in ttk.data._Metal:
                idx = mol.AddAtom(Chem.Atom(symbol))
                conf.SetAtomPosition(
                    idx, atom.get_position(ttk.unit.angstrom, magnitude=True)
                )
                mol.GetAtomWithIdx(idx).SetProp("atomLabel", "*")
                mol.AddBond(source, idx, bond.type.to_rdkit())
    mol = mol.GetMol()
    mol.AddConformer(conf)
    if removeHs:
        mol = Chem.RemoveHs(mol)
    if template:
        # smiles template
        if isinstance(template, str):
            template = Chem.MolFromSmiles(template)
        if isinstance(template, Chem.rdchem.Mol):
            mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    return mol
