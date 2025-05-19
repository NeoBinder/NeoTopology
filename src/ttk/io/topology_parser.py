import gzip
from collections import Counter

import numpy as np

import ttk
from ttk import Topology
from ttk.core import UnitCell, element_from_symbol, from_bond_float, from_rdkit_bond
from ttk.io import PDBParser, PDBxParser


def topology_from_openmmm(
    mm_topology, positions=None, keepIdx=True, parse_protein_bond=True
):
    import openmm

    top = Topology()

    top.periodic_box = UnitCell.from_openmm(mm_topology, positions)
    atom_map = {}
    for chain in mm_topology.chains():
        c = top.add_chain(chain.id)
        for index, residue in enumerate(chain.residues()):
            res = top.add_residue(residue.name, c, residue.id)
            for atom in residue.atoms():
                if positions is None:
                    pos = None
                else:
                    pos = positions[atom.index]
                    pos = pos.value_in_unit(openmm.unit.nanometer)
                    pos = np.array(pos) * ttk.unit.nm

                if keepIdx:
                    idx = atom.index
                else:
                    idx = 0
                newatom = top.add_atom(
                    atom.name,
                    element_from_symbol(atom.element.symbol),
                    res,
                    pos,
                    index=idx,
                )
                atom_map[atom] = newatom
    for bond in mm_topology.bonds():
        top.add_bond(atom_map[bond[0]], atom_map[bond[1]], from_bond_float(bond.order))
    if not parse_protein_bond:
        for residue in top.residues:
            if residue.is_protein:
                for atom in residue.atoms:
                    atom.bonds = []
    return top


def load_from_file(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as f:
            content = f.read().decode().split("\n")
    else:
        with open(path, "r") as f:
            content = f.readlines()
    return content


def topology_from_pdb_content(content, **kwargs):
    parser = PDBParser(**kwargs)
    parser.load_from_content_list(content)
    top = parser.topologies[-1]
    return top


def topology_from_pdb(pdbpath, **kwargs):
    content = load_from_file(pdbpath)
    parser = PDBParser(**kwargs)
    parser.load_from_content_list(content)
    top = parser.topologies[-1]
    return top


def topology_from_pdbx_content(pdbx_path, **kwargs):
    content = load_from_file(pdbx_path)
    parser = PDBxReader(**kwargs)
    parser.parse_from_content(content)
    top = parser.topologies[-1]
    return top


def topology_from_rdkitmol(mol, res_name):
    top = Topology()
    chain = top.add_chain()
    res = top.add_residue(res_name, chain, res_seq=None, segment_id="")
    atom_counter = Counter()

    atoms = mol.GetAtoms()
    positions = mol.GetConformer().GetPositions()
    atom_map = {}

    for atom, pos in zip(atoms, positions):
        symbol = atom.GetSymbol()
        atom_counter.update(symbol)
        newatom = top.add_atom(
            symbol + str(atom_counter[symbol]),
            element_from_symbol(symbol),
            res,
            (pos * ttk.unit.angstrom).to("nm"),
        )
        atom_map[atom.GetIdx()] = newatom
    for bond in mol.GetBonds():
        top.add_bond(
            atom_map[bond.GetBeginAtomIdx()],
            atom_map[bond.GetEndAtomIdx()],
            from_rdkit_bond(bond),
        )
    return top
