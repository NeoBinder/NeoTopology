import numpy as np


def find_interface(chain1, chain2, threshold=0.3):
    pos_c1 = chain1.positions.magnitude
    pos_c2 = chain2.positions.magnitude

    c1_interface_residues = set()
    c1_interface_residues_index = set()
    c2_interface_residues = set()
    c2_interface_residues_index = set()
    for atom in chain1.atoms:
        if (np.linalg.norm(pos_c2 - atom.position.magnitude, axis=1) < threshold).any():
            if atom.residue.index not in c1_interface_residues_index:
                c1_interface_residues_index.add(atom.residue.index)
                c1_interface_residues.add(atom.residue)

    for atom in chain2.atoms:
        if (np.linalg.norm(pos_c1 - atom.position.magnitude, axis=1) < threshold).any():
            if atom.residue.index not in c2_interface_residues_index:
                c2_interface_residues_index.add(atom.residue.index)
                c2_interface_residues.add(atom.residue)
    return c1_interface_residues, c2_interface_residues


def have_interface_ligands(
    topology, mass_threshold=100.0, distance_threshold=5.0, sequence_length_threshold=50
):

    def filter_ligand(res, mass_threshold):
        if res.is_protein or res.is_water or res.mass < mass_threshold:
            return False
        return True

    filtered_ligands = list(
        filter(lambda x: filter_ligand(x, mass_threshold), topology.residues)
    )
    if len(filtered_ligands) == 0:
        return False
    filtered_chains = [
        chain
        for chain in topology.chains
        if chain.is_protein_chain() and len(chain.residues) > sequence_length_threshold
    ]
    if len(filtered_chains) != 2:
        # print(filtered_chains)
        # print("{} is not dimer".format(pdbid))
        return False
    for ligand in filtered_ligands:
        if filtered_chains[0].is_close(ligand, distance_threshold) and filtered_chains[
            1
        ].is_close(ligand, distance_threshold):
            return True
