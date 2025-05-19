"""
process deepalign outputs
"""

import argparse

import numpy as np

import ttk
from ttk.io import topology_parser


def get_pocket(top, threshold):
    # get HETATM
    hetatms = [
        atom
        for atom in top.atoms
        if (atom.residue.is_protein is False) and (atom.residue.is_water is False)
    ]
    print("Found a total of {} hetatms".format(len(hetatms)))

    # select by distance
    atoms = []
    for atm in hetatms:
        # coord = atm.position.magnitude
        atoms.extend(
            top.select_by_dist(atm.position, threshold * ttk.unit.angstrom, "heavy")
        )

    residues = set([atom.residue.res_seq for atom in atoms if atom.residue.is_protein])
    print("Found a total of {} residues".format(len(residues)))
    print(residues)
    return residues


def calculate_pocket_rmsd(top, align_res_dict, selections):
    t_chain, s_chain = top.chains[0], top.chains[1]
    selected_t_res = []
    selected_s_res = []
    s_selections = []
    for t_res in t_chain.residues:
        if t_res.res_seq in selections:
            if align_res_dict[t_res.res_seq][2] != "-":
                selected_t_res.append(t_res)
                s_selections.append(align_res_dict[t_res.res_seq][1])

    s_res_dict = {}
    for s_res in s_chain.residues:
        s_res_dict[s_res.res_seq] = s_res

    selected_s_res = [s_res_dict[x] for x in s_selections]

    assert len(selected_t_res) == len(selected_s_res)

    rmsd_ca = {}
    for i, t_res in enumerate(selected_t_res):
        t_ca = (
            list(t_res.atoms_by_name("CA"))[0].position.to(ttk.unit.angstrom).magnitude
        )
        s_ca = (
            list(selected_s_res[i].atoms_by_name("CA"))[0]
            .position.to(ttk.unit.angstrom)
            .magnitude
        )
        rmsd_ca[t_res.res_seq] = np.linalg.norm(t_ca - s_ca)

    # calculate all atom rmsd for same residues
    rmsd_aa = {}
    same_res = []
    for i, t_res in enumerate(selected_t_res):
        if selected_t_res[i].name != selected_s_res[i].name:
            print("Warning: residues do not match")
            continue
        same_res.append(t_res.res_seq)
        t_atoms = [
            atom.position.to(ttk.unit.angstrom).magnitude for atom in t_res.atoms
        ]
        s_atoms = [
            atom.position.to(ttk.unit.angstrom).magnitude
            for atom in selected_s_res[i].atoms
        ]
        rmsd_aa[t_res.res_seq] = np.linalg.norm(np.array(t_atoms) - np.array(s_atoms))
    print("Same residues: ", same_res)

    return rmsd_ca, rmsd_aa


def get_alignment_from_local(localf):
    """
    return target idx: (target residue, source index, source residue)
    """
    with open(localf, "r") as f:
        lines = f.readlines()

    target, source = [], []
    for line in lines:
        if line == "\n" or line.startswith("T cle") or line.startswith("S cle"):
            continue
        if line.startswith("T "):
            target.append(line)
        elif line.startswith("S "):
            source.append(line)

    align_res_dict = {}
    for t, s in zip(target, source):
        seq_t = t.split()[3]
        seq_s = s.split()[3]
        assert len(seq_t) == len(seq_s)
        start_t = int(t.split()[2])
        start_s = int(s.split()[2])
        t_idx = 0
        s_idx = 0
        for i in range(len(seq_t)):
            align_res_dict[start_t + t_idx] = (seq_t[i], start_s + s_idx, seq_s[i])
            if seq_t[i] != "-":
                t_idx += 1
            if seq_s[i] != "-":
                s_idx += 1
    return align_res_dict


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "-alignfile", type=str, required=True, help="deepalign output prefix"
    )
    p.add_argument(
        "-pocket_threshold",
        type=float,
        default=6,
        help="threshold distance to define pocket",
    )
    p.add_argument(
        "-pdb_for_pocket",
        type=str,
        required=True,
        help="crystal structure to find pocket",
    )
    args = p.parse_args()
    print(args)

    original_top = topology_parser.topology_from_pdb(args.pdb_for_pocket)

    # process input pdb to get pocket residues
    pocket_residues = get_pocket(original_top, args.pocket_threshold)
    print("pocket residues", pocket_residues)

    # process local file to get alignment between residues
    align_res_dict = get_alignment_from_local(args.alignfile + ".local")
    # print("return target idx: (target residue, source index, source residue)")
    # print(align_res_dict)

    # load aligned structure
    aligned_top = topology_parser.topology_from_pdb(args.alignfile + ".pdb")

    rmsd_ca, rmsd_aa = calculate_pocket_rmsd(
        aligned_top, align_res_dict, pocket_residues
    )

    print("=======Pocket stats=======")
    print("Mean RMSD between carbon alpha: {}".format(np.mean(list(rmsd_ca.values()))))
    print(
        "Mean RMSD between all atoms (matched same residues only): {}".format(
            np.mean(list(rmsd_aa.values()))
        )
    )
