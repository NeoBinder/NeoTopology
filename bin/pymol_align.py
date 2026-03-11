"""
use pymol align
---
return
RMSD after refinement
Number of aligned atoms after refinement
Number of refinement cycles
RMSD before refinement
Number of aligned atoms before refinement
Raw alignment score
Number of residues aligned
"""

import argparse
import os

from pymol import cmd

p = argparse.ArgumentParser()
p.add_argument("-m", type=str, required=True, help="mobile pdbf")
p.add_argument("-t", type=str, required=True, help="target pdbf")
p.add_argument("-o", type=str, required=False, help="output filename (no need for suffix)")
p.add_argument(
    "-mode",
    type=str,
    default="align",
    help="align:works when sequence similarity is relatively high, cealign, super:sequence-independent",
)

args = p.parse_args()

if __name__ == "__main__":
    cmd.load(args.m)
    cmd.load(args.t)

    # mobile_pdb = PDB.from_fpath(args.m)

    mobile = os.path.basename(args.m).rstrip(".pdb")
    target = os.path.basename(args.t).rstrip(".pdb")

    if args.o is None:
        if args.mode == "cealign":
            res = cmd.cealign(mobile, target)
        elif args.mode == "super":
            res = cmd.super(mobile, target)
        else:
            res = cmd.align(mobile, target)
    else:
        if args.mode == "cealign":
            res = cmd.cealign(mobile, target, object=args.o)
        elif args.mode == "super":
            res = cmd.super(mobile, target, object=args.o)
        else:
            res = cmd.align(mobile, target, object=args.o)

        cmd.save(args.o + ".pdb", selection=args.o)
        cmd.save(args.o + ".aln", selection=args.o)
        print(f"Saving to {args.o}")

    if args.mode == "cealign":
        print(res)
    else:
        (
            rmsd_after_refine,
            num_aln_atoms_after_refinement,
            num_cycles,
            rmsd_before_refine,
            num_aln_atoms_before_refinement,
            raw_aln_score,
            num_aln_res,
        ) = res
        print(
            f"----RMSD after refinement: {rmsd_after_refine} \n\
        Number of aligned atoms after refinement: {num_aln_atoms_after_refinement} \n\
        Number of refinement cycles: {num_cycles} \n\
        RMSD before refinement: {rmsd_before_refine} \n\
        Number of aligned atoms before refinement: {num_aln_atoms_before_refinement} \n\
        Raw alignment score: {raw_aln_score} \n\
        Number of residues aligned: {num_aln_res} "
        )
