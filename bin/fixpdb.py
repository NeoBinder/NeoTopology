import argparse

from openmm.app import PDBFile
from pdbfixer import PDBFixer

parser = argparse.ArgumentParser(description="fix protein pdb")
parser.add_argument("input", type=str, help="input pdb")
parser.add_argument("output", type=str, help="output pdb")
parser.add_argument(
    "--padding",
    type=float,
    default=1.0,
    help="padding box vectors in nanometer, default 1.0 nanometer",
)
parser.add_argument("--addH", type=bool, default=True, help="Need addHydrogens?")
parser.add_argument(
    "--pH",
    type=float,
    default=7.4,
    help="if addHydrogens, use pH values defines protonations. Default 7.4",
)
args = parser.parse_args()


def fix_protein(protein_path, padding=1.0 * unit.nanometer, pH_value=7.4, addH=True):
    protein_pdb = PDBFixer(filename=protein_path)
    protein_pdb.findMissingResidues()
    protein_pdb.findMissingAtoms()
    protein_pdb.findNonstandardResidues()
    protein_pdb.replaceNonstandardResidues()
    protein_pdb.addMissingAtoms()
    if addH:
        protein_pdb.addMissingHydrogens(pH_value)
    protein_pdb.removeHeterogens(False)
    print("Residues:", protein_pdb.missingResidues)
    print("Atoms:", protein_pdb.missingAtoms)
    print("Terminals:", protein_pdb.missingTerminals)
    print("Non-standard:", protein_pdb.nonstandardResidues)

    positions = []
    for pos in protein_pdb.positions:
        positions.append(pos.value_in_unit(unit.nanometer))

    positions = np.array(positions)
    box_vec = np.eye(3) * (
        (positions.max(0) - positions.min(0)).max()
        + padding.value_in_unit(unit.nanometer) * 2
    )

    protein_pdb.topology.setPeriodicBoxVectors(box_vec)

    print("Uses Periodic box:", protein_pdb.topology.getPeriodicBoxVectors())
    return protein_pdb


if __name__ == "__main__":
    fix_protein(
        args.input, padding=args.padding, pH_value=args.pH, addH=args.addHydrogens
    )
