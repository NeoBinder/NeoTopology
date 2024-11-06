import argparse

from openmm.app import PDBFile
from pdbfixer import PDBFixer

parser = argparse.ArgumentParser(description='fix protein pdb')
parser.add_argument('input', type=str, help='input pdb')
parser.add_argument('output', type=str, help='output pdb')

args = parser.parse_args()

fixer = PDBFixer(filename=args.input)
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(False)

fixer.findMissingAtoms()
fixer.addMissingHydrogens(7.4)
PDBFile.writeFile(fixer.topology, fixer.positions, open(args.output, 'w'))
