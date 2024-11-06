import argparse
import re
from collections import defaultdict

from openmm.app import PDBFile

parser = argparse.ArgumentParser(description='fix protein heteratom index')
parser.add_argument('input', type=str, help='input pdb')
parser.add_argument('output', type=str, help='output pdb')

args = parser.parse_args()
counter = defaultdict(int)
with open(args.input,"r") as inputf:
    with open(args.output,"w") as outputf:
        for line in inputf:
            if line.startswith("HETATM"):
                content = line.strip().split()
                element = content[-1]
                element = re.sub("\d+","",element)
                element = re.sub("\+|\-","",element)
                counter[element] +=1
                idx = "{:<2}".format(counter[element])
                line = line[:14] + idx + line[16:]
            outputf.write(line)


