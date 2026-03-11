import argparse
import re
from collections import defaultdict

parser = argparse.ArgumentParser(description="fix protein heteratom index")
parser.add_argument("input", type=str, help="input pdb")
parser.add_argument("output", type=str, help="output pdb")

args = parser.parse_args()
counter = defaultdict(int)
with open(args.input) as inputf, open(args.output, "w") as outputf:
    for line in inputf:
        if line.startswith("HETATM"):
            content = line.strip().split()
            element = content[-1]
            element = re.sub(r"\d+", "", element)
            element = re.sub(r"\+|\-", "", element)
            counter[element] += 1
            idx = f"{counter[element]:<2}"
            line = line[:14] + idx + line[16:]
        outputf.write(line)
