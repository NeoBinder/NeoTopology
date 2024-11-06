# NeoTopology

NeoTopology is NeoBinder's First open soruce project for Protein Topology Toolkit All-In-One SDK.

NeoTopology is heavily designed for Protein Design, OpenMM Preparation & Study

This Package contains:
- IO from files including PDB/PDBx/rdkit
- Mathatics related with Rotation and Alignment
- Protein/Sequence/Ligand Calculations.

## Installation
NeoTopology can be installed with:

* source code installation
```bash
mkdir -p /path/to/project
cd /path/to/project
git clone git@github.com:NeoBinder/NeoTopology.git
cd /path/to/project/NeoTopology
# installation mode
pip install ./
# or devmode
pip install -e ./
```

* pypi
to be published

## Document
### NeoTopology Topology Object IO
```python
from ttk.io import topology_parser
# Load From PDB
top = topology_parser.topology_from_pdb(pdbpath)

# Load From PDB content
with open(pdbpath,"r") as f:
  content = f.read()
top = topology_parser.topology_from_pdb_content(content)

# Load From Openmm
top=topology_parser.topology_from_openmmm(openmmtop)

print(top.chains)
print(top.n_chains)
print(top.residues)
print(top.n_residues)
print(top.atoms)
print(top.bonds)

# Load Molecule From rdkit
top = topology_parser.topology_from_rdkitmol(mol,res_name)
print(top.residues)
print(top.bonds)

# Export to PDB
from ttk.io import PDBFile
from ttk.io import topology_export

# topology to pdb content
content = PDBFile().to_content(top)
if fname:
    with open(fname, "w") as f:
        f.write(content)

# topology to file directly
topology_export.topology_to_pdb(top,fname)
```

### Topology modification
```python
# topology add residue
top.add_residue(name,top.chains[0])
# chain add residue
top.chains[0].add_residue(res)
```


### Topology calculations
```python
from ttk.calculators import get_center_of_mass
res = top.get_residues_by_name(res_name)
com = get_center_of_mass(res.atoms)
```
