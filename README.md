# NeoTopology

NeoTopology is NeoBinder's First open soruce project for Protein Topology Toolkit All-In-One SDK.

NeoTopology is heavily designed for Protein Design, OpenMM Preparation & Study

This Package contains:
- IO from files including PDB/PDBx/rdkit
- Mathematics related sdk like Rotation and Alignment
- Protein/Sequence/Ligand Calculations

## Installation
NeoTopology can be installed with:

* uv (recommended)
1. Install [uv](https://docs.astral.sh/uv/getting-started/installation/)
2. Install package
```bash
git clone git@github.com:NeoBinder/NeoTopology.git
cd NeoTopology
uv venv -p 3.11
source ./.venv/bin/activate
uv sync
```

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
```
pip install neotopology
```

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
res = top.add_residue(name,top.chains[0],res)
# do some modifications on the resdiue
res.add_atom(atom)

# chain add residue
top.chains[0].add_residue(res)
```


### Topology calculations
```python
from ttk.calculators import get_center_of_mass
res = top.get_residues_by_name(res_name)
com = get_center_of_mass(res.atoms)
```

### 查找二聚体接口
```python
from ttk.calculator.dimer import find_interface

# 获取两个链之间的接口残基
chain1, chain2 = top.chains[0], top.chains[1]
chain1_interface, chain2_interface = find_interface(chain1, chain2, threshold=0.3)

print(f"Chain 1 interface residues: {len(chain1_interface)}")
print(f"Chain 2 interface residues: {len(chain2_interface)}")
```

### 计算键角和二面角
```python
from ttk.math.vector import angle_between, calc_dihedral
import ttk.unit as unit

# 计算两个向量之间的角度
import numpy as np
v1 = np.array([1.0, 0.0, 0.0])
v2 = np.array([0.0, 1.0, 0.0])
angle = angle_between(v1, v2)
print(f"Angle: {angle} radians ({angle * 180 / np.pi} degrees)")

# 计算二面角（例如蛋白质主链二面角）
p0 = np.array([0.0, 0.0, 0.0])
p1 = np.array([1.0, 0.0, 0.0])
p2 = np.array([2.0, 1.0, 0.0])
p3 = np.array([3.0, 1.0, 1.0])
dihedral = calc_dihedral(p0, p1, p2, p3)
print(f"Dihedral angle: {dihedral} radians ({dihedral * 180 / np.pi} degrees)")
```

### 删除溶剂和水分子
```python
# 只删除水分子
top.remove_water()
print(f"After removing water: {top.n_residues} residues")

# 删除所有溶剂（包括离子等）
top.remove_solvent()
print(f"After removing solvent: {top.n_residues} residues")
```

### 按名称获取残基
```python
# 获取单个残基名称
ala_residues = list(top.get_residues_by_name("ALA"))
print(f"Found {len(ala_residues)} ALA residues")

# 获取多个残基名称
protein_residues = list(top.get_residues_by_name(["ALA", "GLY", "VAL"]))
print(f"Found {len(protein_residues)} protein residues")
```

### 创建和索引拓扑
```python
# 为所有原子和残基创建索引
top.create_index()

# 检查索引是否正确
for atom in top.atoms:
    print(f"Atom {atom.name} has index {atom.index}")

for residue in top.residues:
    print(f"Residue {residue.name} has index {residue.index}")
```
