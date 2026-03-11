# NeoTopology

<p align="center">
  <a href="https://github.com/NeoBinder/NeoTopology/actions/workflows/test.yml">
    <img src="https://github.com/NeoBinder/NeoTopology/actions/workflows/test.yml/badge.svg" alt="Tests">
  </a>
  <a href="https://pypi.org/project/neotopology/">
    <img src="https://img.shields.io/pypi/v/neotopology.svg" alt="PyPI version">
  </a>
  <a href="https://pypi.org/project/neotopology/">
    <img src="https://img.shields.io/pypi/pyversions/neotopology.svg" alt="Python Versions">
  </a>
  <a href="https://github.com/NeoBinder/NeoTopology/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="License">
  </a>
</p>

**NeoTopology** is NeoBinder's all-in-one Python SDK for protein topology — purpose-built for
protein design, OpenMM preparation, and structural biology research.

## Features

- **File I/O** — Read and write PDB, PDBx (mmCIF), and RDKit molecule formats
- **Hierarchical data model** — `Topology → Chain → Residue → Atom` with full traversal APIs
- **Mathematics** — Rotation matrices, rigid-body alignment (SVD), dihedral/bond angle calculations
- **Calculators** — Center of mass, RMSD, dimer interface detection
- **Unit-safe coordinates** — All positions carry physical units via [pint](https://pint.readthedocs.io/) (default: nanometers)
- **OpenMM integration** — Convert topologies to/from OpenMM for molecular dynamics
- **CLI utilities** — Fix PDB files, PyMOL-based structural alignment

---

## Installation

### From PyPI (recommended)

```bash
pip install neotopology
```

### With [uv](https://docs.astral.sh/uv/) (fastest)

```bash
git clone https://github.com/NeoBinder/NeoTopology.git
cd NeoTopology
uv venv -p 3.11
source .venv/bin/activate
uv sync
```

### From source

```bash
git clone https://github.com/NeoBinder/NeoTopology.git
cd NeoTopology
pip install -e .          # editable install
# or
pip install .             # regular install
```

> **Requirements:** Python >= 3.11, numpy, pint, rdkit, mendeleev

---

## Quick Start

```python
from ttk.io import topology_parser
from ttk.io.topology_export import topology_to_pdb

# Load a protein from a PDB file
top = topology_parser.topology_from_pdb("protein.pdb")

print(f"Chains:   {top.n_chains}")
print(f"Residues: {top.n_residues}")
print(f"Atoms:    {len(top.atoms)}")

# Export back to PDB
topology_to_pdb(top, "output.pdb")
```

---

## Usage Guide

### Loading Topologies

```python
from ttk.io import topology_parser

# From a PDB file path
top = topology_parser.topology_from_pdb("protein.pdb")

# From a gzip-compressed PDB file
top = topology_parser.topology_from_pdb("protein.pdb.gz")

# From PDB text content
with open("protein.pdb") as f:
    content = f.read()
top = topology_parser.topology_from_pdb_content(content)

# From an RDKit molecule
from rdkit import Chem
mol = Chem.MolFromSmiles("c1ccccc1")          # benzene
top = topology_parser.topology_from_rdkitmol(mol, res_name="BNZ")

# From an OpenMM topology (requires openmm)
import openmm.app as app
pdb = app.PDBFile("protein.pdb")
top = topology_parser.topology_from_openmmm(pdb.topology, pdb.positions)
```

### Exporting Topologies

```python
from ttk.io.topology_export import topology_to_pdb
from ttk.io import PDBFile

# Write a PDB file
topology_to_pdb(top, "output.pdb")

# Get PDB text content (no file written)
content = topology_to_pdb(top)

# Using PDBFile directly
content = PDBFile().to_content(top)
with open("output.pdb", "w") as f:
    f.write(content)
```

### Traversing the Hierarchy

```python
# Chains
for chain in top.chains:
    print(f"Chain {chain.id}: {chain.n_residues} residues, {chain.n_atoms} atoms")

# Residues
for res in top.residues:
    if res.is_protein:
        print(f"  {res.name} ({res.code}) seq={res.res_seq}")
    elif res.is_water:
        print(f"  Water HOH at seq={res.res_seq}")

# Atoms
for atom in top.atoms:
    print(f"  {atom.name} ({atom.element.symbol}) index={atom.index}")
    if atom.position is not None:
        print(f"    position: {atom.position}")

# Bonds
for bond in top.bonds:
    print(f"  {bond.atom1.name} -- {bond.atom2.name}  type={bond.type.name}")
```

### Building Topologies Programmatically

```python
from ttk import Topology, unit
from ttk.core import element_from_symbol

top = Topology()

# Add a chain (auto-assigned ID "A", "B", ...)
chain = top.add_chain()
# Or specify an ID explicitly
chain_b = top.add_chain("B")

# Add a residue to the chain
res = top.add_residue("ALA", chain, res_seq=1)

# Add atoms with 3-D coordinates (positions in nanometers)
n_atom = top.add_atom(
    "N",
    element_from_symbol("N"),
    res,
    is_hetero=False,
    position=unit.Quantity([-0.227, 0.000, 0.000], unit.nanometer),
)
ca_atom = top.add_atom(
    "CA",
    element_from_symbol("C"),
    res,
    is_hetero=False,
    position=unit.Quantity([-0.171, 0.140, 0.000], unit.nanometer),
)

# Add a covalent bond
top.add_bond(n_atom, ca_atom)

# Assign sequential indices to all atoms and residues
top.create_index()
```

### Filtering and Querying

```python
# Filter residues by name (single name, list, or set)
ala_residues = list(top.get_residues_by_name("ALA"))
polar_residues = list(top.get_residues_by_name(["SER", "THR", "ASN", "GLN"]))

# Get atoms by index (call create_index() first)
top.create_index()
atoms = top.get_atom_by_indices([0, 1, 5])

# Access an atom in a residue by clean name
ca = res["CA"]          # raises KeyError if not found
cb = res.get("CB")      # returns None if not found
```

### Removing Solvent

```python
# Remove only water molecules
top.remove_water()

# Remove all solvent (water + ions + simple molecules)
top.remove_solvent()
```

### Merging Topologies

```python
top1 = topology_parser.topology_from_pdb("chain_A.pdb")
top2 = topology_parser.topology_from_pdb("chain_B.pdb")

top1.add_topology(top2)   # modifies top1 in place
print(top1.n_chains)      # combined chain count
```

---

## Calculators

### Center of Mass

```python
from ttk.calculator import get_center_of_mass

# Compute center of mass for a list of atoms (returns pint Quantity)
com = get_center_of_mass(res.atoms)
print(com)                     # e.g. [1.234, 5.678, 0.012] nanometer
print(com.to("angstrom"))      # convert to Angstrom
```

### Dimer Interface Detection

```python
from ttk.calculator.dimer import find_interface, have_interface_ligands

chain_a, chain_b = top.chains[0], top.chains[1]

# Returns two lists of interface residues (one per chain)
# threshold is the contact distance cutoff in nanometers
iface_a, iface_b = find_interface(chain_a, chain_b, threshold=0.5)

print(f"Chain A interface residues: {len(iface_a)}")
print(f"Chain B interface residues: {len(iface_b)}")

# Check whether any ligands sit at the protein-protein interface
has_ligands = have_interface_ligands(top)
```

---

## Mathematics

### Vector Operations

```python
import numpy as np
from ttk.math.vector import angle_between, calc_dihedral, unit_vector

# Normalize a vector
v = np.array([3.0, 4.0, 0.0])
uv = unit_vector(v)           # [0.6, 0.8, 0.0]

# Angle between two vectors (radians)
v1 = np.array([1.0, 0.0, 0.0])
v2 = np.array([0.0, 1.0, 0.0])
angle = angle_between(v1, v2)           # pi/2
print(f"{np.degrees(angle):.1f} deg")  # 90.0 deg

# Dihedral / torsion angle for four points (result in radians, range (-pi, pi])
p0 = np.array([0.0, 0.0, 0.0])
p1 = np.array([1.0, 0.0, 0.0])
p2 = np.array([1.0, 1.0, 0.0])
p3 = np.array([1.0, 1.0, 1.0])
dihedral = calc_dihedral(p0, p1, p2, p3)
```

### Rotation Matrices

```python
import numpy as np
from ttk.math.rotation import RotationMatrix, rigid_transform_3D

# Identity rotation (4x4 homogeneous matrix internally)
rm = RotationMatrix()

# Apply a pure translation to a single coordinate vector
rm.matrix[:3, 3] = [1.0, 2.0, 3.0]
shifted = rm.apply(np.array([0.0, 0.0, 0.0]))   # [1, 2, 3]

# Apply to a stack of coordinates (N x 3 array)
coords = np.random.rand(10, 3)
transformed = rm.apply(coords)

# SVD-based optimal rigid body alignment of two point clouds
# A and B must be 3 x N arrays (each column is one corresponding point)
A = np.random.rand(3, 5)
B = A + 0.1   # slightly displaced copy
transform_matrix = rigid_transform_3D(A, B)  # returns 4 x 4 homogeneous matrix

# Build a rotation from two spanning vectors
rm2 = RotationMatrix.from_bivec(
    np.array([1.0, 0.0, 0.0]),
    np.array([0.0, 1.0, 0.0]),
)
```

### Backbone Dihedral Angles (phi / psi)

```python
# Calculate phi/psi for every residue in a chain
chain = top.chains[0]
chain.calculate_psi_phi()

for res in chain.residues:
    phi = res.property_computed.get("phi")
    psi = res.property_computed.get("psi")
    if phi is not None and psi is not None:
        print(
            f"{res.name}{res.res_seq}: "
            f"phi={np.degrees(phi):.1f} deg  "
            f"psi={np.degrees(psi):.1f} deg"
        )
```

---

## Unit System

NeoTopology uses [pint](https://pint.readthedocs.io/) for physical units.
All atomic coordinates are stored in **nanometers** by default.

```python
from ttk import unit

# Create a quantity
pos = unit.Quantity([1.0, 2.0, 3.0], unit.nanometer)

# Convert units
pos_angstrom = pos.to("angstrom")
print(pos_angstrom.magnitude)   # [10.0, 20.0, 30.0]

# Strip units (get raw numpy array)
print(pos.magnitude)            # [1.0, 2.0, 3.0]
```

> Warning: Always attach units to coordinate arrays.
> Passing a plain list or numpy array as an atom position will raise an AssertionError.

---

## Command-Line Tools

### Fix a PDB File

Uses [PDBFixer](https://github.com/openmm/pdbfixer) + OpenMM to repair missing residues,
missing atoms, non-standard residues, and optionally add hydrogens.

```bash
python bin/fixpdb.py input.pdb output.pdb [--padding 1.0] [--addH True] [--pH 7.4]
```

### Align Structures with PyMOL

Requires a working PyMOL installation.

```bash
python bin/pymol_align.py -m mobile.pdb -t target.pdb [-o aligned] [-mode align]
```

| Mode | Use case |
|------|----------|
| `align` | High sequence similarity |
| `cealign` | Sequence-independent CE alignment |
| `super` | Sequence-independent superposition |

---

## Project Layout

```
NeoTopology/
├── src/ttk/
│   ├── __init__.py          # package entry-point (version, unit registry, Topology export)
│   ├── topology.py          # Topology class
│   ├── core/                # Atom, Bond, BondType, Chain, Residue, Entity, Element, UnitCell
│   ├── io/                  # PDB/PDBx parsing, PDB export, OpenMM/RDKit converters
│   ├── calculator/          # center of mass, RMSD, dimer interface
│   ├── math/                # rotation matrices, vector ops, rigid alignment
│   ├── data/                # chemical constants (amino acids, solvents, metals)
│   └── utils/               # sequence utilities
├── bin/                     # CLI tools
├── tests/                   # pytest test suite
└── pyproject.toml           # build config and tool settings
```

---

## Development

```bash
# Install with dev extras
pip install -e ".[dev]"

# Run tests
pytest tests/

# Run tests with coverage
pytest tests/ --cov=src/ttk --cov-report=html

# Lint
ruff check src/ tests/

# Format
ruff format src/ tests/

# Type-check
mypy src/
```

---

## CI / CD

| Workflow | Trigger | Jobs |
|----------|---------|------|
| **Tests and Code Quality** | Push / PR to `main`, `dev` | `test` (Python 3.11 / 3.12 / 3.13), `lint`, `type-check` |
| **Upload Python Package** | GitHub Release published | Build wheel + sdist, publish to PyPI |

Coverage reports are uploaded to [Codecov](https://codecov.io/).

---

## Contributing

1. Fork the repository and create a branch from `dev`.
2. Write tests for any new functionality.
3. Ensure `pytest`, `ruff`, and `mypy` all pass locally.
4. Open a Pull Request against the `dev` branch.

---

## License

Released under the [MIT License](LICENSE).
