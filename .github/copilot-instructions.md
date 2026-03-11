# Copilot Instructions for NeoTopology

## Project Overview

NeoTopology is NeoBinder's all-in-one SDK for protein topology, designed for protein design, OpenMM preparation, and research. It provides PDB/PDBx/rdkit file I/O, rotation and alignment math, and protein/sequence/ligand calculations.

**Key identifiers:**
- PyPI package name: `neotopology`
- Python import name: `ttk`
- Import example: `from ttk import Topology`

## Repository Structure

```
NeoTopology/
├── src/ttk/                 # Main package (imported as `ttk`)
│   ├── __init__.py          # Version definition and core exports
│   ├── topology.py          # Top-level Topology class (chains/residues/atoms/bonds)
│   ├── core/                # Core data models: Atom, Bond, Chain, Residue, Entity
│   ├── io/                  # PDB/PDBx parsing, export, visualization
│   ├── calculator/          # RMSD, center of mass, dimer interface calculations
│   ├── math/                # Rotation matrix, vector operations
│   ├── data/                # Chemical data constants (solvents, amino acids, metals)
│   └── utils/               # Sequence processing utilities
├── bin/                     # CLI tools (fixpdb.py, pymol_align.py)
├── tests/                   # pytest test suite
├── pyproject.toml           # Build config (setuptools, ruff, mypy)
└── .github/workflows/       # CI (test.yml) and PyPI publish (push_pypi.yaml)
```

## Key File Locations

| Task | Location | Notes |
|------|----------|-------|
| Main entry point | `src/ttk/__init__.py` | Version and core exports |
| Topology class | `src/ttk/topology.py` | chains/residues/atoms/bonds |
| PDB parsing | `src/ttk/io/topology_parser.py` | `topology_from_pdb()`, `topology_from_pdb_content()` |
| PDB export | `src/ttk/io/topology_export.py` | `topology_to_pdb()` |
| OpenMM conversion | `src/ttk/io/topology_parser.py` | `topology_from_openmmm()` |
| Atom/Bond/Chain/Residue | `src/ttk/core/` | Entity base class hierarchy |
| Element data | `src/ttk/core/element.py` | `element_from_symbol()` |
| Rotation matrix | `src/ttk/math/rotation.py` | `RotationMatrix` class |
| Center of mass | `src/ttk/calculator/__init__.py` | `get_center_of_mass()` |
| Solvent detection | `src/ttk/core/residue.py` | `residue.is_solvent`, `residue.is_water` |
| Dimer interface | `src/ttk/calculator/dimer.py` | `find_interface()` |
| Chemical constants | `src/ttk/data/` | `_AMINO_ACID_CODES`, `_SOLVENT_TYPES`, `_WATER_RESIDUES`, `_Metal` |

## Development Setup

```bash
# Recommended: uv package manager
uv venv -p 3.11
source ./.venv/bin/activate
uv sync

# Alternative: pip
pip install -e ./
```

## How to Build, Lint, and Test

```bash
# Run all tests
pytest tests/

# Run a specific test file
pytest tests/test_topology.py

# Lint with ruff
ruff check src/

# Type check with mypy
mypy src/

# Build the package
python -m build
```

## Conventions

### Units (Critical)
- All positions/distances use **nanometers** (`ttk.unit.nanometer`) by default.
- Use `pint.UnitRegistry` quantities (`ttk.unit.Quantity`) for all coordinates — **never** raw numpy arrays.
- When importing from OpenMM, convert units explicitly: `pos.to("nm")`.

### Object Hierarchy
- `Entity` is the base class for `Topology`, `Chain`, `Residue`, and `Atom`.
- Every level exposes an `atoms` property for hierarchical access.
- Filtering by name/index/type is supported at every level.

### Indexing
- Never manually maintain atom/bond indices; always call `topology.create_index()` after structural changes.

### Element Symbols
- Never hardcode element symbols in code; always retrieve them via `ttk.core.element.element_from_symbol()`.

### Data Constants
- All chemical constants (amino acid codes, solvent types, metal sets) live in `src/ttk/data/`.
- Use these constants rather than creating new string literals or sets elsewhere.

## Anti-Patterns to Avoid

- ❌ Using raw `numpy` arrays for coordinates — use `ttk.unit.Quantity` instead.
- ❌ Hardcoding atom element symbols — use `element_from_symbol()` from `ttk.core.element`.
- ❌ Manually managing atom/bond indices — call `topology.create_index()`.
- ❌ Using `setup.py` — the project uses `pyproject.toml`.

## Important Notes

- `topology.remove_water()` removes only water molecules; `topology.remove_solvent()` removes all solvents (water + ions, etc.).
- PDB parsing auto-detects gzip compression (`.gz` suffix).
- `atom.n_bonds` is disabled (commented out) due to performance cost; avoid re-enabling without profiling.
- CI uses Python 3.9; local development targets Python 3.11.
- PyPI publish is triggered only on `release: published` events.
