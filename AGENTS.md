# NeoTopology 知识库

**生成时间**: 2026-01-29
**项目**: NeoBinder 蛋白质拓扑工具包
**核心栈**: Python 3.7+, setuptools, numpy, rdkit, mendeleev, pint

## 项目概述

NeoTopology 是专为蛋白质设计、OpenMM 准备和研究打造的全合一 SDK。提供 PDB/PDBx/rdkit 文件 IO、旋转对齐等数学计算、以及蛋白质/序列/配体计算功能。

## 结构

```
NeoTopology/
├── src/ttk/           # 主包
│   ├── core/          # 核心数据模型（原子/键/链/残基）
│   ├── io/            # PDB/PDBx 解析、导出、可视化
│   ├── calculator/    # RMSD、质心计算
│   ├── math/          # 旋转矩阵、向量运算
│   ├── data/          # 溶剂、氨基酸、金属元素常量
│   └── utils/         # 序列处理工具
├── bin/               # 命令行工具（修复 PDB、对齐）
├── pyproject.toml     # 现代构建配置
└── .github/workflows/ # CI/CD（发布到 PyPI）
```

## 寻找位置

| 任务 | 位置 | 说明 |
|------|------|------|
| 主入口 | `src/ttk/__init__.py` | 版本定义、核心导出 |
| Topology 核心类 | `src/ttk/topology.py` | 顶层拓扑结构，chains/residues/atoms/bonds |
| PDB 解析 | `src/ttk/io/topology_parser.py` | `topology_from_pdb()`, `topology_from_pdb_content()` |
| PDB 导出 | `src/ttk/io/topology_export.py` | `topology_to_pdb()` |
| OpenMM 转换 | `src/ttk/io/topology_parser.py` | `topology_from_openmmm()` |
| 原子/键/链/残基 | `src/ttk/core/` | Atom, Bond, Chain, Residue, Entity 类 |
| 元素数据 | `src/ttk/core/element.py` | `element_from_symbol()` |
| 旋转矩阵 | `src/ttk/math/rotation.py` | RotationMatrix 类 |
| 质心计算 | `src/ttk/calculator/__init__.py` | `get_center_of_mass()` |
| 溶剂判断 | `src/ttk/core/residue.py` | `residue.is_solvent`, `residue.is_water` |
| 二聚体接口 | `src/ttk/calculator/dimer.py` | `find_interface()` |

## 约定（偏离标准）

**单位系统**：
- 位置默认单位：纳米（nanometer, `ttk.unit.nanometer`）
- 使用 `pint.UnitRegistry` 统一处理物理单位
- 所有涉及距离/坐标的计算必须带单位

**包名差异**：
- PyPI 发布名：`neotopology`
- Python 包名：`ttk`
- 导入语句：`from ttk import Topology`

**现代工具链**：
- 推荐使用 `uv` 包管理器（比 pip 快）
- 构建配置：`pyproject.toml`（非 setup.py）
- 版本动态读取：`[tool.setuptools.dynamic] version = {attr = "ttk.__version__}`

## 反模式（本项目）

- ❌ 不要直接使用 numpy 数组表示坐标，必须使用 `ttk.unit.Quantity`
- ❌ 不要在代码中硬编码原子元素符号，应从 `ttk.core.element` 获取
- ❌ 不要手动维护原子/键索引，应调用 `topology.create_index()`

## 独特样式

**数据常量集中管理**：
- `ttk.data` 模块包含大量化学数据常量（1283 行）
- 氨基酸代码映射：`_AMINO_ACID_CODES`（PDB 3字母 → 1字母）
- 溶剂类型：`_SOLVENT_TYPES`, `_WATER_RESIDUES`
- 金属元素集合：`_Metal`

**层级结构设计**：
- Entity → Topology/Chain/Residue/Atom（继承关系）
- 每层都有 `atoms` 属性，提供层级访问
- 支持按名称/索引/类型筛选

## 命令

```bash
# 开发环境（推荐使用 uv）
uv venv -p 3.11
source ./.venv/bin/activate
uv sync

# 标准安装
pip install -e ./

# 构建
python -m build

# PDB 文件修复
python bin/fixpdb.py input.pdb output.pdb

# PyMOL 对齐
python bin/pymol_align.py ref.pdb mobile.pdb
```

## 注意事项

**单位转换陷阱**：
- 从 OpenMM 导入时需要手动转换单位：`pos.to("nm")`
- 未检测到单位时，`Entity.positions` 默认使用纳米单位但会警告

**性能提示**：
- `atom.n_bonds` 属性已禁用（注释代码），计算性能开销大
- 建议使用 `topology.create_index()` 后再访问索引属性

**溶剂清理**：
- `topology.remove_water()` 仅删除水分子
- `topology.remove_solvent()` 删除所有溶剂（包括离子等）

**PDB 解析细节**：
- 支持直接从文件路径、文件内容、OpenMM topology、rdkit molecule 解析
- PDB 文件会自动检测 gzip 压缩（`.gz` 后缀）

**CI/CD**：
- 仅在 `release: published` 时触发发布到 PyPI
- CI 使用 Python 3.9，开发推荐 Python 3.11（版本不一致）
- 使用 `python -m build` 而非 `setup.py sdist bdist_wheel`
