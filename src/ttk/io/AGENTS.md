# IO 模块

文件解析和导出，支持多种分子结构格式。

## 结构

| 文件 | 用途 |
|------|------|
| `topology_parser.py` | PDB/PDBx/OpenMM/rdkit 解析为 Topology |
| `topology_export.py` | Topology 导出为 PDB/OpenMM/几何字符串 |
| `PDBFile.py` | PDB 文件解析和写入（ATOM/HETATM/CONECT） |
| `PDBxParser.py` | PDBx（mmCIF）格式解析 |
| `PDBUtil.py` | PDB 头信息解析、对称矩阵、元数据提取 |
| `visualize.py` | py3Dmol 可视化、序列对齐展示 |
| `align.py` | MUSCLE 多序列对齐 |
| `residue_export.py` | 残基转 rdkit Mol 对象 |

## 主要函数

| 函数 | 文件 | 用途 |
|------|------|------|
| `topology_from_pdb(pdbpath)` | topology_parser.py | 从 PDB 文件加载 |
| `topology_from_pdb_content(content)` | topology_parser.py | 从 PDB 字符串加载 |
| `topology_from_openmmm(mm_topology)` | topology_parser.py | 从 OpenMM 对象加载 |
| `topology_from_rdkitmol(mol, res_name)` | topology_parser.py | 从 rdkit 分子加载 |
| `topology_to_pdb(ttk_topology, fname)` | topology_export.py | 导出为 PDB 文件 |
| `topology_to_openmm(ttk_topology)` | topology_export.py | 转换为 OpenMM 对象 |
| `res2mol(res, removeHs)` | residue_export.py | 残基转为 rdkit Mol |
| `visualise_topology(topology)` | visualize.py | py3Dmol 交互式可视化 |

## 支持格式

| 格式 | 解析 | 导出 | 说明 |
|------|------|------|------|
| PDB | ✅ | ✅ | 支持 gzip 压缩（.gz）、多模型、CONECT 键 |
| PDBx/mmCIF | ✅ | ❌ | 完整解析 atom_site、struct_conn |
| OpenMM | ✅ | ✅ | 自动单位转换（nm ↔ Å） |
| rdkit | ✅ | ✅ | 用于配体编辑、SMILES 模板 |
| 几何字符串 | ❌ | ✅ | 用于量子化学输入 |

## 使用模式

```python
from ttk.io import topology_parser, topology_export, PDBFile

# 从不同来源加载
top = topology_parser.topology_from_pdb("protein.pdb")
top = topology_parser.topology_from_pdb_content(pdb_string)
top = topology_parser.topology_from_openmmm(openmm_top)
top = topology_parser.topology_from_rdkitmol(rdkit_mol, "LIG")

# 导出 PDB
content = PDBFile().to_content(top)  # 获取字符串
topology_export.topology_to_pdb(top, "output.pdb")  # 直接写入文件

# 转换为 OpenMM
omm_top, positions = topology_export.topology_to_openmm(top)

# 残基转 rdkit（用于 SMILES 模板）
res = top.residues[5]
mol = residue_export.res2mol(res, template="CCO")  # 带模板
```

## 注意事项

**解析细节**：
- PDB 坐标自动转换为纳米（Å → nm）
- PDB 解析器支持 `config={'extraParticleIdentifier': 'EP'}` 处理虚拟粒子
- OpenMM 导入时需手动转换单位：`pos.to("nm")`

**导出细节**：
- `PDBFile.to_content()` 自动生成 CRYST1（存在周期盒）
- 标准残基使用 `ATOM` 记录，非标准使用 `HETATM`
- 仅导出非标准残基的 `CONECT` 键

**可视化**：
- `visualise_topology(show_interface=True)` 自动高亮链间接口
- 序列对齐需安装 MUSCLE 可执行文件

**单位系统**：
- 所有位置默认使用纳米（`ttk.unit.nanometer`）
- 导出时自动转为埃（PDB 标准要求）
