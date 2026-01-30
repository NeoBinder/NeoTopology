# Core 模块

核心数据模型层，定义拓扑层级结构的基础类。

## 结构

| 文件 | 内容 |
|------|------|
| entity.py | Entity 基类，提供原子集合的通用操作 |
| atom.py | Atom 原子类 |
| residue.py | Residue 残基类（继承 Entity） |
| chain.py | Chain 链类（继承 Entity） |
| bond.py | Bond 键类、BondType 键类型枚举 |
| element.py | 元素符号解析（基于 mendeleev） |
| unit_cell.py | 晶胞参数转换 |
| utils.py | 原子筛选、单位转换工具 |

## 核心类

| 类 | 文件 | 职责 |
|------|------|------|
| Entity | entity.py | 顶层基类，管理原子集合和坐标 |
| Atom | atom.py | 原子对象，包含元素、位置、键连接 |
| Residue | residue.py | 残基对象，包含原子列表和残基属性 |
| Chain | chain.py | 链对象，包含残基列表和二级结构计算 |
| Bond | bond.py | 键对象，连接两个原子 |
| BondType | bond.py | 键类型枚举（Single/Double/Triple/Aromatic） |
| UnitCell | unit_cell.py | 晶胞参数，支持 OpenMM 互转 |

## 关键属性/方法

| 名称 | 类 | 用途 |
|------|------|------|
| positions | Entity | 坐标数组（带单位，setter 自动分配到原子） |
| heavy_atoms | Entity | 重原子列表（排除氢） |
| select_by_dist() | Entity | 按距离筛选原子 |
| name/element/residue | Atom | 原子名称、元素、所属残基 |
| is_backbone/is_sidechain | Atom | 判断主链/侧链原子 |
| name/res_seq/chain | Residue | 残基名称、序号、所属链 |
| is_protein/is_water/is_solvent | Residue | 残基类型判断 |
| calculate_psi_phi() | Chain | 计算 phi/psi 二面角 |
| to_topology() | Chain | 链转换为独立 Topology |
| connect() | Bond | 获取键的另一个原子 |
| from_parameter()/from_openmm() | UnitCell | 从参数/OpenMM 创建晶胞 |
| element_from_symbol() | - | 元素符号解析 |

## 约定

**层级继承**：
- Entity → Topology/Chain/Residue（单继承）
- 每层都有 `atoms` 属性，支持层级遍历

**坐标单位**：
- Atom.position 必须是 `ttk.unit.Quantity` 类型
- Entity.positions setter 自动检测单位，无单位时警告并默认为纳米

**命名规范**：
- Atom 使用 `clean_name`（去除空格）用于快速查找
- Residue 支持 `residue["CA"]` 语法访问原子

## 注意事项

**性能陷阱**：
- `atom.n_bonds` 已禁用（计算开销大），改用 `atom.bonds` 列表
- Chain 的 `calculate_psi_phi()` 需要完整骨架原子（N/CA/C）

**残基判断**：
- `is_protein()` 基于 `_PROTEIN_RESIDUES` 集合判断
- `is_water()` 基于 `_WATER_RESIDUES` 集合判断
- `is_solvent()` 基于 `_SIMPLE_SOLVENT` 集合判断

**键类型转换**：
- BondType 支持与 rdkit.BondType 互转（`to_rdkit()`/`from_rdkit_bond()`）
- 未知键类型使用 `BondTypeUnknown`

**元素缓存**：
- `element_from_symbol()` 会缓存解析结果到 `symbol_element_dict`
