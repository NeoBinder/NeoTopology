"""测试 ttk.topology 模块"""

import pytest
from ttk.topology import Topology, expand_symmetry
from ttk.core import element_from_symbol, Atom, Bond, Chain, Residue


class TestTopologyInit:
    """测试 Topology 初始化"""

    def test_topology_init(self):
        """测试创建空的 Topology"""
        topo = Topology()
        assert topo.chains == []
        assert topo.periodic_box is None
        assert topo.n_chains == 0


class TestTopologyProperties:
    """测试 Topology 属性"""

    def test_n_chains(self):
        """测试 n_chains 属性"""
        topo = Topology()
        assert topo.n_chains == 0

        topo.add_chain()
        assert topo.n_chains == 1

        topo.add_chain()
        assert topo.n_chains == 2

    def test_n_residues(self):
        """测试 n_residues 属性"""
        topo = Topology()
        chain = topo.add_chain()
        assert topo.n_residues == 0

        topo.add_residue("ALA", chain)
        assert topo.n_residues == 1

        topo.add_residue("GLY", chain)
        assert topo.n_residues == 2

    def test_residues_property(self):
        """测试 residues 属性"""
        topo = Topology()
        chain = topo.add_chain()

        res1 = topo.add_residue("ALA", chain)
        res2 = topo.add_residue("GLY", chain)

        residues = topo.residues
        assert len(residues) == 2
        assert res1 in residues
        assert res2 in residues

    def test_atoms_property(self):
        """测试 atoms 属性"""
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)

        atom1 = topo.add_atom("CA", element_from_symbol("C"), res, is_hetero=False)
        atom2 = topo.add_atom("N", element_from_symbol("N"), res, is_hetero=False)

        atoms = topo.atoms
        assert len(atoms) == 2
        assert atom1 in atoms
        assert atom2 in atoms


class TestTopologyAddOperations:
    """测试 Topology 添加操作"""

    def test_add_chain(self):
        """测试添加链"""
        topo = Topology()

        chain1 = topo.add_chain()
        assert chain1 is not None
        assert topo.n_chains == 1

        chain2 = topo.add_chain("B")
        assert chain2 is not None
        assert topo.n_chains == 2
        assert chain2.id == "B"

    def test_add_chain_auto_id(self):
        """测试自动生成链 ID"""
        topo = Topology()

        chain1 = topo.add_chain()
        assert chain1.id == "A"

        chain2 = topo.add_chain()
        assert chain2.id == "B"

        chain3 = topo.add_chain("Z")
        assert chain3.id == "Z"

    def test_add_residue(self):
        """测试添加残基"""
        topo = Topology()
        chain = topo.add_chain()

        res = topo.add_residue("ALA", chain)
        assert res is not None
        assert res.name == "ALA"
        assert res.chain == chain
        assert chain.residues[-1] == res

    def test_add_residue_with_seq(self):
        """测试添加带序号的残基"""
        topo = Topology()
        chain = topo.add_chain()

        res = topo.add_residue("ALA", chain, res_seq=100)
        assert res.res_seq == 100

    def test_add_atom(self):
        """测试添加原子"""
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)

        atom = topo.add_atom("CA", element_from_symbol("C"), res, is_hetero=False)

        assert atom is not None
        assert atom.name == "CA"
        assert atom.element.symbol == "C"
        assert atom.residue == res
        assert atom in res.atoms

    def test_add_atom_with_position(self):
        """测试添加带位置的原子"""
        from ttk import unit

        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)

        position = unit.Quantity([1.0, 2.0, 3.0], unit.nanometer)
        atom = topo.add_atom(
            "CA", element_from_symbol("C"), res, is_hetero=False, position=position
        )

        assert atom.position.magnitude.tolist() == position.magnitude.tolist()

    def test_add_bond(self):
        """测试添加键"""
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)

        atom1 = topo.add_atom("CA", element_from_symbol("C"), res, is_hetero=False)
        atom2 = topo.add_atom("N", element_from_symbol("N"), res, is_hetero=False)

        bond = topo.add_bond(atom1, atom2)

        assert bond is not None
        assert bond.atom1 == atom1
        assert bond.atom2 == atom2
        assert bond in atom1.bonds
        assert bond in atom2.bonds

    def test_add_bond_with_type(self):
        """测试添加带类型的键"""
        from ttk.core import BondTypeDouble

        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)

        atom1 = topo.add_atom("CA", element_from_symbol("C"), res, is_hetero=False)
        atom2 = topo.add_atom("C", element_from_symbol("C"), res, is_hetero=False)

        bond = topo.add_bond(atom1, atom2, bondtype=BondTypeDouble)

        assert bond.type == BondTypeDouble

    def test_add_topology(self):
        """测试合并两个 Topology"""
        topo1 = Topology()
        chain1 = topo1.add_chain()
        res1 = topo1.add_residue("ALA", chain1)
        atom1 = topo1.add_atom("CA", element_from_symbol("C"), res1, is_hetero=False)

        topo2 = Topology()
        chain2 = topo2.add_chain()
        res2 = topo2.add_residue("GLY", chain2)
        atom2 = topo2.add_atom("CA", element_from_symbol("C"), res2, is_hetero=False)

        # 合并前
        assert topo1.n_chains == 1
        assert topo1.n_residues == 1

        # 合并
        topo1.add_topology(topo2)

        # 合并后
        assert topo1.n_chains == 2
        assert topo1.n_residues == 2


class TestTopologyGetOperations:
    """测试 Topology 获取操作"""

    def test_get_atom_by_indices(self):
        """测试按索引获取原子"""
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)

        atom1 = topo.add_atom("CA", element_from_symbol("C"), res, is_hetero=False, index=0)
        atom2 = topo.add_atom("N", element_from_symbol("N"), res, is_hetero=False, index=1)

        atoms = topo.get_atom_by_indices([0, 1])
        assert len(atoms) == 2
        assert atom1 in atoms
        assert atom2 in atoms

    def test_get_chains(self):
        """测试获取链列表"""
        topo = Topology()
        chain1 = topo.add_chain()
        chain2 = topo.add_chain()

        chains = topo.get_chains()
        assert len(chains) == 2
        assert chain1 in chains
        assert chain2 in chains

    def test_get_residues(self):
        """测试获取残基生成器"""
        topo = Topology()
        chain = topo.add_chain()

        res1 = topo.add_residue("ALA", chain)
        res2 = topo.add_residue("GLY", chain)

        residues = list(topo.get_residues())
        assert len(residues) == 2
        assert res1 in residues
        assert res2 in residues

    def test_get_atoms(self):
        """测试获取原子生成器"""
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)

        atom1 = topo.add_atom("CA", element_from_symbol("C"), res, is_hetero=False)
        atom2 = topo.add_atom("N", element_from_symbol("N"), res, is_hetero=False)

        atoms = list(topo.get_atoms())
        assert len(atoms) == 2
        assert atom1 in atoms
        assert atom2 in atoms

    def test_get_atoms_with_index(self):
        """测试获取原子并创建索引"""
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)

        atom1 = topo.add_atom("CA", element_from_symbol("C"), res, is_hetero=False)
        atom2 = topo.add_atom("N", element_from_symbol("N"), res, is_hetero=False)

        # 获取原子并创建索引
        atoms = list(topo.get_atoms(create_index=True))

        assert len(atoms) == 2
        assert atom1.index == 0
        assert atom2.index == 1

    def test_get_bonds(self):
        """测试获取键生成器"""
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)

        atom1 = topo.add_atom("CA", element_from_symbol("C"), res, is_hetero=False)
        atom2 = topo.add_atom("N", element_from_symbol("N"), res, is_hetero=False)

        bond = topo.add_bond(atom1, atom2)

        bonds = list(topo.get_bonds())
        assert len(bonds) == 1
        assert bond in bonds

    def test_get_residues_by_name_single(self):
        """测试按名称获取残基（单个名称）"""
        topo = Topology()
        chain = topo.add_chain()

        res1 = topo.add_residue("ALA", chain)
        res2 = topo.add_residue("GLY", chain)

        residues = list(topo.get_residues_by_name("ALA"))
        assert len(residues) == 1
        assert res1 in residues
        assert res2 not in residues

    def test_get_residues_by_name_set(self):
        """测试按名称获取残基（集合）"""
        topo = Topology()
        chain = topo.add_chain()

        res1 = topo.add_residue("ALA", chain)
        res2 = topo.add_residue("GLY", chain)
        res3 = topo.add_residue("VAL", chain)

        residues = list(topo.get_residues_by_name(set(["ALA", "GLY"])))
        assert len(residues) == 2
        assert res1 in residues
        assert res2 in residues
        assert res3 not in residues

    def test_get_residues_by_name_list(self):
        """测试按名称获取残基（列表）"""
        topo = Topology()
        chain = topo.add_chain()

        res1 = topo.add_residue("ALA", chain)
        res2 = topo.add_residue("GLY", chain)

        residues = list(topo.get_residues_by_name(["ALA"]))
        assert len(residues) == 1
        assert res1 in residues


class TestTopologyDeleteOperations:
    """测试 Topology 删除操作"""

    def test_remove_water(self):
        """测试删除水分子"""
        topo = Topology()
        chain = topo.add_chain()

        topo.add_residue("HOH", chain)  # 水分子
        topo.add_residue("ALA", chain)  # 氨基酸

        assert topo.n_residues == 2

        result = topo.remove_water()

        assert result is True
        assert topo.n_residues == 1
        assert chain.residues[0].name == "ALA"

    def test_remove_solvent(self):
        """测试删除所有溶剂"""
        topo = Topology()
        chain = topo.add_chain()

        topo.add_residue("HOH", chain)  # 水分子
        topo.add_residue("NA", chain)  # 离子
        topo.add_residue("ALA", chain)  # 氨基酸

        assert topo.n_residues == 3

        result = topo.remove_solvent()

        assert result is True
        assert topo.n_residues == 1
        assert chain.residues[0].name == "ALA"


class TestTopologyModifyOperations:
    """测试 Topology 修改操作"""

    def test_create_index(self):
        """测试创建索引"""
        topo = Topology()
        chain = topo.add_chain()

        res1 = topo.add_residue("ALA", chain)
        atom1 = topo.add_atom("CA", element_from_symbol("C"), res1, is_hetero=False)
        atom2 = topo.add_atom("N", element_from_symbol("N"), res1, is_hetero=False)

        res2 = topo.add_residue("GLY", chain)
        atom3 = topo.add_atom("CA", element_from_symbol("C"), res2, is_hetero=False)

        topo.create_index()

        assert res1.index == 0
        assert res2.index == 1
        assert atom1.index == 0
        assert atom2.index == 1
        assert atom3.index == 2

    def test_create_index_with_start(self):
        """测试从指定索引开始创建索引"""
        topo = Topology()
        chain = topo.add_chain()

        res1 = topo.add_residue("ALA", chain)
        atom1 = topo.add_atom("CA", element_from_symbol("C"), res1, is_hetero=False)

        topo.create_index(atom_start_idx=10, res_start_idx=5)

        assert res1.index == 5
        assert atom1.index == 10


class TestExpandSymmetry:
    """测试 expand_symmetry 函数"""

    def test_expand_symmetry_empty(self):
        """测试展开空对称性"""
        topo = Topology()
        topo.symmetry = {}  # 空对称性

        result = expand_symmetry(topo)

        assert result.n_chains == 0

    def test_expand_symmetry_basic(self):
        """测试基本对称性展开（需要实现 symmetry 属性）"""
        # 注意：这个测试需要 Topology 类实现 symmetry 属性
        # 当前代码中可能没有实现，暂时跳过
        pytest.skip("需要 Topology.symmetry 属性实现")
