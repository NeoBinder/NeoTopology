"""测试 ttk.core 模块"""

import numpy as np
import pytest

from ttk import unit
from ttk.core import (
    Atom,
    Bond,
    BondType,
    BondTypeUnknown,
    BondTypeSingle,
    BondTypeDouble,
    BondTypeTriple,
    BondTypeAromatic,
    BondTypeAmide,
    Chain,
    Entity,
    Residue,
    element_from_symbol,
    from_bond_float,
    from_bond_name,
)


class TestElement:
    """测试元素相关函数"""

    def test_element_from_symbol(self):
        """测试从符号获取元素"""
        carbon = element_from_symbol("C")
        assert carbon.symbol == "C"
        assert carbon.mass > 0

    def test_element_from_symbol_unknown(self):
        """测试未知元素符号"""
        unknown = element_from_symbol("X")
        assert unknown.symbol == "X"
        assert unknown.mass == 0.0


class TestAtom:
    """测试 Atom 类"""

    def test_atom_init(self):
        """测试原子初始化"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)

        assert atom.name == "CA"
        assert atom.element.symbol == "C"
        assert atom.residue == residue
        assert atom.is_hetero is False
        assert atom.bonds_dict == {}
        assert atom._position is None

    def test_atom_position(self):
        """测试原子位置"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)

        position = unit.Quantity([1.0, 2.0, 3.0], unit.nanometer)
        atom.position = position

        assert np.allclose(atom.position.magnitude, position.magnitude)

    def test_atom_position_without_unit(self):
        """测试不带单位的原子位置"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)

        with pytest.raises(AssertionError):
            atom.position = [1.0, 2.0, 3.0]

    def test_atom_index(self):
        """测试原子索引"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)

        assert atom.index == 0

        atom.index = 5
        assert atom.index == 5

    def test_atom_bonds(self):
        """测试原子键"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom1 = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        atom2 = Atom("N", element_from_symbol("N"), residue, is_hetero=False)

        bond = Bond(atom1, atom2)
        atom1.bonds_dict[bond.connection_hash] = bond

        assert len(atom1.bonds) == 1
        assert bond in atom1.bonds

    def test_atom_symbol(self):
        """测试原子符号属性"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)

        assert atom.symbol == "C"

    def test_atom_clean_name(self):
        """测试原子清理名称"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom = Atom("CA ", element_from_symbol("C"), residue, is_hetero=False)

        assert atom.clean_name == "CA"

    def test_atom_is_backbone(self):
        """测试是否为主链原子"""
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")
        residue.name = "ALA"  # 确保是蛋白质残基

        ca_atom = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        cb_atom = Atom("CB", element_from_symbol("C"), residue, is_hetero=False)

        assert ca_atom.is_backbone is True
        assert cb_atom.is_backbone is False

    def test_atom_is_sidechain(self):
        """测试是否为侧链原子"""
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")
        residue.name = "ALA"  # 确保是蛋白质残基

        cb_atom = Atom("CB", element_from_symbol("C"), residue, is_hetero=False)
        n_atom = Atom("N", element_from_symbol("N"), residue, is_hetero=False)

        assert cb_atom.is_sidechain is True
        assert n_atom.is_sidechain is False

    def test_atom_occupancy(self):
        """测试原子占据率"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)

        assert atom.occupancy == 1.00

    def test_atom_str_and_repr(self):
        """测试原子字符串表示"""
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")
        atom = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)

        str_repr = str(atom)
        assert "ALA" in str_repr
        assert "CA" in str_repr

        repr_str = repr(atom)
        assert "Atom:" in repr_str


class TestBondType:
    """测试 BondType 类"""

    def test_bondtype_single(self):
        """测试单键"""
        single = BondType.Single()
        assert single.name == "Single"
        assert float(single) == 1.0

    def test_bondtype_double(self):
        """测试双键"""
        double = BondType.Double()
        assert double.name == "Double"
        assert float(double) == 2.0

    def test_bondtype_triple(self):
        """测试三键"""
        triple = BondType.Triple()
        assert triple.name == "Triple"
        assert float(triple) == 3.0

    def test_bondtype_aromatic(self):
        """测试芳香键"""
        aromatic = BondType.Aromatic()
        assert aromatic.name == "Aromatic"
        assert float(aromatic) == 1.5

    def test_bondtype_amide(self):
        """测试酰胺键"""
        amide = BondType.Amide()
        assert amide.name == "Amide"
        assert float(amide) == 1.25

    def test_bondtype_unknown(self):
        """测试未知键类型"""
        unknown = BondType.Unknown()
        assert unknown.name == "Unknown"
        assert float(unknown) == 0.0

    def test_bondtype_eq(self):
        """测试键类型相等"""
        single1 = BondType.Single()
        single2 = BondType.Single()

        assert single1 == single2

    def test_bondtype_hash(self):
        """测试键类型哈希"""
        single = BondType.Single()
        double = BondType.Double()

        assert hash(single) != hash(double)


class TestBond:
    """测试 Bond 类"""

    def test_bond_init(self):
        """测试键初始化"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom1 = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        atom2 = Atom("N", element_from_symbol("N"), residue, is_hetero=False)

        bond = Bond(atom1, atom2)

        assert bond.atom1 == atom1
        assert bond.atom2 == atom2
        assert bond.type == BondTypeUnknown

    def test_bond_with_type(self):
        """测试带类型的键"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom1 = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        atom2 = Atom("C", element_from_symbol("C"), residue, is_hetero=False)

        bond = Bond(atom1, atom2, bondtype=BondTypeDouble)

        assert bond.type == BondTypeDouble

    def test_bond_connection_hash(self):
        """测试键连接哈希"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom1 = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        atom2 = Atom("N", element_from_symbol("N"), residue, is_hetero=False)

        bond1 = Bond(atom1, atom2)
        bond2 = Bond(atom2, atom1)

        assert bond1.connection_hash == bond2.connection_hash

    def test_bond_eq(self):
        """测试键相等"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom1 = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        atom2 = Atom("N", element_from_symbol("N"), residue, is_hetero=False)

        bond1 = Bond(atom1, atom2)
        bond2 = Bond(atom2, atom1)

        assert bond1 == bond2

    def test_bond_connect(self):
        """测试键连接"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom1 = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        atom2 = Atom("N", element_from_symbol("N"), residue, is_hetero=False)

        bond = Bond(atom1, atom2)

        assert bond.connect(atom1) == atom2
        assert bond.connect(atom2) == atom1

    def test_bond_connect_invalid(self):
        """测试无效的键连接"""
        residue = Residue("ALA", 0, Chain(0, None), "1")
        atom1 = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        atom2 = Atom("N", element_from_symbol("N"), residue, is_hetero=False)
        atom3 = Atom("C", element_from_symbol("C"), residue, is_hetero=False)

        bond = Bond(atom1, atom2)

        with pytest.raises(Exception):
            bond.connect(atom3)


class TestResidue:
    """测试 Residue 类"""

    def test_residue_init(self):
        """测试残基初始化"""
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")

        assert residue.name == "ALA"
        assert residue.index == 0
        assert residue.chain == chain
        assert residue.res_seq == "1"
        assert residue.segment_id == ""
        assert residue.atoms == []

    def test_residue_str_and_repr(self):
        """测试残基字符串表示"""
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")

        str_repr = str(residue)
        assert "ALA" in str_repr
        assert "1" in str_repr

        repr_str = repr(residue)
        assert "Residue:" in repr_str

    def test_residue_add_atom(self):
        """测试残基添加原子"""
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")
        atom = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)

        residue.add_atom(atom)

        assert len(residue.atoms) == 1
        assert atom in residue.atoms

    def test_residue_getitem(self):
        """测试残基索引访问"""
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")
        atom = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        residue.add_atom(atom)

        retrieved_atom = residue["CA"]
        assert retrieved_atom == atom

    def test_residue_getitem_invalid(self):
        """测试无效的残基索引访问"""
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")

        with pytest.raises(KeyError):
            residue["CA"]

    def test_residue_get(self):
        """测试残基 get 方法"""
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")
        atom = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        residue.add_atom(atom)

        assert residue.get("CA") == atom
        assert residue.get("CB") is None
        assert residue.get("CB", "default") == "default"

    def test_residue_atoms_by_name(self):
        """测试按名称获取原子"""
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")
        atom1 = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        atom2 = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        atom3 = Atom("N", element_from_symbol("N"), residue, is_hetero=False)
        residue.add_atom(atom1)
        residue.add_atom(atom2)
        residue.add_atom(atom3)

        atoms = list(residue.atoms_by_name("CA"))
        assert len(atoms) == 2
        assert atom1 in atoms
        assert atom2 in atoms

    def test_residue_get_atom_indices(self):
        """测试获取原子索引"""
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")
        atom1 = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        atom2 = Atom("N", element_from_symbol("N"), residue, is_hetero=False)
        residue.add_atom(atom1)
        residue.add_atom(atom2)

        atom1.index = 0
        atom2.index = 1

        indices = residue.get_atom_indices()
        assert indices == [0, 1]


class TestEntity:
    """测试 Entity 类"""

    def test_entity_init(self):
        """测试实体初始化"""
        entity = Entity()

        assert entity.atoms == []

    def test_entity_positions(self):
        """测试实体位置"""
        entity = Entity()
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")
        atom1 = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        atom2 = Atom("N", element_from_symbol("N"), residue, is_hetero=False)

        position1 = unit.Quantity([1.0, 2.0, 3.0], unit.nanometer)
        position2 = unit.Quantity([4.0, 5.0, 6.0], unit.nanometer)

        atom1.position = position1
        atom2.position = position2

        entity.atoms = [atom1, atom2]

        positions = entity.positions
        assert positions.shape == (2, 3)
        assert positions[0][0] == position1[0]
        assert positions[1][0] == position2[0]

    def test_entity_positions_setter(self):
        """测试实体位置设置器"""
        entity = Entity()
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")
        atom1 = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        atom2 = Atom("N", element_from_symbol("N"), residue, is_hetero=False)

        entity.atoms = [atom1, atom2]

        positions = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        entity.positions = positions

        assert np.allclose(atom1.position.magnitude, positions[0])
        assert np.allclose(atom2.position.magnitude, positions[1])

    def test_entity_heavy_atoms(self):
        """测试重原子"""
        entity = Entity()
        chain = Chain(0, None)
        residue = Residue("ALA", 0, chain, "1")
        c_atom = Atom("CA", element_from_symbol("C"), residue, is_hetero=False)
        h_atom = Atom("H", element_from_symbol("H"), residue, is_hetero=False)

        entity.atoms = [c_atom, h_atom]

        heavy_atoms = entity.heavy_atoms
        assert len(heavy_atoms) == 1
        assert c_atom in heavy_atoms
        assert h_atom not in heavy_atoms


class TestBondUtilityFunctions:
    """测试键工具函数"""

    def test_from_bond_float(self):
        """测试从浮点数创建键类型"""
        assert from_bond_float(1.0) == BondTypeSingle
        assert from_bond_float(2.0) == BondTypeDouble
        assert from_bond_float(3.0) == BondTypeTriple
        assert from_bond_float(1.5) == BondTypeAromatic
        assert from_bond_float(1.25) == BondTypeAmide

    def test_from_bond_name(self):
        """测试从名称创建键类型"""
        assert from_bond_name("covale") == BondTypeSingle


class TestChain:
    """测试 Chain 类（如果有必要，可以添加更多测试）"""

    def test_chain_init(self):
        """测试链初始化"""
        from ttk.topology import Topology

        topo = Topology()
        chain = Chain(0, topo, "A")

        assert chain.index == 0
        assert chain.topology == topo
        assert chain.id == "A"
        assert chain.residues == []
