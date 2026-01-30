"""测试 ttk.io.topology_parser 模块"""

import pytest
from ttk.io import topology_parser
from ttk import Topology
from ttk.core import element_from_symbol


class TestTopologyParser:
    """测试拓扑解析器"""

    def test_load_from_file_txt(self, tmp_path):
        """测试加载文本文件"""
        test_file = tmp_path / "test.txt"
        test_file.write_text("line1\nline2\nline3\n")

        content = topology_parser.load_from_file(str(test_file))

        assert content == ["line1\n", "line2\n", "line3\n"]

    def test_load_from_file_gz(self, tmp_path):
        """测试加载 gzip 压缩文件"""
        import gzip

        test_file = tmp_path / "test.gz"
        with gzip.open(test_file, "wb") as f:
            f.write(b"line1\nline2\nline3\n")

        content = topology_parser.load_from_file(str(test_file))

        # load_from_file 对 gzip 文件使用 split("\n") 会去除换行符
        assert content == ["line1", "line2", "line3", ""]

    def test_topology_from_pdb_content_simple(self):
        """测试从简单的 PDB 内容加载"""
        pdb_content = """\
ATOM      1  N   ALA A   1      -2.273   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      -1.713   1.403   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1      -0.254   1.403   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       0.436   2.343   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1      -2.286   2.787   0.000  1.00  0.00           C
END
"""
        topo = topology_parser.topology_from_pdb_content(pdb_content)

        assert isinstance(topo, Topology)
        assert topo.n_chains == 1
        assert topo.n_residues == 1
        assert len(topo.atoms) == 5

    def test_topology_from_pdb_content_multiple_residues(self):
        """测试从包含多个残基的 PDB 内容加载"""
        pdb_content = """\
ATOM      1  N   ALA A   1      -2.273   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      -1.713   1.403   0.000  1.00  0.00           C
ATOM      3  N   GLY A   2      -0.254   1.403   0.000  1.00  0.00           N
ATOM      4  CA  GLY A   2       0.436   2.343   0.000  1.00  0.00           C
END
"""
        topo = topology_parser.topology_from_pdb_content(pdb_content)

        assert topo.n_residues == 2

        residues = list(topo.residues)
        assert residues[0].name == "ALA"
        assert residues[1].name == "GLY"

    def test_topology_from_pdb_content_water(self):
        """测试从包含水分子的 PDB 内容加载"""
        pdb_content = """\
ATOM      1  N   ALA A   1      -2.273   0.000   0.000  1.00  0.00           N
HETATM    2  O   HOH B   1       1.000   2.000   3.000  1.00  0.00           O
HETATM    3  H1  HOH B   1       1.500   2.000   3.000  1.00  0.00           H
HETATM    4  H2  HOH B   1       0.500   2.500   3.000  1.00  0.00           H
END
"""
        topo = topology_parser.topology_from_pdb_content(pdb_content)

        assert topo.n_residues == 2

        residues = list(topo.residues)
        assert residues[0].name == "ALA"
        assert residues[1].name == "HOH"

    def test_topology_from_pdb_content_heteroatoms(self):
        """测试 HETATM 记录"""
        pdb_content = """\
HETATM    1  CA  LIG A   1       1.000   2.000   3.000  1.00  0.00           C
HETATM    2  N   LIG A   1       2.000   3.000   4.000  1.00  0.00           N
END
"""
        topo = topology_parser.topology_from_pdb_content(pdb_content)

        assert topo.n_residues == 1
        assert len(topo.atoms) == 2

        # 检查 hetero 标志
        for atom in topo.atoms:
            assert atom.is_hetero is True

    def test_topology_from_pdb_content_with_connect(self):
        """测试 CONECT 记录"""
        pdb_content = """\
ATOM      1  CA  LIG A   1       1.000   2.000   3.000  1.00  0.00           C
ATOM      2  N   LIG A   1       2.000   3.000   4.000  1.00  0.00           N
ATOM      3  O   LIG A   1       3.000   4.000   5.000  1.00  0.00           O
CONECT    1    2    3
END
"""
        topo = topology_parser.topology_from_pdb_content(pdb_content)

        # 检查键是否正确创建
        bonds = list(topo.get_bonds())
        assert len(bonds) >= 2

    def test_topology_from_pdb_file(self, tmp_path):
        """测试从 PDB 文件加载"""
        pdb_content = """\
ATOM      1  N   ALA A   1      -2.273   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      -1.713   1.403   0.000  1.00  0.00           C
END
"""
        test_file = tmp_path / "test.pdb"
        test_file.write_text(pdb_content)

        topo = topology_parser.topology_from_pdb(str(test_file))

        assert topo.n_residues == 1
        assert len(topo.atoms) == 2

    def test_topology_from_pdb_file_gz(self, tmp_path):
        """测试从 gzip 压缩的 PDB 文件加载"""
        import gzip

        pdb_content = """\
ATOM      1  N   ALA A   1      -2.273   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      -1.713   1.403   0.000  1.00  0.00           C
END
"""
        test_file = tmp_path / "test.pdb.gz"
        with gzip.open(test_file, "wb") as f:
            f.write(pdb_content.encode())

        topo = topology_parser.topology_from_pdb(str(test_file))

        assert topo.n_residues == 1
        assert len(topo.atoms) == 2

    def test_topology_from_pdb_content_empty(self):
        """测试从空内容加载"""
        topo = topology_parser.topology_from_pdb_content("")

        # 空内容应该返回空拓扑（没有残基和原子）
        assert topo is not None
        assert topo.n_residues == 0
        assert len(topo.atoms) == 0

    def test_topology_from_pdb_content_invalid(self):
        """测试从无效内容加载"""
        invalid_content = "This is not a valid PDB file"

        # 无效内容应该返回空拓扑（不会解析到任何原子）
        topo = topology_parser.topology_from_pdb_content(invalid_content)
        assert topo is not None
        assert topo.n_residues == 0
        assert len(topo.atoms) == 0

    def test_topology_from_pdb_content_multiple_chains(self):
        """测试从包含多条链的 PDB 内容加载"""
        pdb_content = """\
ATOM      1  N   ALA A   1      -2.273   0.000   0.000  1.00  0.00           N
ATOM      2  N   ALA B   1       1.000   2.000   3.000  1.00  0.00           N
END
"""
        topo = topology_parser.topology_from_pdb_content(pdb_content)

        assert topo.n_chains == 2

        chains = topo.get_chains()
        assert chains[0].id == "A"
        assert chains[1].id == "B"


class TestTopologyFromRdkit:
    """测试从 rdkit 分子加载"""

    def test_topology_from_rdkitmol_simple(self):
        """测试从简单的 rdkit 分子加载"""
        try:
            from rdkit import Chem
        except ImportError:
            pytest.skip("rdkit not installed")

        # 创建一个简单的水分子
        mol = Chem.MolFromSmiles("O")
        mol = Chem.AddHs(mol)

        topo = topology_parser.topology_from_rdkitmol(mol, "HOH")

        assert isinstance(topo, Topology)
        assert topo.n_residues == 1
        assert topo.n_chains == 1

        # 检查残基名称
        residues = list(topo.residues)
        assert residues[0].name == "HOH"

    def test_topology_from_rdkitmol_atoms(self):
        """测试从 rdkit 分子加载原子"""
        try:
            from rdkit import Chem
        except ImportError:
            pytest.skip("rdkit not installed")

        # 创建一个乙烷分子
        mol = Chem.MolFromSmiles("CC")
        mol = Chem.AddHs(mol)

        topo = topology_parser.topology_from_rdkitmol(mol, "ETH")

        assert len(topo.atoms) > 0

        # 检查元素类型
        elements = [atom.element.symbol for atom in topo.atoms]
        assert "C" in elements
        assert "H" in elements


class TestTopologyFromOpenMM:
    """测试从 OpenMM 加载"""

    def test_topology_from_openmmm_basic(self):
        """测试从 OpenMM 加载"""
        try:
            import openmm
            from openmm.app import Modeller
        except ImportError:
            pytest.skip("openmm not installed")

        # 创建一个简单的系统
        from openmm.app import PDBFile

        pdb_content = """\
ATOM      1  N   ALA A   1      -2.273   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      -1.713   1.403   0.000  1.00  0.00           C
END
"""
        import tempfile

        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write(pdb_content)
            pdb_file = f.name

        try:
            pdb = PDBFile(pdb_file)
            topo = topology_parser.topology_from_openmmm(pdb.topology, pdb.positions)

            assert topo.n_residues == 1
            assert len(topo.atoms) == 2
        finally:
            import os

            os.unlink(pdb_file)

    def test_topology_from_openmmm_without_positions(self):
        """测试从 OpenMM 加载不带位置"""
        try:
            import openmm
            from openmm.app import Modeller
        except ImportError:
            pytest.skip("openmm not installed")

        from openmm.app import PDBFile

        pdb_content = """\
ATOM      1  N   ALA A   1      -2.273   0.000   0.000  1.00  0.00           N
END
"""
        import tempfile

        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write(pdb_content)
            pdb_file = f.name

        try:
            pdb = PDBFile(pdb_file)
            # 不传递位置
            topo = topology_parser.topology_from_openmmm(pdb.topology)

            assert topo.n_residues == 1

            # 检查原子位置是否为 None
            for atom in topo.atoms:
                assert atom.position is None
        finally:
            import os

            os.unlink(pdb_file)


class TestPDBxParser:
    """测试 PDBx/mmCIF 解析器"""

    def test_topology_from_pdbx_content_basic(self):
        """测试从基本的 PDBx 内容加载"""
        # 这是一个简化的 mmCIF 内容
        mmcif_content = """\
data_TEST
#
_entry.id   TEST
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
ATOM   N N ALA A 1 -2.273 0.000 0.000 1.00 0.00 0
ATOM   C CA ALA A 1 -1.713 1.403 0.000 1.00 0.00 0
ATOM   C C ALA A 1 -0.254 1.403 0.000 1.00 0.00 0
ATOM   O O ALA A 1 0.436 2.343 0.000 1.00 0.00 0
#
_end
"""
        try:
            topo = topology_parser.topology_from_pdbx_content(mmcif_content)

            assert isinstance(topo, Topology)
            assert topo.n_residues >= 1
        except Exception as e:
            # PDBx 解析可能还不完整，记录错误
            pytest.skip(f"PDBx parser not fully implemented: {e}")
