"""Tests for ttk.io.topology_export and related residue/chain properties."""


from ttk import Topology, unit
from ttk.core import element_from_symbol
from ttk.io.topology_export import topology_to_pdb


class TestTopologyToPdb:
    """Tests for topology_to_pdb function."""

    def test_topology_to_pdb_returns_string(self):
        """topology_to_pdb should return a non-empty PDB string."""
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain, res_seq=1)
        topo.add_atom(
            "N",
            element_from_symbol("N"),
            res,
            is_hetero=False,
            position=unit.Quantity([-2.273, 0.000, 0.000], unit.nanometer),
        )
        topo.add_atom(
            "CA",
            element_from_symbol("C"),
            res,
            is_hetero=False,
            position=unit.Quantity([-1.713, 1.403, 0.000], unit.nanometer),
        )
        topo.create_index()

        content = topology_to_pdb(topo)

        assert isinstance(content, str)
        assert len(content) > 0

    def test_topology_to_pdb_contains_atom_records(self):
        """PDB output should contain ATOM records."""
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain, res_seq=1)
        topo.add_atom(
            "CA",
            element_from_symbol("C"),
            res,
            is_hetero=False,
            position=unit.Quantity([1.0, 2.0, 3.0], unit.nanometer),
        )
        topo.create_index()

        content = topology_to_pdb(topo)

        assert "ATOM" in content or "HETATM" in content

    def test_topology_to_pdb_writes_file(self, tmp_path):
        """topology_to_pdb should write a file when fname is given."""
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain, res_seq=1)
        topo.add_atom(
            "CA",
            element_from_symbol("C"),
            res,
            is_hetero=False,
            position=unit.Quantity([1.0, 2.0, 3.0], unit.nanometer),
        )
        topo.create_index()

        output_file = tmp_path / "output.pdb"
        topology_to_pdb(topo, str(output_file))

        assert output_file.exists()
        assert output_file.stat().st_size > 0

    def test_topology_to_pdb_roundtrip(self):
        """A topology exported to PDB and re-imported should preserve basic structure."""
        from ttk.io.topology_parser import topology_from_pdb_content

        topo = Topology()
        chain = topo.add_chain()
        res1 = topo.add_residue("ALA", chain, res_seq=1)
        topo.add_atom(
            "N",
            element_from_symbol("N"),
            res1,
            is_hetero=False,
            position=unit.Quantity([-0.227, 0.000, 0.000], unit.nanometer),
        )
        topo.add_atom(
            "CA",
            element_from_symbol("C"),
            res1,
            is_hetero=False,
            position=unit.Quantity([-0.171, 0.140, 0.000], unit.nanometer),
        )
        res2 = topo.add_residue("GLY", chain, res_seq=2)
        topo.add_atom(
            "N",
            element_from_symbol("N"),
            res2,
            is_hetero=False,
            position=unit.Quantity([0.100, 0.200, 0.000], unit.nanometer),
        )
        topo.create_index()

        pdb_content = topology_to_pdb(topo)
        reloaded = topology_from_pdb_content(pdb_content)

        assert reloaded.n_residues == 2
        assert reloaded.n_chains == 1

    def test_topology_to_pdb_heteroatom(self):
        """HETATM atoms should appear as HETATM in PDB output."""
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("LIG", chain, res_seq=1)
        topo.add_atom(
            "C1",
            element_from_symbol("C"),
            res,
            is_hetero=True,
            position=unit.Quantity([1.0, 2.0, 3.0], unit.nanometer),
        )
        topo.create_index()

        content = topology_to_pdb(topo)

        assert "HETATM" in content


class TestResidueProperties:
    """Tests for Residue property methods."""

    def _make_residue(self, name):
        topo = Topology()
        chain = topo.add_chain()
        return topo.add_residue(name, chain)

    def test_is_protein_for_ala(self):
        res = self._make_residue("ALA")
        assert res.is_protein is True

    def test_is_protein_for_gly(self):
        res = self._make_residue("GLY")
        assert res.is_protein is True

    def test_is_protein_false_for_water(self):
        res = self._make_residue("HOH")
        assert res.is_protein is False

    def test_is_protein_false_for_ligand(self):
        res = self._make_residue("LIG")
        assert res.is_protein is False

    def test_is_water_for_hoh(self):
        res = self._make_residue("HOH")
        assert res.is_water is True

    def test_is_water_for_wat(self):
        res = self._make_residue("WAT")
        assert res.is_water is True

    def test_is_water_false_for_ala(self):
        res = self._make_residue("ALA")
        assert res.is_water is False

    def test_is_solvent_for_na(self):
        res = self._make_residue("NA")
        assert res.is_solvent is True

    def test_is_solvent_false_for_ala(self):
        res = self._make_residue("ALA")
        assert res.is_solvent is False

    def test_code_for_ala(self):
        res = self._make_residue("ALA")
        assert res.code == "A"

    def test_code_for_gly(self):
        res = self._make_residue("GLY")
        assert res.code == "G"

    def test_code_none_for_water(self):
        res = self._make_residue("HOH")
        assert res.code is None

    def test_mass_increases_with_atoms(self):
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)
        topo.add_atom("N", element_from_symbol("N"), res, is_hetero=False)
        topo.add_atom("CA", element_from_symbol("C"), res, is_hetero=False)

        assert res.mass > 0

    def test_n_atoms(self):
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)
        topo.add_atom("N", element_from_symbol("N"), res, is_hetero=False)
        topo.add_atom("CA", element_from_symbol("C"), res, is_hetero=False)

        assert res.n_atoms == 2

    def test_is_valid_backbone_true(self):
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)
        topo.add_atom("N", element_from_symbol("N"), res, is_hetero=False)
        topo.add_atom("CA", element_from_symbol("C"), res, is_hetero=False)
        topo.add_atom("C", element_from_symbol("C"), res, is_hetero=False)

        assert res.is_valid_backbone is True

    def test_is_valid_backbone_false(self):
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)
        topo.add_atom("CB", element_from_symbol("C"), res, is_hetero=False)

        assert res.is_valid_backbone is False


class TestChainProperties:
    """Tests for Chain property methods."""

    def _make_protein_chain(self, n_residues=5):
        topo = Topology()
        chain = topo.add_chain()
        for i in range(n_residues):
            topo.add_residue("ALA", chain, res_seq=i + 1)
        return chain

    def test_n_residues(self):
        chain = self._make_protein_chain(3)
        assert chain.n_residues == 3

    def test_n_atoms(self):
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)
        topo.add_atom("N", element_from_symbol("N"), res, is_hetero=False)
        topo.add_atom("CA", element_from_symbol("C"), res, is_hetero=False)

        assert chain.n_atoms == 2

    def test_atoms_property(self):
        topo = Topology()
        chain = topo.add_chain()
        res = topo.add_residue("ALA", chain)
        atom1 = topo.add_atom("N", element_from_symbol("N"), res, is_hetero=False)
        atom2 = topo.add_atom("CA", element_from_symbol("C"), res, is_hetero=False)

        assert atom1 in chain.atoms
        assert atom2 in chain.atoms

    def test_is_protein_chain_true(self):
        chain = self._make_protein_chain(10)
        assert chain.is_protein_chain() is True

    def test_is_protein_chain_false_too_short(self):
        chain = self._make_protein_chain(2)
        assert chain.is_protein_chain() is False

    def test_is_protein_chain_false_all_water(self):
        topo = Topology()
        chain = topo.add_chain()
        for i in range(10):
            topo.add_residue("HOH", chain, res_seq=i + 1)
        assert chain.is_protein_chain() is False

    def test_is_valid_empty_chain(self):
        topo = Topology()
        chain = topo.add_chain()
        assert chain.is_valid is False

    def test_is_valid_chain_with_residues(self):
        chain = self._make_protein_chain(3)
        assert chain.is_valid is True

    def test_residue_accessor(self):
        chain = self._make_protein_chain(3)
        assert chain.residue(0) == chain.residues[0]
        assert chain.residue(1) == chain.residues[1]

    def test_delete_residues(self):
        topo = Topology()
        chain = topo.add_chain()
        res1 = topo.add_residue("ALA", chain, res_seq=1)
        res2 = topo.add_residue("GLY", chain, res_seq=2)
        res3 = topo.add_residue("VAL", chain, res_seq=3)

        chain.delete_residues([res2])

        assert res1 in chain.residues
        assert res2 not in chain.residues
        assert res3 in chain.residues

    def test_select_residues_by_atoms(self):
        topo = Topology()
        chain = topo.add_chain()
        res1 = topo.add_residue("ALA", chain, res_seq=1)
        atom1 = topo.add_atom("CA", element_from_symbol("C"), res1, is_hetero=False)
        res2 = topo.add_residue("GLY", chain, res_seq=2)
        topo.add_atom("CA", element_from_symbol("C"), res2, is_hetero=False)

        result = chain.select_residues_by_atoms([atom1])

        assert res1 in result
        assert res2 not in result
