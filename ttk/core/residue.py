import numpy as np

from ttk.core import (
    _AMINO_ACID_CODES,
    _PROTEIN_RESIDUES,
    _SIMPLE_SOLVENT,
    _WATER_RESIDUES,
    Entity,
)


class Residue(Entity):
    """A Residue object represents a residue within a Topology.
    Attributes
    ----------
    name : str
        The name of the Residue
    index : int
        The index of the Residue within its Topology
    chain : mdtraj.topology.Chain
        The chain within which this residue belongs
    res_seq : int
        The residue sequence number
    segment_id : str, optional
        A label for the segment to which this residue belongs
    """

    def __init__(self, name, index, chain, res_seq, segment_id=""):
        """Construct a new Residue.  You should call add_residue()
        on the Topology instead of calling this directly."""
        self.name = name
        self.index = index
        self.chain = chain
        self.res_seq = res_seq
        self.segment_id = segment_id
        self.atoms = []
        self.property_computed = {}

    def __str__(self):
        return "%s-%s" % (self.name, self.res_seq)

    def __repr__(self):
        return "Residue:{} at {}".format(str(self), hex(id(self)))

    def __getitem__(self, key):
        for atom in self.atoms:
            if atom.clean_name == key:
                return atom
        raise KeyError

    def get(self, key, defalut=None):
        for atom in self.atoms:
            if atom.clean_name == key:
                return atom
        return defalut

    def add_atom(self, atom):
        self.atoms.append(atom)

    def atoms_by_name(self, name):
        """Iterator over all Atoms in the Residue with a specified name
        Parameters
        ----------
        name : str
            The particular atom name of interest.
        Examples
        --------
        >>> for atom in residue.atoms_by_name('CA'):
        ...     print(atom)
        Returns
        -------
        atomiter : generator
        """
        for atom in self.atoms:
            if atom.name.strip() == name:
                yield atom

    def get_atom_indices(self):
        return [atom.index for atom in self.atoms]

    @property
    def is_valid_backbone(self):
        atoms_set = set([atom.clean_name for atom in self.atoms])
        backbone_set = set(["N", "CA", "C"])
        return len(backbone_set.intersection(atoms_set)) == 3

    @property
    def n_atoms(self):
        """Get the number of atoms in this Residue"""
        return len(self.atoms)

    @property
    def is_protein(self):
        """Whether the residue is one found in proteins."""
        return self.name in _PROTEIN_RESIDUES

    @property
    def code(self):
        """Get the one letter code for this Residue"""
        if self.is_protein:
            return _AMINO_ACID_CODES[self.name]
        else:
            return None

    @property
    def is_water(self):
        """Whether the residue is water.
        Residue names according to VMD
        References
        ----------
        http://www.ks.uiuc.edu/Research/vmd/vmd-1.3/ug/node133.html
        """
        return self.name in _WATER_RESIDUES

    @property
    def is_solvent(self):
        return self.name in _SIMPLE_SOLVENT

    @property
    def mass(self):
        return sum([atom.element.mass for atom in self.atoms])
