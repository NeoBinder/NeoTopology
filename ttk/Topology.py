import copy

from ttk.core import Atom, Bond, Chain, Entity, Residue
from ttk.data import _PROTEIN_LETTERS


def expand_symmetry(top):
    newtop = Topology()
    symmetry_matrices = top.symmetry
    for k, rmatrix in symmetry_matrices.items():
        current_top = copy.deepcopy(top)
        positions = current_top.positions
        current_top.positions = rmatrix.apply(positions.magnitude) * positions.units
        newtop.add_topology(current_top)
    return newtop


class Topology(Entity):

    def __init__(self):
        self.chains = []
        self.periodic_box = None

    #########################################
    # get operation
    #########################################
    @property
    def n_chains(self):
        """Get the number of chains in the Topology"""
        return len(self.chains)

    @property
    def n_residues(self):
        """Get the number of residues in the Topology."""
        n_res = sum([len(chain.residues) for chain in self.chains])
        return n_res

    @property
    def residues(self):
        residues = []
        for chain in self.chains:
            for res in chain.residues:
                residues.append(res)
        return residues

    @property
    def atoms(self):
        atoms = []
        for chain in self.chains:
            for res in chain.residues:
                for atom in res.atoms:
                    atoms.append(atom)
        return atoms

    # get operation
    def get_atom_by_indices(self, indices):
        atoms = []
        for atom in self.atoms:
            if atom.index in indices:
                atoms.append(atom)
        return atoms

    def get_chains(self):
        return self.chains

    def get_residues(self):
        for chain in self.chains:
            for res in chain.residues:
                yield res

    def get_atoms(self, create_index=False):
        if create_index:
            self.create_index()
        for chain in self.chains:
            for res in chain.residues:
                for atom in res.atoms:
                    yield atom

    def get_bonds(self):
        bonds = set()
        for chain in self.chains:
            for res in chain.residues:
                for atom in res.atoms:
                    for bond in atom.bonds:
                        if bond not in bonds:
                            bonds.add(bond)
                            yield bond

    def get_residues_by_name(self, name_set):
        if isinstance(name_set, str):
            name_set = set([name_set])
        if isinstance(name_set, list):
            name_set = set(name_set)
        for residue in self.get_residues():
            if residue.name in name_set:
                yield (residue)

    #########################################
    # add operation
    #########################################
    def add_atom(
        self, name, element, residue, is_hetero, position=None, index=0, **kwargs
    ):
        """Create a new Atom and add it to the Topology.
        Parameters
        ----------
        name : str
            The name of the atom to add
        element : ttk.core.Element
            The element of the atom to add
        residue : ttk.core.Residue
            The Residue to add it to
        Returns
        -------
        atom : ttk.core.Atom
            the newly created Atom
        """
        atom = Atom(name, element, residue, is_hetero, **kwargs)
        residue.add_atom(atom)
        if position is not None:
            atom.position = position
        atom.index = index
        return atom

    def add_bond(self, atom1, atom2, bondtype=None):
        """Create a new bond and add it to the Topology.
        Parameters
        ----------
        atom1 : ttk.core.Atom
            The first Atom connected by the bond
        atom2 : ttk.core.Atom
            The second Atom connected by the bond
        bondtype : ttk.core.Bond.BondType or None, Default: None, Optional
            Bond type of the bond, or None if not known/provided
        """
        bond = Bond(atom1, atom2, bondtype=bondtype)
        if bond.connection_hash not in atom1.bonds_dict:
            atom1.bonds_dict[bond.connection_hash] = bond
        if bond.connection_hash not in atom1.bonds_dict:
            atom2.bonds_dict[bond.connection_hash] = bond
        return bond

    def add_chain(self, chain_id=None):
        if chain_id is None:
            exist_ids = set([c.id for c in self.chains])
            for p_id in _PROTEIN_LETTERS:
                if p_id not in exist_ids:
                    chain_id = p_id
                    break
        chain = Chain(len(self.chains), self, chain_id)
        self.chains.append(chain)
        return chain

    def add_residue(self, name, chain, res_seq=None, segment_id=""):
        """Create a new Residue and add it to the Topology.
        Parameters
        ----------
        name : str
            The name of the residue to add
        chain : ttk.core.Chain
            The Chain to add it to
        res_seq : int, optional
            Residue sequence number, such as from a PDB record. These sequence
            numbers are arbitrary, and do not necessarily start at 0 (or 1).
            If not supplied, the res Seq attribute will be set to the
            residue's sequential (0 based) index.
        segment_id : str, optional
            A label for the segment to which this residue belongs
        Returns
        -------
        residue : ttk.core.Residue
            The newly created Residue
        """
        if res_seq is None:
            res_seq = str(self.n_residues)
        res = Residue(name, len(chain.residues), chain, res_seq, segment_id)
        chain.add_residue(res)
        return res

    def add_topology(self, add_topology):
        #  if not isinstance(add_topology, Topology):
        #  raise Exception("add_topology only accept ttk Topology")
        atom_map = {}
        for chain in add_topology.chains:
            newchain = self.add_chain()
            for res in chain.residues:
                newres = self.add_residue(
                    res.name, newchain, res_seq=res.res_seq, segment_id=res.segment_id
                )
                for atom in res.atoms:
                    newatom = self.add_atom(
                        atom.name,
                        atom.element,
                        newres,
                        atom.is_hetero,
                        position=atom.position,
                    )
                    atom_map[atom] = newatom
        for atom in add_topology.get_atoms():
            newatom = atom_map[atom]
            for bond in atom.bonds:
                self.add_bond(newatom, atom_map[bond.connect(atom)], bond.type)

        self.create_index()

        #  self.chains += new_topology.chains

    #########################################
    # delete operation
    #########################################
    def remove_water(self):
        for chain in self.chains:
            delete_res = [res for res in chain.residues if res.is_water]
            chain.delete_residues(delete_res)
        self.chains = [chain for chain in self.chains if chain.is_valid]
        return True

    def remove_solvent(self):
        for chain in self.chains:
            delete_res = [res for res in chain.residues if res.is_solvent]
            chain.delete_residues(delete_res)
        self.chains = [chain for chain in self.chains if chain.is_valid]
        return True

    #########################################
    # compute operation
    #########################################

    #########################################
    # modify operation
    #########################################
    def create_index(self, atom_start_idx=0, res_start_idx=0):
        atom_idx = atom_start_idx
        res_idx = res_start_idx
        for chain in self.chains:
            for res in chain.residues:
                res.index = res_idx
                res_idx += 1
                for atom in res.atoms:
                    atom.index = atom_idx
                    atom_idx += 1
