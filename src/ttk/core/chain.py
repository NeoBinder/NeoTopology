import numpy as np
import itertools
from ttk.core import Entity
from ttk.math.vector import calc_dihedral
import ttk


class Chain(Entity):
    """A Chain object represents a chain within a Topology.
    Attributes
    ----------
    index : int
        The index of the Chain within its Topology
    topology : mdtraj.Topology
        The Topology this Chain belongs to
    chain_id : str
        The PDB chainID of the Chain.
    residues : generator
        Iterator over all Residues in the Chain.
    atoms : generator
        Iterator over all Atoms in the Chain.
    """

    def __init__(self, index, topology, chain_id=None):
        """Construct a new Chain.  You should call add_chain() on the Topology instead of calling this directly."""
        # The index of the Chain within its Topology
        self.index = index
        # The Topology this Chain belongs to
        self.topology = topology
        self.residues = []
        # PDB format chainID
        self.chain_id = chain_id

    def __repr__(self):
        return "Chain {} with {} residue at {}".format(
            self.index, len(self.residues), hex(id(self))
        )

    #########################################
    # get methods
    #########################################
    def select_residues_by_atoms(self, atoms):
        residues = set()
        for atom in atoms:
            if atom.residue in self.residues:
                residues.add(atom.residue)
        return list(residues)

    def get_atom_indices(self):
        return [atom.index for atom in self.atoms]

    def residue(self, index):
        """Get a specific residue in this Chain.
        Parameters
        ----------
        index : int
            The index of the residue to select.
        Returns
        -------
        residue : Residue
        """
        return self.residues[index]

    @property
    def id(self):
        return self.chain_id

    @property
    def n_residues(self):
        """Get the number of residues in this Chain."""
        return len(self.residues)

    @property
    def atoms(self):
        """
        Returns all atoms
        -------
        """
        atoms = []
        for residue in self.residues:
            for atom in residue.atoms:
                atoms.append(atom)
        return atoms

    @property
    def n_atoms(self):
        """Get the number of atoms in this Chain"""
        return sum([r.n_atoms for r in self.residues])

    @property
    def is_valid(self):
        valid = True
        # check number of residues
        if len(self.residues) == 0:
            valid = False
        return valid

    #########################################
    # add methods
    #########################################
    def add_residue(self, residue):
        self.residues.append(residue)

    #########################################
    # delete methods
    #########################################
    def delete_residues(self, target_residues):
        self.residues = [res for res in self.residues if res not in target_residues]
        return True

    #########################################
    # compute methods
    #########################################
    def is_protein_chain(self):
        filtered_residues = [res for res in self.residues if not res.is_water]
        if len(filtered_residues) < 3:
            return False
        protein_res_num = sum([int(res.is_protein) for res in filtered_residues])
        if protein_res_num == 0:
            return False
        return protein_res_num / len(filtered_residues) > 0.9

    def calculate_psi_phi(self):
        if not self.is_protein_chain():
            return
        for index in range(self.n_residues):
            current_res = self.residues[index]
            if not current_res.is_valid_backbone:
                current_res.property_computed["psi"] = None
                current_res.property_computed["phi"] = None
                continue
            prev_res = self.residues[index - 1]

            n = current_res["N"].position.magnitude
            ca = current_res["CA"].position.magnitude
            c = current_res["C"].position.magnitude

            if current_res.res_seq - prev_res.res_seq == 1:
                prev_ca = prev_res["CA"].get_position(magnitude=True)
                phi = calc_dihedral(prev_ca, n, ca, c)
            else:
                phi = None

            if index == self.n_residues - 1:
                psi = None
            else:
                next_res = self.residues[index + 1]
                if next_res.res_seq - current_res.res_seq == 1:
                    next_natom = next_res.get("N")
                    if next_natom:
                        next_n = next_natom.position.magnitude
                        psi = calc_dihedral(n, ca, c, next_n)
                    else:
                        psi = None
                else:
                    psi = None

            current_res.property_computed["psi"] = psi
            current_res.property_computed["phi"] = phi
        pass

    def to_topology(self):
        atom_map = {}
        # create_topology
        top = ttk.Topology()
        # create the chain
        newchain = top.add_chain()

        for res in self.residues:
            newres = top.add_residue(
                res.name, newchain, res_seq=res.res_seq, segment_id=res.segment_id
            )
            for atom in res.atoms:
                newatom = top.add_atom(atom.name, atom.element, newres, atom.position)
                atom_map[atom] = newatom
        for atom in self.atoms:
            newatom = atom_map[atom]
            for bond in atom.bonds:
                top.add_bond(newatom, atom_map[bond.connect(atom)], bond.type)
        return top
