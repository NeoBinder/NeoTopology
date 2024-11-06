import hashlib


class BondType(object):

    def __init__(self, name, order):
        self.name = name
        self.order = order

    def __float__(self):
        return self.order

    def __repr__(self):
        return "BondType:{} at {}".format(self.name, hex(id(self)))

    def __eq__(self, bondtype):
        return self.name == bondtype.name or self.order == bondtype.order

    def __hash__(self):
        hashv = hashlib.sha256(self.name.encode()).hexdigest()
        return int(hashv, 16)

    def to_rdkit(self):
        from rdkit.Chem import BondType as rdkit_BondType

        name2rdkit = {
            "Single": rdkit_BondType.SINGLE,
            "Double": rdkit_BondType.DOUBLE,
            "Aromatic": rdkit_BondType.AROMATIC,
        }
        return name2rdkit[self.name]

    @classmethod
    def Single(cls):
        return cls("Single", 1.0)

    @classmethod
    def Double(cls):
        return cls("Double", 2.0)

    @classmethod
    def Triple(cls):
        return cls("Triple", 3.0)

    @classmethod
    def Aromatic(cls):
        return cls("Aromatic", 1.5)

    @classmethod
    def Amide(cls):
        return cls("Amide", 1.25)

    @classmethod
    def Unknown(cls):
        return cls("Unknown", 0)


Single = BondType.Single()
Double = BondType.Double()
Triple = BondType.Triple()
Aromatic = BondType.Aromatic()
Amide = BondType.Amide()
BondTypeUnknown = BondType.Unknown()

float2bondtype = {
    None: BondTypeUnknown,
    1.0: Single,
    1.25: Amide,
    1.5: Aromatic,
    2.0: Double,
    3.0: Triple,
}
bondname2bondtype = {"covale": Single}


def from_rdkit_bond(bond):
    from rdkit.Chem import BondType as rdkit_BondType

    bond_dict = {
        rdkit_BondType.SINGLE: Single,
        rdkit_BondType.DOUBLE: Double,
        rdkit_BondType.AROMATIC: Aromatic,
    }
    return bond_dict[bond.GetBondType()]


def from_bond_float(bond_float):
    """
    Convert a float to known bond type class, or None if no matched class if found
    Parameters
    ----------
    bond_float : float
        Representation of built in bond types as a float,
        Maps this float to specific bond class, if the float has no map, None is returned instead
    Returns
    -------
    bond_type : mdtraj.topology.Singleton subclass or None
        Bond type matched to the float if known
        If no match is found, returns None (which is also a valid type for the Bond class)
    """
    return float2bondtype[bond_float]


def from_bond_name(bond_name):
    return bondname2bondtype[bond_name]


class Bond:
    """A Bond representation of a bond between two Atoms within a Topology
    Attributes
    ----------
    atom1 : mdtraj.topology.Atom
        The first atom in the bond
    atom2 : mdtraj.topology.Atom
        The second atom in the bond
    order : instance of mdtraj.topology.Singleton or None
    type : int on [1,3] domain or None
    """

    def __init__(self, atom1, atom2, bondtype=None):
        self.atom1 = atom1
        self.atom2 = atom2
        if bondtype is None or isinstance(bondtype, float):
            bondtype = from_bond_float(bondtype)
        self.type = bondtype

    def __hash__(self):
        return hash(self.atom1) ^ hash(self.atom2) ^ hash(self.type)

    @property
    def connection_hash(self):
        return hash(self.atom1) ^ hash(self.atom2)

    def __eq__(self, bond):
        if bond.type == self.type:
            if self.atom1 in [bond.atom1, bond.atom2] and self.atom2 in [
                bond.atom1,
                bond.atom2,
            ]:
                return True
        return False

    def __repr__(self):
        output = "Bond({}, {}".format(self.atom1, self.atom2)
        if self.type is not None:
            output += ", type={}".format(self.type)
        output += ")"
        output += "at {}".format(hex(id(self)))
        return output

    def connect(self, atom):
        if atom is self.atom1:
            return self.atom2
        elif atom is self.atom2:
            return self.atom1
        else:
            raise Exception("Bond and Atom have no correspondence")

    #  def __getnewargs__(self):
    #  """
    #  Support for pickle protocol 2:
    #  http://docs.python.org/2/library/pickle.html#pickling-and-unpickling-normal-class-instances
    #  """
    #  return self[0], self[1], self.type, self.order

    #  @property
    #  def _equality_tuple(self):
    #  # Hierarchy of parameters: Atom1 index -> Atom2 index -> type -> order
    #  return (self[0].index, self[1].index,
    #  float(self.type) if self.type is not None else 0.0,
    #  self.order if self.order is not None else 0)

    #  def __deepcopy__(self, memo):
    #  return Bond(self[0], self[1], self.type, self.order)

    #  def __eq__(self, other):
    #  if not isinstance(other, Bond):
    #  return False
    #  return self._equality_tuple == other._equality_tuple

    #  def __hash__(self):
    #  # Set of atoms making up bonds, the type, and the order
    #  return hash((self[0], self[1], self.type, self.order))

    #  def __gt__(self, other):
    #  # Cannot use total_ordering because namedtuple
    #  # has its own __gt__, __lt__, etc. methods, which
    #  # supersede total_ordering
    #  self._other_is_bond(other)
    #  return self._equality_tuple > other._equality_tuple

    #  def __ge__(self, other):
    #  self._other_is_bond(other)
    #  return self._equality_tuple >= other._equality_tuple

    #  def __lt__(self, other):
    #  self._other_is_bond(other)
    #  return self._equality_tuple < other._equality_tuple

    #  def __le__(self, other):
    #  self._other_is_bond(other)
    #  return self._equality_tuple <= other._equality_tuple

    #  def __ne__(self, other):
    #  return not self.__eq__(other)

    #  def __getstate__(self):
    #  # This is required for pickle because the parent class
    #  # does not properly return a state
    #  return self.__dict__
