import ttk


class Atom(object):
    """An Atom object represents a residue within a Topology.
    Attributes
    ----------
    name : str
        The name of the Atom
    element : mendeleev.element.Element
        The element of the Atoms
    index : int
        The index of the Atom within its Topology
    residue : Residue
        The Residue this Atom belongs to
    serial : int
        The serial number from the PDB specification. Unlike index,
        this may not be contiguous or 0-indexed.
    """

    def __init__(self, name, element, residue, is_hetero, **kwargs):
        """Construct a new Atom.  You should call add_atom() on the Topology instead of calling this directly."""
        # The name of the Atom
        self.name = name
        # That Atom's element
        self.element = element
        # The Residue this Atom belongs to
        self.residue = residue
        self.bonds_dict = {}

        # position in nm
        self._position = None
        self.is_hetero = is_hetero
        self.properties = {**kwargs, "clean_name": name.replace(" ", "")}

    @property
    def bonds(self):
        return list(self.bonds_dict.values())

    @property
    def clean_name(self):
        return self.properties["clean_name"]

    def __str__(self):
        return "{}-{}".format(self.residue, self.name)

    def __repr__(self):
        return "<Atom:{} at {}>".format(str(self), hex(id(self)))

    def __eq__(self, value):
        return id(self) == id(value)

    def __hash__(self):
        return id(self)

    @property
    def index(self):
        return self.properties.get("index", 0)

    @index.setter
    def index(self, value):
        self.properties["index"] = value

    @property
    def symbol(self):
        return self.element.symbol

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, value):
        if value is None:
            return
        assert isinstance(value, ttk.unit.Quantity)
        self._position = value

    def get_position(self, unit=ttk.unit.nanometer, magnitude=False):
        res = self.position.to(unit)
        if magnitude:
            res = res.magnitude
        return res

    def get_bond_dict(self):
        return {bond.connection_hash: bond for bond in self.bonds}

    @property
    def occupancy(self):
        return self.properties.get("occupancy", 1.00)

    @property
    def bfactor(self):
        return self.properties.get("occupancy", 0.00)

    #  @property
    #  def n_bonds(self):
    #  """Number of bonds in which the atom participates."""
    #  # TODO: this info could be cached.
    #  return ilen(bond for bond in self.residue.chain.topology.bonds
    #  if self in bond)

    @property
    def is_backbone(self):
        """Whether the atom is in the backbone of a protein residue"""
        return self.name in set(["C", "CA", "N", "O"]) and self.residue.is_protein

    @property
    def is_sidechain(self):
        """Whether the atom is in the sidechain of a protein residue"""
        return (
            self.name not in set(["C", "CA", "N", "O", "HA", "H"])
            and self.residue.is_protein
        )

    @property
    def segment_id(self):
        """User specified segment_id of the residue to which this atom belongs"""
        return self.residue.segment_id
