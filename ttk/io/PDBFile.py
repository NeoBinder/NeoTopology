import datetime
from collections import Counter
import numpy as np

import ttk
from ttk import unit
from ttk.core import element_from_symbol, UnitCell, from_bond_float

from .PDBUtil import parse_pdb_header_list, split_pdb_content


class PDBParser:

    def __init__(self, config=None, **kwargs):
        if config is None:
            self.config = {}
        else:
            self.config = config
        self.config.update(kwargs)
        self.header = None
        self._current_topology = None
        self.topologies = []
        self._filetype = kwargs.get("filetype", "pdb")
        self._current_topology_id = 0

    @property
    def filetype(self):
        return "pdb"

    def load_from_content_list(self, content):
        if isinstance(content, str):
            content = content.split("\n")
        assert isinstance(content, list)
        header, models = split_pdb_content(content)
        self.header = parse_pdb_header_list(header)
        self.parse_models(models)
        if self.header["periodic_box"]:
            periodic_box = UnitCell.from_content(self.header["periodic_box"])
            for top in self.topologies:
                top.periodic_box = periodic_box
        if len(self.header["symmetry_rotation_matrix"].keys()):
            top.symmetry = self.header["symmetry_rotation_matrix"]

    def add_model(self, line=None):
        if line is None:
            self._current_topology = ttk.Topology()
            self._current_topology.id = 1
        else:
            self._current_topology = ttk.Topology()
            self._current_topology.id = int(line[10:14])
            self._current_topology_id = int(line[10:14])
        self._parse_dict = {}
        self.topologies.append(self._current_topology)

    def parse_chain(self, chainid):
        current_chain = self._parse_dict.get("current_chain")
        if current_chain is None or current_chain.id != chainid:
            current_chain = self._current_topology.add_chain(chainid)
            self._parse_dict["current_chain"] = current_chain
        return current_chain

    def parse_residue(self, resname, chain, resseq):
        current_residue = self._parse_dict.get("current_residue")
        if (
            current_residue is None
            or current_residue.name != resname
            or current_residue.res_seq != resseq
        ):
            current_residue = self._current_topology.add_residue(
                resname, chain, res_seq=resseq
            )
            self._parse_dict["current_residue"] = current_residue
        return current_residue

    def parse_atom(self, line, current_residue):
        fullname = line[12:16]
        # get rid of whitespace in atom names
        split_list = fullname.split()
        if len(split_list) != 1:
            # atom name has internal spaces, e.g. " N B ", so
            # we do not strip spaces
            name = fullname
        else:
            # atom name is like " CA ", so we can strip spaces
            name = split_list[0]
        altloc = line[16]
        try:
            serial_number = int(line[6:11])
        except Exception:
            serial_number = 0
        icode = line[26]  # insertion code
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        try:
            occupancy = float(line[54:60])
        except:
            occupancy = 1.0
        try:
            temperature_factor = float(line[60:66])
        except:
            temperature_factor = 0.0
        element_symbol = line[76:78].strip()
        if element_symbol == self.config.get("extraParticleIdentifier", "EP"):
            element = "EP"
        else:
            atom_element = element_from_symbol(element_symbol)
        try:
            formal_charge = int(line[78:80])
        except ValueError:
            formal_charge = None

        is_hetero = line[:6] == "HETATM"

        self._current_topology.add_atom(
            fullname,
            atom_element,
            current_residue,
            position=unit.Quantity(np.array([x, y, z]) / 10, "nm"),
            index=serial_number,
            formal_charge=formal_charge,
            temperature_factor=temperature_factor,
            occupancy=occupancy,
            is_hetero=is_hetero,
        )

    def parse_bonds(self, bonds):
        top = self.topologies[-1]
        atom_maps = {atom.index: atom for atom in top.atoms}
        bond_dict = {}
        for bond in bonds:
            bond = " ".join(
                [bond[6:11], bond[11:16], bond[16:21], bond[21:26], bond[26:31]]
            ).split()
            bond = list(map(int, bond))
            if len(bond) > 1:
                bond_dict[bond[0]] = Counter(bond[1:]) + bond_dict.get(
                    bond[0], Counter()
                )

        for source_atom_idx, atom_bond_dict in bond_dict.items():
            for target_atom_idx, bond_order in atom_bond_dict.items():
                top.add_bond(
                    atom_maps[source_atom_idx],
                    atom_maps[target_atom_idx],
                    from_bond_float(float(bond_order)),
                )
        pass

    def parse_models(self, models):
        bonds = []
        UNKNOWN_RECORD_TYPE = []
        for line in models:
            record_type = line[0:6].strip()
            if not line.strip():
                continue  # skip empty lines
            elif record_type in ("ATOM", "HETATM"):
                # Initialize the Model - there was no explicit MODEL record
                if self._current_topology is None:
                    self.add_model()
                chainid = line[21]
                current_chain = self.parse_chain(chainid)

                resname = line[17:20].strip()
                resseq = int(line[22:26].split()[0])  # sequence identifier
                current_residue = self.parse_residue(resname, current_chain, resseq)
                # atom name,element,residue,position
                current_atom = self.parse_atom(line, current_residue)

            elif record_type == "TER":
                self._parse_dict = {}
            elif record_type == "MODEL":
                self.add_model(line)
            elif record_type in ("END", "CONECT"):
                bonds.append(line)
            elif record_type == "ENDMDL":
                self._current_topology = None
                self._parse_dict = {}
            elif record_type in ("SIGUIJ", "SIGATM"):
                raise NotImplemented("valid value without implementation")
            elif record_type not in {
                "ATOM",
                "HETATM",
                "MODEL",
                "ENDMDL",
                "TER",
                "ANISOU",
                # These are older 2.3 format specs:
                "SIGATM",
                "SIGUIJ",
                # bookkeeping records after coordinates:
                "MASTER",
            }:
                UNKNOWN_RECORD_TYPE.append(record_type)
                #  raise ValueError("invalid value {}".format(line))
                print(
                    "There are unknown record type:{}\n Please pay attention or ignore it".format(
                        set(UNKNOWN_RECORD_TYPE)
                    )
                )
        self.parse_bonds(bonds)


class PDBFile:
    standardResidues = {
        "ALA",
        "ASN",
        "CYS",
        "GLU",
        "HIS",
        "HID",
        "LEU",
        "MET",
        "PRO",
        "THR",
        "TYR",
        "ARG",
        "ASP",
        "GLN",
        "GLY",
        "ILE",
        "LYS",
        "PHE",
        "SER",
        "TRP",
        "VAL",
        "A",
        "G",
        "C",
        "U",
        "I",
        "DA",
        "DG",
        "DC",
        "DT",
        "DI",
        "HOH",
    }

    def __init__(self, config=None, **kwargs):
        if config is None:
            self.config = {}
        else:
            self.config = config

        self.config.update(kwargs)

    @staticmethod
    def generate_header():
        return "REMARK CREATED WITH ttk {}, {} \n".format(
            ttk.__version__, datetime.datetime.now().strftime("%Y-%m-%d")
        )

    @staticmethod
    def generate_cryst(periodic_box):
        cryst1 = periodic_box.to_cryst1()
        content = "CRYST1"
        content += "{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}".format(*cryst1)
        content += " P 1           1\n"
        return content

    @staticmethod
    def topology2model(topology, config):
        content = ""
        index = 0
        for chain in topology.chains:
            for res in chain.residues:
                if res.name in PDBFile.standardResidues:
                    recordName = "ATOM  "
                else:
                    recordName = "HETATM"
                for atom in res.atoms:
                    desc = ""
                    index += 1
                    atom.index = index
                    if len(atom.symbol) > 3:
                        raise Exception("atom symbol greater than 3 letter")
                    desc += recordName
                    desc += "{:>5} ".format(index)
                    if len(atom.name) < 4:
                        desc += " {:<3} ".format(atom.name)
                    elif len(atom.name) == 4:
                        desc += "{:<4} ".format(atom.name)
                    else:
                        raise Exception(
                            "atom name '{}' is longger than 4 !!".format(atom.name)
                        )
                    if len(res.name) > 3:
                        resname = res.name[:3]
                    else:
                        resname = res.name
                    desc += "{:>3} ".format(resname)
                    desc += chain.id
                    if config.get("use_res_index"):
                        desc += "{:>4} ".format(res.index)
                    else:
                        desc += "{:>4} ".format(res.res_seq)
                    desc += "   "
                    desc += "{:8.3f}{:8.3f}{:8.3f}".format(
                        *(atom.position.to("angstrom").magnitude)
                    )
                    desc += "{:6.2f}{:6.2f}".format(atom.occupancy, atom.bfactor)
                    desc += "      "
                    desc += res.segment_id.rjust(4, " ")
                    desc += atom.symbol.rjust(2, " ")
                    content += desc + "\n"
            index += 1
            desc = "TER   {:>5}      {} {} {:>3}\n".format(
                index, resname, chain.id, res.res_seq
            )
            content += desc
        return content

    @staticmethod
    def topology2connect(topology):
        bond_atoms = []
        content = ""
        for atom in topology.get_atoms():
            if atom.residue.name in PDBFile.standardResidues:
                continue
            bond_atoms.append(atom)
        for atom in bond_atoms:
            bonds = [bond for bond in atom.bonds]
            idx = 0
            desc = ""
            while len(bonds):
                if idx % 4 == 0:
                    if idx != 0:
                        desc += "\n"
                    desc += "CONECT" + str(atom.index).rjust(5, " ")
                bond = bonds.pop()
                connect_atom = bond.connect(atom)
                desc += str(connect_atom.index).rjust(5, " ")
                idx += 1
            content += desc
            content += "\n"
        return content

    def to_content(self, topology):
        content = ""
        content += self.generate_header()
        if topology.periodic_box is not None:
            content += self.generate_cryst(topology.periodic_box)
        content += self.topology2model(topology, self.config)
        content += self.topology2connect(topology)
        return content

    def write_file(self, topology, filehandler):
        content = self.to_content(topology)
        filehandler.write(content)
