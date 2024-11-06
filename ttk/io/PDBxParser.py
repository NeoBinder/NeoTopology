import re

import numpy as np
import pandas as pd

import ttk
from ttk import unit
from ttk.core import UnitCell, element_from_symbol, from_bond_name

reserved_state_categories = {"loop", "global", "save", "stop"}
# match the line starting with _ and split the line into key: value.
# value can also be quoted within '' or ""
# ?: means non-capturing group:
# \S+: words
# \s+: white spaces
content_re = re.compile(r"(?:_(.+?)[.](\S+))\s+['\"]?([^'\"]*)['\"]?")
category_re = re.compile(r"(?:_(.+?)[.](\S+))")
table_re = re.compile(r"(?:'([^']*)'|(\S+))")


class DefinitionContainer:

    def __init__(self):
        self.data_type = "definition"


class DataCategoriesContainer:

    def __init__(self, group):
        self.data_type = "data_categories"
        self.group = group
        self.categories = []
        self.data = []

    def add_category(self, category):
        self.categories.append(category)

    def add_value(self, value):
        self.data.append(value)

    def to_dataframe(self, index_name="id"):
        df = pd.DataFrame(self.data, columns=self.categories)
        if index_name:
            df[index] = df[index_name]
            df = df.set_index(index_name)
        return df

    @property
    def name(self):
        return self.group


class DataContainer:

    def __init__(self, category):
        self.data_type = "data"
        self.category = category
        self.data = {}

    def setvalue(self, key, value):
        self.data[key] = value

    @property
    def name(self):
        return self.category


def split_chunks(content):
    chunks = []
    chunk = []
    for line in content:
        if line.startswith("#"):
            if len(chunk):
                chunks.append(chunk)
            chunk = []
        else:
            chunk.append(line)
    return chunks


def process_group(chunk, current_group):
    if current_group == "LOOP":
        group, category = category_re.findall(chunk[1])[0]
        current_container = DataCategoriesContainer(group)
        for line in chunk[1:]:
            if line.startswith("_"):
                # categories
                group, category = category_re.findall(line)[0]
                current_container.add_category(category)
            else:
                array = ["".join(match) for match in table_re.findall(line)]
                current_container.add_value(array)
    else:
        raise NotImplementedError(
            "current_category {} not implemented".format(current_group)
        )
    return current_container


def process_category(chunk):
    line = chunk[0]
    category, key, value = content_re.findall(line)[0]
    current_container = DataContainer(category[1:])

    new_chunk = []
    for line in chunk:
        if not line.startswith("_"):
            new_chunk[-1] = new_chunk[-1] + " " + line
        else:
            new_chunk.append(line)
    for line in new_chunk:
        # parse DataContainer
        _, key, value = content_re.findall(line)[0]
        assert _ == category
        current_container.setvalue(key, value)
    return current_container


class PDBxParser:

    def __init__(self, config=None, **kwargs):
        if config is None:
            self.config = {"remove_water": False}
        else:
            config = config
        self.config.update(**kwargs)
        self.unknown_data = []
        self.data_dict = {}
        self.data_categories_dict = {}

        self.topologies = []
        self.current_topology = None
        self.filetype = kwargs.get("filetype", "pdbx")

    def parse_from_content(self, content):
        if isinstance(content, str):
            content = content.split("\n")
        assert isinstance(content, list)

        chunks = split_chunks(content)

        for chunk in chunks:
            line = chunk[0]
            idx = line.find("_")
            if idx == -1:
                current_category = "UNKNOWN"
            current_category = line[:idx].upper()
            if current_category.lower() in reserved_state_categories:
                container = process_group(chunk, current_category)
                self.data_categories_dict[container.group] = container
            elif line.startswith("_"):
                container = process_category(chunk)
                self.data_dict[container.category] = container
            else:
                self.unknown_data.append(chunk)
        self.parse_data()

    def parse_data(self):
        # parse chains
        self.pares_model()
        # parse bond
        self.parse_bond()

        # parse unitcell last step
        cell = self.data_dict.get("cell", None)
        if cell:
            cell_parameter = cell.data
            periodic_box = UnitCell.from_parameter(
                a_length=cell_parameter["length_a"],
                b_length=cell_parameter["length_b"],
                c_length=cell_parameter["length_c"],
                alpha=cell_parameter["angle_alpha"],
                beta=cell_parameter["angle_beta"],
                gamma=cell_parameter["angle_gamma"],
            )

    def parse_model(self):
        # assume 1 top per file
        top = ttk.Topology()
        chains_map = {}
        residues_map = {}
        atoms_map = {}

        df = self.data_categories_dict["atom_site"].to_dataframe()
        for idx, row in df.iterrows():
            chainid = row["label_asym_id"]
            if chainid not in chains_map:
                chains_map[chainid] = top.add_chain(chainid)

            resname = row["label_comp_id"]
            res_seq = row["label_seq_id"]
            res = top.add_residue(resname, chains_map[chainid], res_seq=res_seq)
            residues_map["{}-{}".format(chainid, res_seq)] = res

            is_hetero = row["group_PDB"] == "ATOM"
            atom = top.add_atom(
                row["label_atom_id"],
                element_from_symbol(row["type_symbol"]),
                res,
                position=unit.Quantity(
                    np.array(row["Cartn_x"], row["Cartn_y"], row["Cartn_x"]).astype(
                        float
                    )
                    / 10,
                    "nm",
                ),
                index=int(idx),
                formal_charge=row["pdbx_formal_charge"],
                occupancy=row["occupancy"],
                is_hetero=is_hetero,
            )
            atoms_map["{}-{}".format(chainid, atom.index)] = atom

        for idx, row in self.data_categories_dict["struct_conn"].to_dataframe(
            index_name=None
        ):
            src_atom = atoms_map[
                "{}-{}".format(row["ptnr1_label_asym_id"], row["ptnr1_label_atom_id"])
            ]
            trg_atom = atoms_map[
                "{}-{}".format(row["ptnr2_label_asym_id"], row["ptnr2_label_atom_id"])
            ]
            top.add_bond(src_atom, trg_atom, from_bond_name(row["conn_type_id"]))

    def parse_bond_template(self):
        templates = {}
        for idx, row in self.data_categories_dict["chem_comp_bond"].to_dataframe(
            index_name="comp_id"
        ):
            pass
