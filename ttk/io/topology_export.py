import ttk
from ttk.io import PDBFile


def topology_from_indices(topology, indices, start_idx=1):
    if start_idx != 0:
        indices = [idx - start_idx for idx in indices]
    top = ttk.Topology()
    top.periodic_box = topology.periodic_box
    indices_set = set(indices)
    atom_map = {}
    for chain in topology.chains:
        c = top.add_chain(chain.id)
        for residue in chain.residues:
            residue_atom_indices = set([atom.index for atom in residue.atoms])
            intersection = residue_atom_indices.intersection(indices_set)
            if len(intersection) == 0:
                continue
            r = top.add_residue(residue.name, c, residue.res_seq, residue.segment_id)
            for atom in residue.atoms:
                if atom.index in intersection:
                    newatom = top.add_atom(
                        atom.name,
                        atom.element,
                        r,
                        atom.position,
                        origin_index=atom.index,
                    )
                    atom_map[atom] = newatom
    for bond in topology.get_bonds():
        a1, a2 = bond.atom1, bond.atom2
        if a1.index in indices_set and a2.index in indices_set:
            top.add_bond(atom_map[a1], atom_map[a2], bond.type)
    return top


def topology_to_pdb(ttk_topology, fname=None):
    content = PDBFile().to_content(ttk_topology)
    if fname:
        with open(fname, "w") as f:
            f.write(content)
    return content


def topology_to_openmm(ttk_topology):
    from openmm.app import Element, Topology
    from openmm.unit import nanometer

    newTopology = Topology()
    newTopology.setPeriodicBoxVectors(ttk_topology.periodic_box.to_openmm())
    newAtoms = {}
    newPositions = [] * nanometer
    for chain in ttk_topology.chains:
        newChain = newTopology.addChain(chain.id)
        for residue in chain.residues:
            newResidue = newTopology.addResidue(
                residue.name, newChain, residue.res_seq, residue.segment_id
            )
            for atom in residue.atoms:
                newAtom = newTopology.addAtom(
                    atom.name,
                    Element.getBySymbol(atom.element.symbol),
                    newResidue,
                    atom.index,
                )
                newAtoms[hash(atom)] = newAtom
                # newPositions.append(atom.position * nanometer)
                newPositions.append(
                    atom.position.to(ttk.unit.nanometer).magnitude * nanometer
                )
    for bond in ttk_topology.get_bonds():
        newTopology.addBond(newAtoms[hash(bond.atom1)], newAtoms[hash(bond.atom2)])
    return newTopology, newPositions


def topology2geom(model, charge):
    topology = model
    is_open_shelled = False
    geom = ""
    line = " {:3} {: > 7.3f} {: > 7.3f} {: > 7.3f} \n"
    total_elec = 0
    for atom in topology.get_atoms():
        #  x, y, z = positions[i][0], positions[i][1], positions[i][2]
        element = atom.element
        symbol = element.symbol
        total_elec += element.atomic_number
        # pos = (atom.position * 10).tolist()
        pos = atom.position.to("angstrom").magnitude
        geom += line.format(symbol, *pos)

    if total_elec % 2 != 0:
        total_elec += charge
        if total_elec % 2 != 0:
            is_open_shelled = True

    return geom, is_open_shelled
