# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 23:25:23 2024

@author: eccn3
"""

import numpy as np

from ase import Atoms, Atom
from ase.data import covalent_radii

import itertools
import re

import networkx as nx


from AtomsHelper.utils import min_img_dist

from rdkit import Chem
from rdkit.Chem.rdchem import BondType


def graph_to_atoms(graph, periodic_boundary_condition = [1,1,1],
                   graph_keys = {'element': 'symbol',
                                 'position': 'position',
                                 'dipole': 'magmom',
                                 'charge': 'charge'}):
    atoms = Atoms()
    atoms.pbc = periodic_boundary_condition
    if 'cell' in graph.graph.keys():
        atoms.cell = graph.graph['cell']
    
    
    for i in graph.nodes:
        atom_properties = dict()
        for j in graph_keys:
            if j in graph.nodes[i]:
                atom_properties[graph_keys[j]] = graph.nodes[i][j]
                
        atoms.append(Atom(**atom_properties))
    
    return atoms

def atoms_to_graph(atoms,
                  distance_tune = 1.2):
    
    graph = nx.Graph()
    graph.graph['cell'] = atoms.cell
    graph.graph['pbc'] = atoms.pbc
    
    
    elements = set([i.number for i in atoms])
    radii = {i: covalent_radii[i] for i in elements}
    bond_distance = {frozenset(i): distance_tune*(radii[i[0]]+ radii[i[1]]) for i in itertools.combinations_with_replacement(elements, 2) }
    
    
    for count, i in enumerate(atoms):
        graph.add_node(count,
                       element = i.symbol,
                       number = i.number,
                       position = i.position,
                       dipole = i.magmom,
                       charge = 'charge',
                       )

    for i, j in itertools.combinations(list(graph.nodes), 2):
        
        
        dist = np.linalg.norm(graph.nodes[i]['position']-graph.nodes[j]['position'])

        if all(graph.graph['pbc']) and dist > 4:
            dist = min_img_dist(atoms.cell, 
                                graph.nodes[i]['position'],
                                graph.nodes[j]['position'])
            
        else:
            dist = np.linalg.norm(graph.nodes[i]['position']-graph.nodes[j]['position'])
        
        
        if dist <= bond_distance[frozenset((graph.nodes[i]['number'], graph.nodes[j]['number']))]:
            graph.add_edge(i,j,distance = dist)
             
    return graph


def ddec6_to_graph(directory = None, 
                   bader = False,
                   surface = None,
                   ase_atoms = False,
                   weight_cutoff = 0.2):

    if directory:
        import os
        thisdir = os.getcwd()
        
        
        if any(['DDEC6' in i for i in directory.rsplit('/', 1)]):
            directory = directory.rsplit('/',1)[0]
        
        os.chdir(directory)
        
    periodic = np.array([p for p in itertools.product([-1, 0, 1], repeat=3)])

    nac = open('DDEC6_even_tempered_net_atomic_charges.xyz')
    
    graph = nx.Graph()
    
    num_atoms = int(next(nac))
    next(nac)
    
    if surface == None:
        surface = ['Pt', 'Ni', 'Pd', 'Ir']
    else:
        surface = [surface]
    
    for count, i in enumerate(itertools.islice(nac, num_atoms)):
        
        element, x, y, z, charge = i.split()
        if count == 0:
            if element in surface:
                graph.graph['surface']=element
            else:
                graph.graph['surface']='gas'
        if any(element == i for i in surface):
            typ = 'slab'
        else:
            typ = 'adsorbate'
        
        graph.add_node(count,
                       element = element,
                       system_type = typ,
                       position = [float(x),float(y),float(z)],
                       charge = float(charge))
        
    record = False
    for i in nac:
        if record:
            idx, _, _, _,_,_, x, y, z, mag, xy, xz, yz, x2y2, z2, ev1, ev2, ev3 = i.split()
            idx = int(idx)-1
            graph.nodes[idx]['dipole'] = [float(x),float(y),float(z)]
            graph.nodes[idx]['dipole_mag'] = mag
            graph.nodes[idx]['quadrupole'] = [float(j) for j in [xy, xz, yz, x2y2, z2]]
            graph.nodes[idx]['tless_q_ev'] = [float(ev1),float(ev2),float(ev3)]
            if idx + 1 == num_atoms:
                break
        if 'atom number, atomic symbol, x, y, z, net_charge, dipole_x,' in i:
            record = True

    bo = open('DDEC6_even_tempered_bond_orders.xyz')
    next(bo)
    
    cell_text = next(bo)
    cell_text_split = re.split('[{}]', cell_text)
    
    cell_x = np.array([float(i) for i in cell_text_split[3].split()])
    cell_y = np.array([float(i) for i in cell_text_split[5].split()])
    cell_z = np.array([float(i) for i in cell_text_split[7].split()])                  
    
    cell_dims = np.vstack([cell_x, cell_y, cell_z])
    
    graph.graph['cell'] = cell_dims
    
    for count, i in enumerate(itertools.islice(bo, num_atoms)):
        element, x, y, z, total = i.split()
        
        position = [float(i) for i in [x,y,z]]
        
        if graph.nodes[count]['position'] == position:
            graph.nodes[count]['bond_order_total'] = float(total)
        else:
            print(count, "Uh oh it doesn't match up")
    
    for i in bo:
        if 'Printing BOs for ATOM #' in i:
            index = int(i.split()[5])-1
            next(bo)
            j=next(bo)
            while 'Bonded to the' in j:
                j_split = j.split()
                
                weight = float(j_split[20])
     
                index2 = int(j_split[12])-1
                

                periodic_img_pos = graph.nodes[index2]['position'] + np.matmul(periodic, cell_dims)
                
                distance = np.min(np.linalg.norm(np.subtract(graph.nodes[index]['position'],
                                       periodic_img_pos), axis = 1))
                
                
                graph.add_edge(index, index2, 
                               weight = weight,
                               distance = distance)

                j = next(bo)
                
    if directory:
        os.chdir(thisdir)
        
        
    return graph


def graph_to_mol(graph, 
              metal = 'Pt', 
              gen_m = False, 
              edge_property = 'weight',
              return_mapping = True):
    
    """
    Converts networkx graph (self.g) into an rdkit mol
    
    Bonds are rounded to nearest half, but decimal values are still 
    kept via bond property
    
    In
    -------
    self: self.g

    Returns
    -------
    self: self.mol

    """
    
    rdkit_bond_chart = {0.: BondType.ZERO,
                        0.5: BondType.SINGLE,
                        1.: BondType.SINGLE,
                        1.5: BondType.ONEANDAHALF,
                        2.: BondType.DOUBLE,
                        2.5: BondType.TWOANDAHALF,
                        3.: BondType.TRIPLE, 
                        3.5: BondType.THREEANDAHALF,
                        4.: BondType.QUADRUPLE
                        }
    
    g = graph.copy()
    
    mol = Chem.RWMol()
    atomic_nums = nx.get_node_attributes(g, 'element')
    node_to_idx = {}
    edge_to_idx = {}
    for node in g.nodes():
        
        a=Chem.Atom(atomic_nums[node])
        idx = mol.AddAtom(a)
        mol.GetAtomWithIdx(idx).SetIntProp('nx_idx', int(node))
        node_to_idx[node] = idx

    bond_weight = nx.get_edge_attributes(g, edge_property)
    bond_length = nx.get_edge_attributes(g, 'distance')
    
    for edge in g.edges():
        
        first, second = edge
        ifirst = node_to_idx[first]
        isecond = node_to_idx[second]
        bond = bond_weight[first, second]
        
        if g.nodes[first]['element'] == g.nodes[second]['element'] == metal:
            bond_type = rdkit_bond_chart[1.]
            bond_index = 1.
        
        else:
            bond_index = float(round(bond*2))/2
            bond_type = rdkit_bond_chart[bond_index]
        if bond_length:
            bond_distance = round(bond_length[first, second], 5)
        
        mol.AddBond(ifirst, isecond, bond_type)
        edge_to_idx[first, second] = (ifirst, isecond)
        mol.GetBondBetweenAtoms(ifirst, isecond).SetDoubleProp('bond_index', bond_index)
        mol.GetBondBetweenAtoms(ifirst, isecond).SetDoubleProp('weight', bond)
        
        if bond_length:
            mol.GetBondBetweenAtoms(ifirst, isecond).SetDoubleProp('distance', bond_distance)
    
    mol = mol
    node_to_idx
    edge_to_idx
    
    if return_mapping:
        return mol, node_to_idx, edge_to_idx
    else:
        return mol

def mol_to_graph(mol, bondids=None):
    """
    Straight up, Himaghna Bhattacharjee's function.
    Check it out here:https://github.com/himaghna
    
    rdkit mol to nx graph
    
    Parameters
    ----------
    mol: RDKIT molecule object
    bondids: tuple
        Bond-ids making the graph.
    Returns
    -------
        networkx Graph object.
    """
    g = nx.Graph()
    
    mol_idx = set()
    
    if not bondids:
        bondids = [i.GetIdx() for i in list(mol.GetBonds())]
    
    for bondid in bondids:
        
        Bond = mol.GetBondWithIdx(bondid)
        
        bond_props = Bond.GetPropsAsDict()
        
        begin_atom_idx = Bond.GetBeginAtomIdx()
        end_atom_idx = Bond.GetEndAtomIdx()
        
        if begin_atom_idx not in mol_idx:
        
            g.add_node(begin_atom_idx,
                        atomic_number=(mol.GetAtomWithIdx(begin_atom_idx))
                                .GetAtomicNum(),
                        element=(mol.GetAtomWithIdx(begin_atom_idx))
                                .GetSymbol())
            
            mol_idx.add(begin_atom_idx)
            
        if end_atom_idx not in mol_idx:
            
            g.add_node(end_atom_idx,
                        atomic_number=(mol.GetAtomWithIdx(end_atom_idx))
                                .GetAtomicNum(),
                        element=(mol.GetAtomWithIdx(end_atom_idx))
                                .GetSymbol())
            mol_idx.add(begin_atom_idx)
            
            
        g.add_edge(begin_atom_idx, end_atom_idx,
                   **bond_props)

    return g
    
if __name__=='__main__':
    from AtomsHelper.atoms_helper import atoms_helper
    from AtomsHelper.plotting import draw_surf_graph as dsg
    from AtomsHelper.graph_ops import remove_small_edges  
    from rdkit.Chem.Draw import IPythonConsole
    
    graph = ddec6_to_graph('Examples/ddec6/CH2CHCH3CH3/DDEC6_even_tempered_bond_orders.xyz')
    
    dsg(graph)
    
    helper = atoms_helper(remove_small_edges(graph))
    
    helper.gen_ads_subgraph(1)
    
    mol = graph_to_mol(helper.ads_subgraph)
    
    mol[0]
    
    mol_to_graph_test = mol_to_graph(mol[0])
    
    dsg(mol_to_graph_test, edge_labels = True)
    
   #  graph = ddec6_read('Examples/ddec6/CHOH')
   #  atoms = graph_to_atoms(graph)
   #  min_img_dist(atoms.cell, atoms[-1].position, atoms[-3].position)
   #  graph2 = atoms_to_graph(atoms)
    
   # # view(atoms)
   #  draw_surf_graph(graph2)
