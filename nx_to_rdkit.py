# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 02:55:45 2023

@author: eccn3
"""

import networkx as nx
from networkx import weisfeiler_lehman_graph_hash as wlgh

from AtomsHelper.utils import draw_surf_graph, remove_small_edges

from rdkit import Chem
from rdkit.Chem.rdmolops import FindAllSubgraphsOfLengthMToN
from rdkit.Chem.rdchem import BondType

from ase import formula
from ase.data import chemical_symbols
from ase.visualize import view

from collections import Counter

class rdkit_helper():
    def __init__(self, nx_graph = None, rdkit_mol = None, substrate = 'Pt', edge_property = 'weight'):
        self.substrate = substrate
        self.gen_m_g = None
        self.g = None
        self.edge_property = edge_property
        
        metals = [i for i in chemical_symbols if i not in formula.non_metals]
        
        if nx_graph:
            gen_m_graph = 'M' in nx.get_node_attributes(nx_graph, 'element').values()
            if gen_m_graph: 
                self.gen_m_g = nx_graph
                self.gen_metal_assign()
            else:
                self.g = nx_graph
            
        self.mol = rdkit_mol
        
        self.rdkit_bond_chart = {0.: BondType.ZERO,
                                 0.5: BondType.SINGLE,
                                 1.: BondType.SINGLE,
                                 1.5: BondType.ONEANDAHALF,
                                 2.: BondType.DOUBLE,
                                 2.5: BondType.TWOANDAHALF,
                                 3.: BondType.TRIPLE, 
                                 3.5: BondType.THREEANDAHALF,
                                 4.: BondType.QUADRUPLE
                                 }

    def gen_metal_assign(self):
        
        """
        If there exists a generalizes metal graph, assigns substrate element
        to generalized metal nodes
        
        In
        -------
        self: self.gen_m_g

        Returns
        -------
        self: self.g

        """
        
        if self.gen_m_g:
            G = self.gen_m_g.copy()
            
            for node in G.nodes():
                if G.nodes[node]['element'] == 'M':
                    G.nodes[node]['element'] = self.substrate
            
            self.g = G
        else:
            print('No generalized metal graph')

    def nx_to_mol(self, gen_m = False):
        
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
        
        if not self.g:
            print('No applicable nx graph')
            
        else:
            
            G = self.g.copy()
            
            mol = Chem.RWMol()
            atomic_nums = nx.get_node_attributes(G, 'element')
            node_to_idx = {}
            edge_to_idx = {}
            for node in G.nodes():
                
                a=Chem.Atom(atomic_nums[node])
                idx = mol.AddAtom(a)
                mol.GetAtomWithIdx(idx).SetIntProp('nx_idx', int(node))
                node_to_idx[node] = idx
        
            bond_weight = nx.get_edge_attributes(G, self.edge_property)
            bond_length = nx.get_edge_attributes(G, 'distance')
            
            for edge in G.edges():
                
                first, second = edge
                ifirst = node_to_idx[first]
                isecond = node_to_idx[second]
                bond = bond_weight[first, second]
                
                if G.nodes[first]['element'] == G.nodes[second]['element'] == self.substrate:
                    bond_type = self.rdkit_bond_chart[1.]
                    bond_index = 1.
                
                else:
                    bond_index = float(round(bond*2))/2
                    bond_type = self.rdkit_bond_chart[bond_index]
                if bond_length:
                    bond_distance = round(bond_length[first, second], 5)
                
                mol.AddBond(ifirst, isecond, bond_type)
                edge_to_idx[first, second] = (ifirst, isecond)
                mol.GetBondBetweenAtoms(ifirst, isecond).SetDoubleProp('bond_index', bond_index)
                mol.GetBondBetweenAtoms(ifirst, isecond).SetDoubleProp('weight', bond)
                
                if bond_length:
                    mol.GetBondBetweenAtoms(ifirst, isecond).SetDoubleProp('distance', bond_distance)
            
            self.mol = mol
            self.g_node_map = node_to_idx
            self.g_edge_map = edge_to_idx
            
            return mol

    def mol_to_nx(self, mol, bondids=None):
        """
        Straight up, Himaghna's function.
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

    def enumerate_subgraphs(self, min_path=1, max_path=2, useHs=True,
                            edge_attr='bond_index', node_attr='element'):
        
        """
        Enumerate subgraphs 
        
        In
        -------
        self: self.mol

        Returns
        -------
        self: self.subgraphs: tuple of a list of tuples

        """
        
        self.mol_subgraphs = FindAllSubgraphsOfLengthMToN(self.mol,
                                                      min = min_path,
                                                      max = max_path,
                                                      useHs = useHs)
        
        nx_subgraphs = []
        nx_subgraph_hashes = []
        # self.mol_subgraphs is split into into a list of lists;
        # Subgraphs are split into lists where each list is all subgraphs with 
        # len(edges) = n

        
        for i in self.mol_subgraphs:
            for j in i:
                subgraph = self.mol_to_nx(self.mol, j)
                nx_subgraphs.append(subgraph)
                nx_subgraph_hashes.append(wlgh(subgraph, iterations=max_path,
                                               edge_attr = edge_attr,
                                               node_attr = node_attr))
        
        elements = Counter(list(nx.get_node_attributes(self.g, 'element').values()))
        
        def atom_graph(element):
            graph =nx.Graph()
            graph.add_node(0)
            graph.nodes[0]['element'] = element
            return graph
        
        element_graphs = {i: atom_graph(i) for i in elements.keys()}
        element_hashes = {i: wlgh(element_graphs[i], 
                                  iterations=max_path,
                                  edge_attr = edge_attr,
                                  node_attr = node_attr) for i in elements.keys()}
        
        for i in elements.keys():
            for j in range(elements[i]):
                nx_subgraphs.insert(0, element_graphs[i])
                nx_subgraph_hashes.insert(0, element_hashes[i])
        
        self.nx_subgraphs = nx_subgraphs
        self.nx_subgraph_hashes = nx_subgraph_hashes
        
        self.nx_subgraph_counter = Counter(self.nx_subgraph_hashes)
        
        hash_index= {i: [] for i in self.nx_subgraph_counter.keys()}
        for count, i in enumerate(self.nx_subgraph_hashes):
            hash_index[i].append(count)
        
        self.nx_subgraph_indices = hash_index
        
        dic = {'nx_subgraphs': nx_subgraphs,
               'nx_subgraph_hashes': nx_subgraph_hashes,
               'nx_subgraph_counter': self.nx_subgraph_counter,
               'nx_subgraph_indices': self.nx_subgraph_indices}
        
        return dic
        
    def remove_small_bonds(self, measure = 'rdkit_bond'):
        
        """
        Enumerate subgraphs 
        
        In
        -------
        self: self.mol 

        Returns
        -------
        self: self.mol (cutoff small bonds)

        """
        
        bondids = [i.GetIdx() for i in list(self.mol.GetBonds())]
        removes = []
        
        if measure == 'rdkit_bond':
            
            for bondid in bondids:
                
                Bond = self.mol.GetBondWithIdx(bondid)
    
                bond_type= Bond.GetBondType()
                
                if bond_type.name == 'ZERO':
                
                    begin_atom_idx = Bond.GetBeginAtomIdx()
                    end_atom_idx = Bond.GetEndAtomIdx()
                    removes.append((begin_atom_idx, end_atom_idx))
                
        for i in removes:
            self.mol.RemoveBond(i[0], i[1])
            
    def print_mol_atoms(self):
        
        """
        Returns self.mol atom properties
        
        In
        -------
        self: self.mol

        Returns
        -------
        print: atom idx, atomic number

        """
        
        if not self.mol:
            print('No rdkit mol')
            return
        else:
            for atom in self.mol.GetAtoms():
                print(atom.GetIdx(),
                      atom.GetAtomicNum(),
                      )
                
    def print_mol_bonds(self):
        
        """
        Returns self.mol bond properties
        
        In
        -------
        self: self.mol

        Returns
        -------
        print: atom idx (begin, end), dictionary of properties
        
        
        """
        if not self.mol:
            print('No rdkit mol')
            return
        else:
            for bond in self.mol.GetBonds():
                print(bond.GetBeginAtomIdx(),
                      bond.GetEndAtomIdx(),
                      bond.GetPropsAsDict())
                
    # def mol_with_atom_index(self):
        
    #         """
    #         Returns self.mol bond properties
            
    #         In
    #         -------
    #         self: self.mol
    
    #         Returns
    #         -------
    #         print: atom idx (begin, end), dictionary of properties
            
    #         """
    #     for atom in self.mol.GetAtoms():
    #         atom.SetAtomMapNum(atom.GetIdx())
    #     return mol

def simplify_graph(graph, substrate = 'Pt', cutoff = 0.01):
        
    graph = remove_small_edges(graph, cutoff = cutoff)
    
    keep = []
    
    for node in graph.nodes:
        if graph.nodes[node]['element'] == substrate:
            continue
        else:
            neighbors = list(graph.neighbors(node))
            keep.append(node)
            for i in neighbors:
                if graph.edges[i, node]['weight'] > cutoff:
                    keep.append(i)
            
    keeps = set(keep)
    simple_graph = graph.subgraph(keeps).copy()
    
    for i in simple_graph.edges:
        if simple_graph.nodes[i[0]]['element'] not in ['C','H','O'] and \
            simple_graph.nodes[i[1]]['element'] not in ['C','H','O']:
            simple_graph.remove_edge(i[0], i[1])
            
    
    return  simple_graph,keeps
    
if __name__ == '__main__':
    import pickle as pk
    
    from rdkit.Chem.Draw import IPythonConsole
    from rdkit.Chem import Draw, AllChem
    from ase.io import read
    
    from AtomsHelper.atoms_helper import geom2graph
    from AtomsHelper.ddec6_graph_gen import gen_graph_ddec6
    from AtomsHelper.utils import remove_small_edges, draw_surf_graph
    import numpy as np
    from collections import Counter

    IPythonConsole.ipython_useSVG=True
    
    # graph = pk.load(open('Examples/nx_mol_helper/OCH3-Pt-hol.pk', 'rb'))
    # ddec6 = gen_graph_ddec6('Examples/CHOH/ddec6')
    # ddec6, _ = simplify_graph(ddec6)
    # ddec6 = remove_small_edges(ddec6)
    # atoms = [read('Examples/ase_to_graph/CONTCAR_Pt_CHOH_ontop'), read('Examples/ase_to_graph/CONTCAR_CHCH_holhol')]
    # geoms = [geom2graph(atoms = atom) for atom in atoms]
    examples = ['Examples/ddec6/CHOH',
                'Examples/ddec6/CH2CHCH3CH3']
    
    cutoff = 0.3
    
    ddec6 = [gen_graph_ddec6(i) for i in examples]
    for count, i in enumerate(ddec6):
        ddec6[count], _  = simplify_graph(i, cutoff = cutoff)
        draw_surf_graph(ddec6[count], cutoff = cutoff)
    rdhelp = [rdkit_helper(nx_graph=i) for i in ddec6]
    
    for i in rdhelp:
        i.nx_to_mol()
        i.remove_small_bonds()
        
        i.smiles = Chem.MolToSmiles(i.mol)
        g = i.mol_to_nx(i.mol)
        
        i.enumerate_subgraphs(min_path=1, max_path = 4)
        i.count = Counter(i.nx_subgraph_hashes)
        
    
    
    # rd_help = rdkit_helper(nx_graph = ddec6)
    # rd_help.nx_to_mol()
    # rd_help.remove_small_bonds()
    # # rd_help.mol.UpdatePropertyCache()
    # rd_help.smiles = Chem.MolToSmiles(rd_help.mol)
    # rd_help.enumerate_subgraphs(min_path= 1, max_path = 4)
    # count = Counter(rd_help.nx_subgraph_hashes)
    # hash_index= {i: [] for i in count.keys()}
    # for count, i in enumerate(rd_help.nx_subgraph_hashes):
        # hash_index[i].append(count)
    
    # new_rd_mol = Chem.MolFromSmiles(rd_help.smiles)
    # san_rd_mol = Chem.SanitizeMol(new_rd_mol)

    # benz = Chem.MolFromSmiles('C1CCCCC1')
    # benz_morgan = np.array(AllChem.GetMorganFingerprintAsBitVect(benz,2))
    # benz_locs = np.where(benz_morgan !=0 )[0]
    # octz = Chem.MolFromSmiles('C1CCCCCCC1')
    # octz_morgan = np.array(AllChem.GetMorganFingerprintAsBitVect(octz,2))
    # octz_locs = np.where(octz_morgan !=0 )[0]
    