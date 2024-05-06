# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 02:55:45 2023

@author: eccn3
"""

import networkx as nx
from networkx import weisfeiler_lehman_graph_hash as wlgh


from AtomsHelper.conversions import graph_to_mol, mol_to_graph

from rdkit import Chem
from rdkit.Chem.rdmolops import FindAllSubgraphsOfLengthMToN
from rdkit.Chem.rdchem import BondType

from collections import Counter

class nx_rdkit_helper():
    def __init__(self, nx_graph = None, rdkit_mol = None, substrate = 'Pt', edge_property = 'weight'):
        self.substrate = substrate
        self.gen_m_g = None
        self.g = None
        self.edge_property = edge_property
        
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

    def graph_to_mol(self, **kwargs):
        
        self.mol, self.node_to_idx, self.edge_to_idx = graph_to_mol(self.g, **kwargs)
        

    def mol_to_graph(self, **kwargs):
        
        self.g = mol_to_graph(self.mol, **kwargs)

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
                subgraph = mol_to_graph(self.mol, j)
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
        
    def remove_small_bonds_in_mol(self, measure = 'rdkit_bond', max_bond = 0.5):
        
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


    
if __name__ == '__main__':
    
    from rdkit.Chem.Draw import IPythonConsole

    from AtomsHelper.conversions import ddec6_to_graph
    from AtomsHelper.plotting import draw_surf_graph
    from AtomsHelper.graph_ops import simplify_graph


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
    
    ddec6 = [ddec6_to_graph(i) for i in examples]
    for count, i in enumerate(ddec6):
        ddec6[count]  = simplify_graph(i, cutoff = cutoff)
        draw_surf_graph(ddec6[count], cutoff = cutoff)
    rdhelp = [nx_rdkit_helper(nx_graph=i) for i in ddec6]
    
    for i in rdhelp:
        i.graph_to_mol()
        i.remove_small_bonds_in_mol()
        
        i.smiles = Chem.MolToSmiles(i.mol)
        g = i.mol_to_graph()
        
        i.enumerate_subgraphs(min_path=1, max_path = 4)
        i.count = Counter(i.nx_subgraph_hashes)
        
  