# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 12:45:27 2022

@author: eccn3
"""

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
import numpy as np


import networkx as nx

from AtomsHelper.atoms_helper import atoms_helper
from AtomsHelper.nx_rdkit_helper import nx_rdkit_helper

class graph_enumerator:
    
    def __init__(self, input_data: list, surface_type = 'Pt', cutoff = 0.3):
        
        if all(isinstance(i, type(input_data[0])) for i in input_data):
            self.input_data = input_data
            
            if type(self.input_data[0]) == nx.graph.Graph:
                self.input_data_type = nx.graph.Graph
                self.helpers = [atoms_helper(i, surface_type = surface_type, bond_cutoff = cutoff, convert = False) for i in self.input_data]
                [i.gen_ads_subgraph() for i in self.helpers]
                self.systems = [nx_rdkit_helper(i.ads_subgraph) for i in self.helpers]
                    
            elif all([isinstance(i, Chem.rdchem.Mol) for i in input_data]):
                self.input_data_type = Chem.rdchem.Mol
                self.systems = [nx_rdkit_helper(rdkit_mol = i) for i in input_data]
            
            elif all([isinstance(i, atoms_helper) for i in input_data]):
                self.input_data_type = nx.graph.Graph
                self.helpers = input_data
                [i.gen_ads_subgraph() for i in self.helpers]
                self.systems = [nx_rdkit_helper(i.ads_subgraph) for i in self.helpers]
            
        for i in self.systems:
            if self.input_data_type == nx.graph.Graph:
                i.graph_to_mol()
            elif self.input_data_type == Chem.rdchem.Mol:
                i.mol_to_graph()
            
    def enumerate_subgraphs(self, min_path = 1, max_path = 3):
        for i in self.systems:
            i.enumerate_subgraphs(min_path = min_path, max_path = max_path)

    
    def collect_system_subgraphs(self):
        
        hashes = set()

        
        for i in self.systems:
            hashes.update(i.nx_subgraph_hashes)

        hashes = list(hashes)
        
        subgraphs = dict()
        
        for i in hashes:
            for j in self.systems:
                if i in j.nx_subgraph_hashes:
                    subgraphs[i] = j.nx_subgraphs[j.nx_subgraph_indices[i][0]]
        

        # self.subgraph_matrix = np.zeros((len(hashes), len(self.systems)))
        self.subgraph_counter = {i: np.zeros(len(self.systems)) for i in hashes}
        
        for i in hashes:
            for count2, j in enumerate(self.systems):
            
                if i in j.nx_subgraph_hashes:
                    self.subgraph_counter[i][count2] = j.nx_subgraph_counter[i]
                    
        self.subgraph_matrix = np.vstack(np.array([self.subgraph_counter[i] for i in hashes]))
        

        return hashes, subgraphs

if __name__=='__main__':
    import pickle as pk
    import pandas as pd
    from rdkit.Chem import Draw
    from AtomsHelper.plotting import draw_surf_graph as dsg

    
    IPythonConsole.ipython_useSVG=True
    
    if not 'df' in locals() or 'df' in globals():
        df = pd.DataFrame(pk.load(open('Examples/graphs/LDAtoPBE_df.pk', 'rb')))
        graph1 = df.loc[86]['Active Graph']
        graph2 = df.loc[87]['Active Graph']
        graph3 = df.loc[88]['Active Graph']
    
    obj = graph_enumerator([graph1, graph2, graph3], surface_type = 'Pt')
    obj.enumerate_subgraphs()
    hashes, subgraphs = obj.collect_system_subgraphs()
