# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 12:45:27 2022

@author: eccn3
"""

import math

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole

import networkx as nx
import networkx.algorithms.isomorphism as iso

from pysmiles import write_smiles

import itertools

class graph_enumerator:
    
    def __init__(self, graphs = [], rd_mol = [], substrates = None):
        
        for i in graphs:
            self.systems = [{'graph': i,
                             'smiles': None,
                             'rdkit_mol': None,
                             'subgraphs': [],
                             'substrate': None
                             }]
        if substrates:
            for i,j in zip(substrates, self.systems):
                j['substrate'] = i
                
        self.graphs = {'graphs': graphs}
                      
        self.descriptors = []
        self.subgraph_matrix = None
        self.smiles = None
         
    def regen_check(self, function, check, regenerate = False):
        
        """
        Adds option to regenerate all entries, or just ones that fail a 
        boolean check for a dictionary entry
        
        in: 
            self, function to perform, dict key to check, and regenerate boolean
            
        out:
            performs function on all entries that do not pass the check
        """
        
        
        if regenerate:
            for i in self.systems:
                function(i)
                
        if not regenerate:
            for i in self.systems:
                if not i[check]:
                    function(i)
    
    def gen_smiles(self, entry):
        
        """
        Generates a SMILES string. Rounds up all fractional bond orders.
        
        in:
            self, dict (requires 'graph')
        
        out: 
            self, dict with SMILES string in dict['smiles']
        
        """
        
        graph = entry['graph']
        
        for edge in graph.edges():
            
            if 'weight' in graph.edges[edge]:
                
                graph.edges[edge]['order'] = graph.edges[edge]['weight']
                
                if isinstance(graph.edges[edge]['order'], float):
                    graph.edges[edge]['order'] = math.ceil(graph.edges[edge]['order'])
                    
        entry['smiles'] = write_smiles(graph)
    


    def gen_rdkitmol(self, entry):

        """
        Generates an rdkit mol object based on SMILES string. 
        
        in:
            self, dict (requires 'smiles')
        
        out: 
            self, dict with rdkit mol in dict['rdkit_mol']
        
        """        

        if not entry['smiles']: 
            entry['smiles'] = self.gen_smiles(entry['graph'])
        smiles = entry['smiles'].replace('[M]', f"[{entry['substrate']}]")
        entry['rdkit_mol'] = Chem.MolFromSmiles(smiles, sanitize=False)
                
        
    def gen_subgraphs(self, smiles = False, regenerate = False):
        
        if smiles:

            for i in self.systems:
                if i['smiles'] in smiles:
                    i['subgraphs'] = self.get_subgraphs(i['graph'])
            
        
    def get_subgraphs(graph, max_graph_length=5):
        
        if not max_graph_length:
            max_graph_length = graph.number_of_nodes()
        
        all_connected_subgraphs = []
        for nb_nodes in range(1, max_graph_length+1):
            for subgraph in (graph.subgraph(selected_nodes) for selected_nodes in itertools.combinations(graph, nb_nodes)):
                if nx.is_connected(subgraph):
                    all_connected_subgraphs.append(subgraph)
                
        return all_connected_subgraphs

    # def get_subgraphs(graph, max_graph_length=5):
        
        

    def check_iso(graph1, graph2, soft_check=True):
        """ Check isomorphism using node attribute {atomic_number} and
        edge attribute {type}.
        Parameters
        ----------
        graph1: networkx Graph object.
        graph2: networkx Graph object.
        Returns
        -------
            bool
            True is graph1 and graph 2 are isomorphic.
        """
        
        if soft_check:
            quick_check = nx.faster_could_be_isomorphic(graph1, graph2)
            if not quick_check:
                return quick_check
        
        return nx.is_isomorphic(graph1, graph2,
                                node_match=iso.categorical_node_match(
                                    'element', 'C'),
                                edge_match=iso.categorical_edge_match('weight', 1))

    def update_X_descriptors(self, X, descriptors, graphs):
        """Updates coefficient matrix with counts of graph descriptors
        and creates new column for descriptor not seen before.
        Parameters
        ---------
        X: list
            Coefficient matrix with each row for  a molecule and
            each column for a coefficient of descriptor in molecule
        descriptor: list
            Unique graphs which form basis set of molecule expansion
        graphs: list
            Subgraphs of a molecule.
        Returns
        -------
            list, list
        (updated) X, descriptors
        """
        row_entry = list()
        row_entry.extend([0] * len(descriptors))
        for graph in (graphs):
            graph_classified = False # if graph has been assigned a descriptor
            for index, descriptor in enumerate(descriptors):
                if not len(descriptor.nodes()) == len(graph.nodes()) or \
                        not len(descriptor.edges()) == len(graph.edges()):
                    continue
                elif self.check_iso(graph, descriptor):
                    row_entry[index] +=1
                    graph_classified = True
                    break  # stop going thru descriptors
            if not graph_classified:
                # graph is a descriptor not seen previously
                descriptors.append(graph)
                row_entry.append(1)
        X.append(row_entry)
        
        return X, descriptors

    def get_graph_descriptors(self, graphs, max_length=5):
        
        X,descriptors = list(), list()
        
        if not max_length:
            max_length = max([len(i.nodes()) for i in graphs])
        
        for count, graph in enumerate(graphs):
            
            print('Enumerating subgraphs in graph ', count)
            
            subgraphs = self.get_subgraphs(graph, max_graph_length = max_length)
            
            X, descriptors = self.update_X_descriptors(X=X, descriptors=descriptors, graphs=subgraphs)
            
            assert len(X[-1]) == max([len(row) for row in X])
            
            def complete_rows(n):
                """ Given a 2D array, make it a square,
                i.e. make all rows of length n by padding
                the missing entries with zero.
                Parameters
                ----------
                n: int
                    Size of max row.
                Returns
                -------
                    list
                    Squared matrix
                """
                X_out = []
                for row in X:
                    row.extend([0] * (n - len(row)))
                    # row_entry.extend([0] * len(descriptors))
                    X_out.append(row)
                return X_out
            
            X = complete_rows(n=len(X[-1]))
    
        return X, descriptors
    
    
if __name__=='__main__':
    import pickle as pk
    import pandas as pd
    from rdkit.Chem import Draw
    IPythonConsole.ipython_useSVG=True
    
    if not 'df' in locals() or 'df' in globals():
        df = pd.DataFrame(pk.load(open('Examples/graphs/LDAtoPBE_df.pk', 'rb')))
        graph = df.loc[86]['Active Graph']
    
    obj = graph_enumerator([graph], substrates = ['Pt'])
    
    obj.regen_check(obj.gen_smiles, 'smiles', regenerate=True)
    obj.regen_check(obj.gen_rdkitmol, 'rdkit_mol', regenerate=True)
    obj.regen_check(obj.gen_subgraphs, 'subgraphs', regenerate = True)
    print(obj.systems[0])
