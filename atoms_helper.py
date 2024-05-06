# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 13:23:21 2022

@author: eccn3
"""

import networkx as nx

from ase.data import chemical_symbols
from ase import formula
from ase.io import read

from AtomsHelper.conversions import graph_to_atoms, atoms_to_graph, ddec6_to_graph
from AtomsHelper import graph_ops

class atoms_helper:
    
    """
    
    Object containing vasp structure - to - networkx conversion
    
    :vasp: file location of structure to convert
    
    """
    
    def __init__(self, input_data, 
                 surface_type = 'Pt', 
                 asp = None, 
                 atoms = None, 
                 solver = 'glpk', 
                 distance_tune = 1.2,
                 H_dist = 2.75, 
                 ddec6 = None,
                 bond_cutoff = 0.3,
                 convert = True):
         
       
        """
        Attribute setter
        """
        for name, value in list(vars().items()):
                if name != 'self':
                    setattr(self,name,value)
        
        if isinstance(input_data, str):
            if 'DDEC6' in input_data or ddec6:
                
                self.full_graph = ddec6_to_graph(input_data)
                self.graph = graph_ops.remove_small_edges(self.full_graph)
                
                self.atoms = graph_to_atoms(self.graph, self.bond_cutoff)
                
            else:
                try:
                    self.atoms = read(input_data)
                    if convert: 
                        self.graph = atoms_to_graph(self.atoms)
                except Exception as e:
                    print("Can't read input")
                    print(e)
                    
        if isinstance(input_data, nx.classes.graph.Graph):
            self.graph = input_data
            if convert:
                self.atoms = graph_to_atoms(self.graph)
        
        self.elements = ['H','O','N','C']
        
        try:
            iter(self.surface_type)
            for i in self.surface_type:
                self.elements.append(i)
        except:
            if isinstance(self.surface_type, str):
                self.elements.append(self.surface_type)
            elif isinstance(self.surface_type, int):
                self.elements.append(chemical_symbols[self.surface_type])

        try:
            self.sort_atoms()
            
        except:
            pass

    def sort_atoms(self):
        """
        Sorts atoms by adsorbate, adsorbate coordinated, surface coordinated,
        and slab, based on the graph.
        """
        
        non_metals = [i for i in chemical_symbols if i in formula.non_metals]
        
        self.atoms_system = {'slab coordinated': [],
                             'slab': [],
                             'adsorbate': [],
                             'adsorbate coordinated': []
                             }

        
        for i in self.graph.nodes:
            
            node_neighbor_elements = [self.graph.nodes[j[1]]['element'] for j in self.graph.edges(i)]
            
            if self.graph.nodes[i]['element'] in self.surface_type:
                self.atoms_system['slab'].append(i)
                if any([j in non_metals for j in node_neighbor_elements]):
                    self.graph.nodes[i]['system'] = 'slab coordinated'
                    self.atoms_system['slab coordinated'].append(i)
                    
                else:
                    self.graph.nodes[i]['system'] = 'slab'
                
                self.graph.nodes[i]['cn'] = sum([j in self.surface_type for j in node_neighbor_elements])
                    
                    
            else:
                self.atoms_system['adsorbate'].append(i)
                if any([j in self.surface_type for j in node_neighbor_elements]):
                    self.graph.nodes[i]['system'] = 'adsorbate coordinated'
                    self.atoms_system['adsorbate coordinated'].append(i)
                else:
                    self.graph.nodes[i]['system'] = 'adsorbate'


    def center_adsorbate(self):
        
        # Routine to center the self.atoms object about the lowest adsorbate atom
        
        lowest = 100

        for i in self.atoms_system['adsorbate']:
            if self.graph.nodes[i]['position'][2] < lowest:
                lowest = self.graph.nodes[i]['position'][2]
                center_about = self.graph.nodes[i]['position'][2]
        
        self.atoms = self.center(self.atoms, center_about)
        
        for count, i in enumerate(self.atoms):
            self.graph.nodes['position'] = i.position
        
    def center(self, atoms, about, wrap = True):
        
        # In:   *atoms: ASE atoms object* to center; 
        #       *about: 3D array* centroid to center around; 
        #       *wrap: boolean" to wrap atoms back within cell
        
        # Out:  *atoms object* centered. 
        
        cell_center = sum(atoms.cell)/2
        
        move = cell_center - about
        
        atoms.positions += [move[0], move[1], 0]
        
        if wrap:    
            atoms.wrap()
        
        return atoms
    
    def get_gcn(self):
        
        max_cn = max(nx.get_node_attributes(self.graph, 'cn').values())
        
        second_nearest_neighbors = []
        for i in self.atoms_system['slab coordinated']:
            for j in self.graph.neighbors(i) :
                if self.graph.nodes[j]['system'] == 'slab':
                    second_nearest_neighbors.append(j)
                    
        second_nearest_neighbors = set(second_nearest_neighbors)
        sum_cn = sum([self.graph.nodes[i]['cn'] for i in second_nearest_neighbors])
        
        self.gcn = sum_cn / max_cn
        
        return self.gcn
        
    def gen_ads_subgraph(self, cluster_adjacency = 1, estimate_bond_index = True):
        
        
        
        nodes = set(i for i in self.graph.nodes if 'ads' in self.graph.nodes[i]['system'])
        adds = []
        
        for i in nodes:
            if self.graph.nodes[i]['system'] == 'adsorbate coordinated':
                adds.append(nx.ego_graph(self.graph, i, cluster_adjacency))
                
        for i in adds:
            nodes = nodes.union(set(i.nodes))
            
        self.ads_subgraph = graph_ops.assign_bond_index(self.graph.subgraph(nodes))
        

    def draw_graph(self, cutoff = 0.3, 
                        node_label_type = 'element', edge_labels = False, edge_label_type = 'weight'):
        
        from plotting import draw_surf_graph
        draw_surf_graph(self.graph, cutoff = cutoff, 
                            node_label_type = node_label_type, edge_labels = edge_labels, edge_label_type = edge_label_type)
        
    def generalize_graph_metal(self):
        if not self.graph:
            print('Generalizing metal error: No initial graph generated')
            return
        
        for n in self.graph.nodes(data=True):
            if self.graph.nodes[n[0]]['element'] == self.surface_type:
                self.graph.nodes[n[0]]['element']='M'
                
    
        
if __name__=="__main__":

    from plotting import draw_surf_graph as dsg
    from AtomsHelper import graph_ops
    

    #vasp_path = 'Examples/ase_to_graph/CONTCARmethane'
    #vasp_path = 'Examples/ase_to_graph/CONTCARalphahdown'
    #vasp_path = 'Examples/ase_to_graph/CONTCAR_OCH2CH3_hollow'
    #vasp_path = 'Examples/ase_to_graph/CONTCARmethylontop'
    vasp_path = 'Examples/ase_to_graph/CONTCARCCH2CH3ont_br'
    vasp_path = 'Examples/ase_to_graph/CONTCAR_CHCH_onon'
    
    ase_obj = read(vasp_path)
    ddec6_obj ='Examples/ddec6/CH2CHCH3CH3/DDEC6_even_tempered_net_atomic_charges.xyz'
    vasp = atoms_helper(ddec6_obj)
    # vasp = atoms_helper(vasp_path)
    vasp.sort_atoms()
    vasp.gen_ads_subgraph(cluster_adjacency = 1)
    # vasp = geom2graph(vasp = 'Examples/ase_to_graph/CONTCAR_CHCH_holhol')
    
    # ddec6 = gen_graph_ddec6('Examples/ddec6/CH2CHCH3CH3')
    # vasp = geom2graph(vasp = 'Examples/ase_to_graph/CONTCARpraneetalphahdown', surface_type = 'Pt')
    # vasp = geom2graph(vasp = 'Examples/ase_to_graph/CONTCAR_CHCH_onon', surface_type = 'Pt')
    
    # vasp.find_active_site_x_nearest_neighbors(2)
    
    # for i in vasp.active_site_nn[2]:
    #     vasp.atoms[i].symbol = 'Al'
    # for i in vasp.active_site_nn[1]:
    #     vasp.atoms[i].symbol = 'Au'
    # vasp.add_h_surf_bonds()
    # vasp.gen_graph(fill_bonds = True, verbose=True)
    #print(vasp.smiles_string)
    vasp.get_gcn()
    vasp.draw_graph(edge_labels=True, edge_label_type = 'bond_index', cutoff = False)
    dsg(vasp.ads_subgraph, edge_labels=True, edge_label_type = 'bond_index')
