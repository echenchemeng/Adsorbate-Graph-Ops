# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 18:14:20 2024

@author: eccn3
"""

import matplotlib.pyplot as plt
import networkx as nx
from ase.data.colors import jmol_colors
from ase.data import chemical_symbols
from ase import formula
import numpy as np
import copy

def draw_surf_graph(graph, cutoff = 0.3, 
                    node_label_type = 'element', edge_labels = False, edge_label_type = 'weight'):
   
    '''
    Draws graph of atoms system
    Input:
        graph: networkx graph of atoms
        cutoff: min edge value to be drawn
        node_label_type: what attribute to print as node label
        edge_labels: boolean to draw edge labels
        edge_label_type: what attribute to print as edge label
    Output:
        prints graph
    '''
    
    if edge_label_type !='weight' and edge_labels == False:
        edge_labels=True
    
    if cutoff:
        try:
            if nx.get_edge_attributes(graph, 'weight'):
                graph = remove_small_edges(graph, cutoff = cutoff)
        except:
                pass    

    plt.figure(figsize=(6,5))
   
    colors = {i:j.reshape(1,-1) for i,j in zip(chemical_symbols, jmol_colors)}  

    metals = [i for i in chemical_symbols if i not in formula.non_metals]
    metals.append('M')
    
    colors['M'] = '#9c58a1'

    # for i in colors.keys():
    #     colors[i] = np.array(colors[i]).reshape(1,-1)              

    element_count = {i: [] for i in colors.keys()}

    for i in nx.get_node_attributes(graph, 'element').items():
        element_count[i[1]].append(i[0])
    element_count = {k: v for k,v in element_count.items() if v}
    init_pos = nx.get_node_attributes(graph, 'element')
   
    for i in init_pos.keys():
        if any([init_pos[i] ==j for j in metals]):
            init_pos[i] = [np.random.uniform(),0.05 +0.1* np.random.uniform()]
        else:
            init_pos[i] = [np.random.uniform(), 0.85 +0.1* np.random.uniform()]
   
    pos = nx.spring_layout(graph,pos=init_pos)    
       
    system_count = {'metal surface coordinated': [], 'adsorbate surface coordinated': []}
    keys = system_count.keys()
    for i in nx.get_node_attributes(graph, 'system_type').items():
        if any([i[1] == j for j in keys]):
            system_count[i[1]].append(i[0])
           
           
    options = {"edgecolors": "tab:cyan"}
    for i in system_count.keys():    
        nx.draw_networkx_nodes(graph, pos, nodelist = system_count[i],node_size=500, **options)
       
    for i in element_count.keys():
    
        if i in metals:
            nx.draw_networkx_nodes(graph, pos, nodelist=element_count[i],
                                   node_size=200, node_color = colors[i],
                                   edgecolors = (0,0,0), alpha = 0.5)
        else:
            nx.draw_networkx_nodes(graph, pos, nodelist=element_count[i],
                                   node_size=400, node_color = colors[i],
                                   edgecolors=(0,0,0), alpha = 1)

       
       
    nx.draw_networkx_edges(graph, pos, alpha = 0.5)
   
    labels = nx.get_node_attributes(graph, node_label_type)
    
    nx.draw_networkx_labels(graph, pos, labels)
    
    if edge_labels:
        ed_labels = nx.get_edge_attributes(graph, edge_label_type)
        for i in ed_labels:
            ed_labels[i] = round(ed_labels[i], 2)

        nx.draw_networkx_edge_labels(graph, pos, ed_labels)
    # plt.close()

def remove_small_edges(graph, cutoff = 0.3):
    small_weights = []

    for i in graph.edges():
        if graph.edges[i]['weight'] < cutoff:
            small_weights.append(i)
            
    copy_graph = copy.deepcopy(graph)
    for i in small_weights:
        copy_graph.remove_edge(*i)

    return copy_graph
