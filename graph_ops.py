# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 20:32:13 2024

@author: eccn3
"""

import copy
import networkx as nx

from ase.data import chemical_symbols
from ase import formula

import importlib
import sys

if importlib.util.find_spec('pyomo') is not None:
    import pyomo.environ as pe
    import pyomo.opt as po
    
def remove_small_edges(graph, edge_attribute = 'weight', cutoff = 0.3):
    small_weights = []

    for i in graph.edges():
        if graph.edges[i][edge_attribute] < cutoff:
            small_weights.append(i)
            
    copy_graph = copy.deepcopy(graph)
    for i in small_weights:
        copy_graph.remove_edge(*i)

    return copy_graph

def assign_bond_index(graph_original, valency_dict = {'H': 1,
                                             'O': 2,
                                             'N': 3,
                                             'C': 4,},
                      additional_elements: dict = None,
                      surface_element: dict = {'Pt': 4},
                      edge_attribute = 'bond_index',
                      verbose = False):
    
    graph = copy.deepcopy(graph_original)
    
    valency_dict.update(surface_element)
    
    if not additional_elements:
        
        non_metals = [i for i in chemical_symbols if i in formula.non_metals]
        
        graph_elements = set(nx.get_node_attributes(graph, 'element').values())
        
        for i in graph_elements.difference(valency_dict.keys()):
            if i not in non_metals:
                valency_dict[i] = 4
            else:
                print('No specified valency cap to new element.')
                return
            
    node_elements = nx.get_node_attributes(graph, 'element')
    
    for i in graph.nodes:
        graph.nodes[i]['valency'] = valency_dict[graph.nodes[i]['element']]
    
    for i in graph.edges:
        graph.edges[i][edge_attribute] = 0.
    
    for i in node_elements:
        if node_elements[i] == 'H':
            edges = len(graph.edges(i))
            for j in graph.edges(i):
                graph.edges[j][edge_attribute]= round(float(1/edges), 3)

    # Check if only hydrogen:
    

    if all([ i == 'H' or i == surface_element for i in nx.get_node_attributes(graph, 'element').values() ]):
        return graph


    # Check if only methane
    carbon_idx = [idx for idx, node in graph.nodes(data=True) if node['element']=='C']
    if len(carbon_idx) == 1:
        carbon_neighbors = graph.neighbors(carbon_idx[0])
        
        if all('H' in i for i in [node['element'] for idx, node in graph.subgraph(list(carbon_neighbors)).nodes(data=True)]):
            return graph
        
      
    if 'pyomo.environ' not in sys.modules:
        print('No Pyomo package found or imported. Returning graph')
        return graph
    
    model, results = fill_bonds(graph, 
                                surface_element = list(surface_element.keys())[0], 
                                edge_attribute = edge_attribute,
                                verbose = verbose)
    
    if not model:
        return
    
    solver_results = {'model': model, 'results': results}
    
    if solver_results['model']:
        
        for index in model.edge_weights_ads:
            graph.edges[index][edge_attribute] = float(model.edge_weights_ads[index].value)
            
        if list(surface_element.keys())[0] in [graph.nodes[i]['element'] for i in graph.nodes]:
            
            edges_pt_dict = {i: float(model.edge_weights_pts[i].value) for i in model.edge_weights_pts}          
            nodes = set([i for sub in model.edges_pts for i in sub])
            not_pt_nodes = [i for i in nodes if graph.nodes[i]['element'] != surface_element]
            
            edge_pt_ads_sort = dict.fromkeys(not_pt_nodes, {})
            
            for index in model.edge_weights_pts:
                intersection = set(index).intersection(not_pt_nodes)
                if intersection:
                    edge_pt_ads_sort[list(intersection)[0]][index] = float(model.edge_weights_pts[index].value)
                    
            for node in edge_pt_ads_sort:
                
                edge_sum = 0
                edges= 0
                
                for edge in edge_pt_ads_sort[node]:
                    edge_sum += edge_pt_ads_sort[node][edge]
                    edges += 1
                edge_sum = edge_sum/edges
                
                for edge in edge_pt_ads_sort[node]:
                    edges_pt_dict[edge] = edge_sum
            
            for index in edges_pt_dict: 
               graph.edges[index][edge_attribute] = round(edges_pt_dict[index],3)
                
    else:
        print('Problem with solver')

    return graph

def fill_bonds(graph, 
               surface_element = 'Pt', 
               edge_attribute = 'bond_index',
               verbose = False):
    
    node_dict = dict(graph.nodes)

    edges_ads, edges_pts, edges_hs, nodes_ads, nodes_pts, nodes_hs= {}, {}, {}, {}, {}, {}
    
    # Sort edges into surface and adsorbate
    
    for edge in graph.edges:
        
        edge_elements = [graph.nodes[edge[0]]['element'], graph.nodes[edge[1]]['element']]
        
        if 'H' in edge_elements:
            edges_hs[edge] = graph.edges[edge][edge_attribute]
            
        elif edge_elements == [surface_element, surface_element]:
            continue
        
        elif surface_element in edge_elements:
            edges_pts[edge] = graph.edges[edge][edge_attribute]
                
        else:
            edges_ads[edge] = graph.edges[edge][edge_attribute]
    
    # Sort nodes into surface and adsorbate
    for node in graph.nodes:
        
        if graph.nodes[node]['element'] == surface_element:
            nodes_pts[node] = node_dict[node]['valency']
            
        elif graph.nodes[node]['element'] == 'H':
            nodes_hs[node] = node_dict[node]['valency']
            
        else:
            nodes_ads[node] = node_dict[node]['valency']
            
    if verbose:
        print('Solving bond indices...')

    if not edges_ads and not edges_pts:
        print('No adsorbate/slab edges; trivial solution returned')
        return None, None
    
    model = pe.ConcreteModel()

    model.nodes_ads = pe.Set(initialize = list(nodes_ads.keys()))
    model.nodes_hs = pe.Set(initialize = list(nodes_hs.keys()))
    model.edges_ads = pe.Set(initialize = list(edges_ads.keys()))
    model.edges_hs = pe.Set(initialize = list(edges_hs.keys()))
    
    model.edge_weights_ads = pe.Var(model.edges_ads, 
                                    within = pe.PositiveIntegers,
                                    initialize = {i: 3 for i in model.edges_ads},
                                    bounds = (1, 4))
    
    model.edge_weights_hs = pe.Param(model.edges_hs,
                                     initialize = {i: 1 for i in model.edges_hs})
    
    model.valency = pe.Param(model.nodes_ads, 
                             initialize = nodes_ads, 
                             within=pe.PositiveIntegers)
    
    
    if len(nodes_pts)>0:
        
        model.nodes_pts = pe.Set(initialize = list(set(nodes_pts.keys())))
        model.edges_pts = pe.Set(initialize = list(set(edges_pts)))
        model.edge_weights_pts = pe.Var(model.edges_pts, 
                                        within = pe.PositiveReals,
                                        initialize ={i: 0.5 for i in model.edges_pts},
                                        bounds = (0.01, 4))
    
        def valency_cap(model, i):
            
            edge_sum = 0
            
            if edges_ads:
                edge_sum += sum([model.edge_weights_ads[j] for j in model.edge_weights_ads if i in j]) 
            if edges_hs:
                edge_sum += sum([model.edge_weights_hs[j] for j in model.edge_weights_hs if i in j])
            if edges_pts:
                edge_sum += sum([model.edge_weights_pts[j] for j in model.edge_weights_pts if i in j])
            
            return model.valency[i] == edge_sum
                   
    
        def edge_sum(model):
            
            total_ads = sum([(model.edge_weights_ads[i]) for i in model.edges_ads])
            total_pts = sum([(model.edge_weights_pts[i]) for i in model.edges_pts])*0.001
            
            total = total_ads + total_pts
            
            return total
        
    else:
        
        def valency_cap(model, i):
        
            node_in_edge_ads = sum([model.edge_weights_ads[j] for j in model.edge_weights_ads if i in j])
            node_in_edge_hs = sum([model.edge_weights_hs[j] for j in model.edge_weights_hs if i in j])
            
            return model.valency[i] == node_in_edge_ads + node_in_edge_hs
                   
        
        def edge_sum(model):
            
            return sum([(model.edge_weights_ads[i]) for i in model.edges_ads])



    model.obj = pe.Objective(sense = pe.maximize, rule = edge_sum)

    constraint_nodes = set([i for sub in model.edges_ads for i in sub])
    
    if edges_pts:
        constraint_nodes_pts = set([i for sub in model.edges_pts for i in sub])
        constraint_nodes = [i for i in constraint_nodes.union(constraint_nodes_pts) if graph.nodes[i]['element'] != surface_element]
    
    model.con = pe.Constraint(constraint_nodes, rule = valency_cap)

    solver = po.SolverFactory('glpk')
    
    try:
        results = solver.solve(model)
            
    except:
        if verbose:
            print('Solver unable to get solution')
        model = None
        results = None
    
    return model, results