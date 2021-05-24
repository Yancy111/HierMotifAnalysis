import math
from collections import defaultdict

import json
import numpy as np
import networkx as nx

def shortgofilter(gofilter):
    if gofilter == 'biological_process':
        return '_bp'
    elif gofilter == 'molecular_function':
        return '_mf'
    elif gofilter == 'cellular_component':
        return '_cc'
    else:
        return '_others'

def similarity(a, b, method='multi-cosine', gofilter='biological_process'):
    '''
    Calculate nodes' similarity.
    ---------------------------
    Input: nodes' feature vector
    Return: nodes' GO similarity
    ---------------------------
    '''

    if method == 'multi-cosine':
        den1 = math.sqrt(sum([pow(i,2) for i in a]))
        den2 = math.sqrt(sum([pow(i,2) for i in b]))
        s = sum([a[i]*b[i] for i in range(len(a))])
        if den1*den2 != 0:
            result = s/(den1*den2)
        else:
            result = 0

    elif method == 'jaccard':
        both, single =  0, 0
        for i in range(len(a)):
            if a[i] + b[i] == 2:
                both += 1
            elif a[i] + b[i] == 1:
                single += 1
        if (both + single) != 0:
            result = both / (both + single)
        else:
            result = 0
    
    elif method == 'multi-cosine-GOIC':
        den1 = math.sqrt(sum([pow(i,2) for i in a]))
        den2 = math.sqrt(sum([pow(i,2) for i in b]))
        s = sum([a[i]*b[i] for i in range(len(a))])
        if den1*den2 != 0:
            result = s/(den1*den2)
        else:
            result = 0

    elif method == 'jaccard-GOIC':
        both, single = 0, 0
        for i in range(len(a)):
            if (a[i] != 0) and (b[i] != 0) and (a[i] == b[i]):
                both += a[i]
            elif (a[i] != 0) and (b[i] == 0.0):
                single += a[i]
            elif (a[i] == 0.0) and (b[i] != 0):
                single += b[i]
        if (both + single) != 0:
            result = both / (both + single)
        else:
            result = 0

    return result

def gen_noweightgraph_from_txt(fdir, ens_dir=None, gofilter=None):
    '''
    Generate non-weight graph with gofilters.
    '''
    G = nx.Graph()
    with open(fdir) as ff:
        for interactors in ff.readlines():
            interactors = interactors[:-1].split()
            if interactors[0] != interactors[1]:    # self-loops are not included
                G.add_edge(interactors[0], interactors[1])   
    a, b = G.number_of_nodes(), G.number_of_edges()
    print('Generate a no-self-loop unweighted graph with', a, 'nodes and', b, 'edges.')

    if gofilter == None:
        return G
    if gofilter != None and ens_dir == None:
        print('Invalid! Please input the valid protein ensembl file.')
        return

    if gofilter != None and ens_dir != None:
        note = shortgofilter(gofilter)
        if note not in ens_dir:
            print('Invalid! Please input the valid protein ensembl file.')
            return

        with open(ens_dir) as fe: 
            ens_go = json.load(fe)
        nodes = list(G.nodes())
        for n in nodes:     # remove proteins that do not have corresponding GO annotations
            if n not in ens_go:
                G.remove_node(n)

        nodes = list(G.nodes())
        for n in nodes:     # remove single node
            if not G[n]:
                G.remove_node(n)  

        print(a-G.number_of_nodes(), 'nodes that do not have', gofilter, 'GO annotations and', b-G.number_of_edges(), 'edges have been removed.')
        print('Generate a no-self-loop unweighted graph with', G.number_of_nodes(), 'nodes and', G.number_of_edges(), 'edges.')

        return G

def gen_weightedgraph_from_txt(fdir, ens_dir=None, ens_feature=None, gofilter=None, method='multi-cosine'):
    '''
    Generate a weighted graph with gofilters.
    Weight can be either nodes' similarity or 1. If there is no gofilter or no feature file, then weight is assigned to 1.
    '''
    G = gen_noweightgraph_from_txt(fdir, ens_dir=ens_dir, gofilter=gofilter)
    graph = nx.Graph()
    weighted_list = []  # save weighted edges
            
    if gofilter == None:
        for edge in G.edges():
            weighted_list.append((edge[0], edge[1], 1))

    if gofilter != None and ens_feature == None:
        print('Invalid! Please input the valid protein feature file.')
        return
    
    if gofilter != None and ens_feature != None:
        note = shortgofilter(gofilter)
        if note not in ens_feature:
            print('Invalid! Please input the valid protein feature file.')
            return 
        if method == None:
            print('Invalid! Please enter a valid method.')
            return

        if 'rank' in method:    # e.g. 'multi-cosine+rank'
            print('METHOD: using similarity rank & softmax value as weight')
            method = method[:-5]    # 避免计算相似度时有问题
            with open(ens_feature) as ff:
                features = json.load(ff)
                weights = []
                for edge in G.edges():
                    s0, s1 = features[edge[0]], features[edge[1]]
                    s = similarity(s0, s1, method=method)
                    weights.append((edge[0], edge[1], s))

            weights.sort(key=lambda w:w[2])  # 根据相似度大小排序     
            s_total = sum([weights.index(i) for i in weights]) + len(weights)
            # s_total = (1 + len(weights)) * len(weights) /2 + len(weights)
            for w in weights:
                weighted_list.append((w[0], w[1], (weights.index(w)+1)/s_total))
                # weighted_list.append((w[0], w[1], weights.index(w)))
        else:
            print('METHOD: using similarity value as weight')
            with open(ens_feature) as ff:
                features = json.load(ff)
                for edge in G.edges():
                    s0, s1 = features[edge[0]], features[edge[1]]
                    s = similarity(s0, s1, method=method)
                    weighted_list.append((edge[0], edge[1], s))

    graph.add_weighted_edges_from(weighted_list)
    print('Generate a no-self-loop weighted graph with', graph.number_of_nodes(), 'nodes and', graph.number_of_edges(), 'edges.')
    return graph

def gen_noweightgraph_gosim(fdir, prosim_dir=None, gofilter=None):
    '''
    Generate a no-weight graph according to GO similarity file which is produced by R package GOSemSim.
    --------------------
    Input:
    Return: a networkx graph of no weight
    '''
    G = nx.Graph()
    with open(fdir) as ff:
        for interactors in ff.readlines():
            interactors = interactors[:-1].split()
            if interactors[0] != interactors[1]:    # self-loops are not included
                G.add_edge(interactors[0], interactors[1])   
    a, b = G.number_of_nodes(), G.number_of_edges()
    print('Generate a no-self-loop unweighted graph with', a, 'nodes and', b, 'edges.')

    if gofilter == None:
        return G
    if gofilter != None and prosim_dir == None:
        print('Invalid! Please input the valid protein GO similarity file.')
        return

    if gofilter != None and prosim_dir != None:
        note = shortgofilter(gofilter)
        if note not in prosim_dir:
            print('Invalid! Please input the valid protein GO similarity file.')
            return

        with open(prosim_dir) as fe: 
            prosim = json.load(fe)  # 注意这里只有一半的关系
            pro_go = list(prosim.keys())
        nodes = list(G.nodes())
        for n in nodes:     # remove proteins that do not have corresponding GO annotations
            if n not in pro_go:
                G.remove_node(n)

        nodes = list(G.nodes())
        for n in nodes:     # remove single node
            if not G[n]:
                G.remove_node(n)

        print(a-G.number_of_nodes(), 'nodes that do not have', gofilter, 'GO annotations and', b-G.number_of_edges(), 'edges have been removed.')
        print('Generate a no-self-loop unweighted graph with', G.number_of_nodes(), 'nodes and', G.number_of_edges(), 'edges.')

        return G

def gen_weightedgraph_gosim(fdir, prosim_dir=None, gofilter=None, method='GOSemSim'):
    '''
    Generate a weighted graph according to GO similarity file which is produced by R package GOSemSim.
    ---------------------
    Input:
    Return: a weighted networkx graph  
    '''
    G = gen_noweightgraph_gosim(fdir, prosim_dir=prosim_dir, gofilter=gofilter)
    graph = nx.Graph()
    weighted_list = []

    if 'rank' in method:    # e.g. 'GOSemSim+rank'
        print('METHOD: using GOSemSim value rank & softmax valueas weight')
        with open(prosim_dir) as ff:
            prosim = json.load(ff)
            weights = []
            for edge in G.edges():
                a, b = edge[0], edge[1]
                if (a in prosim) and (b in prosim[a]):
                    s = prosim[a][b]
                elif (b in prosim) and (a in prosim[b]):
                    s = prosim[b][a]
                weights.append((edge[0], edge[1], s))

            weights.sort(key=lambda w:w[2])     
            s_total = sum([weights.index(i) for i in weights]) + len(weights)
            # s_total = (1 + len(weights)) * len(weights) /2 + len(weights)
            for w in weights:
                weighted_list.append((w[0], w[1], (weights.index(w)+1)/s_total))
                # weighted_list.append((w[0], w[1], weights.index(w)))
    else:
        print('METHOD: using GOSemSim value as weight')
        with open(prosim_dir) as ff:
            prosim = json.load(ff)
            for edge in G.edges():
                a, b = edge[0], edge[1]
                if (a in prosim) and (b in prosim[a]):
                    s = prosim[a][b]
                elif (b in prosim) and (a in prosim[b]):
                    s = prosim[b][a]
                weighted_list.append((a, b, s))

    graph.add_weighted_edges_from(weighted_list)
    print('Generate a no-self-loop weighted graph with', graph.number_of_nodes(), 'nodes and', graph.number_of_edges(), 'edges.')
    return graph