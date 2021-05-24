import math
from collections import defaultdict

import json
import numpy as np
import community as community_louvain

def getModularity(partition, G):
    '''
    Calculate graph's modularity based on its partition.
    ----------------------
    Input: partition and the networkx graph which is decomposed.
    Return: the modularity (float).
    '''
    return community_louvain.modularity(partition, G)

def getT(partition, G):
    '''
    Get graph's density T value.
    -----------------------
    Input: partition and the networkx graph which is decomposed.
    Return: T value.
    '''    
    '''
    # old version
    edges = defaultdict(int)  # 存储不同社区的边的数目
    for node in G.nodes():
        com = partition[node]
        for neighbor, datas in G[node].items():
            if partition[neighbor] == com:
                edges[com] += 1   
    inter_edges = sum(list(edges.values()))/2
    '''
    inter_edges = 0
    for edge in G.edges():
        if partition[edge[0]] == partition[edge[1]]:
            inter_edges += 1

    return inter_edges/len(G.edges())

def getEntropy(partition, G, ens_feature_dir):
    '''
    Calculate graph's entropy.
    -------------------------
    Input: partition, the networkx graph which is decomposed and protein feature directory.
    Return: entropy.
    '''
    def entropy(nodes, feat):
        ent = 0
        # 删除无相应label vector的ensembl id
        for n in nodes:
            if n not in feat:
                print(n, 'does not have feature vector.')
                del nodes[nodes.index(n)]
        # compute entropy
        arr = np.array([feat[n] for n in nodes])
        arr = arr.T
        for row in arr:
            nonzero = np.count_nonzero(row)
            if nonzero > 0:
                p = nonzero/len(list(row))
                ent -= p*math.log(p)   
        return ent

    N = G.number_of_nodes()
    communities = defaultdict(list)
    for node, com in partition.items():     # 得到每个社区内的ensembl ids
        communities[com].append(node)
    with open(ens_feature_dir) as f:
        features = json.load(f)
    
    count = 0       # GOSemSim用的是其他方法得到的ensembl-feature，因此存在部分蛋白质无feature的情况
    for node in partition.keys():
        if node not in features:
            count += 1
    if count:
        print('There are', count, 'nodes that do not have corresponding label vectors in feature file.')

    E = sum([len(communities[com])*entropy(communities[com], features)/N for com in communities])
    return E