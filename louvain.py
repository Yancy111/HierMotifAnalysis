import os
import shutil
from collections import defaultdict

import json
import numpy as np
import networkx as nx
import community as community_louvain

'''
TODO best_partition: 选最后一层Or前面的层次；enrichment中太大的模块的处理
'''

def prepare_for_enrichment(partition, base_dir):
    com_ens = defaultdict(list)    # communities and ensembl ids in them
    for ens, com in zip(partition.keys(), partition.values()):
        com_ens[com].append(ens)
    print('This partition of the graph is consisted of', len(com_ens), 'clusters.')

    shutil.rmtree(base_dir+'results/enrichment')
    os.mkdir(base_dir+'results/enrichment')
    count = 0
    for com in com_ens:
        if len(com_ens[com]) >= 3:
            with open(base_dir+'results/enrichment/ensembls_in_'+str(com)+'.txt', 'w') as f:
                count += 1
                f.write('\n'.join(com_ens[com])+'\n')
                print('community', com, 'has', len(com_ens[com]), 'proteins in it.')

    print(count, 'communities have no less than 3 proteins in all.')

def prepare_for_topology(g, den, base_dir):

    def graph_to_txt(graphh):
        edges_list = []
        nodes = sorted(list(graphh.nodes()))
        node_index = dict([(i, str(nodes.index(i)+1)) for i in nodes])
        for edge in graphh.edges():
            a, b = edge[0], edge[1]
            edges_list.append(node_index[a]+'     '+node_index[b]+'     1')
        return '\n'.join(edges_list)+'\n'

    '''
    # 转换为numpy matrix
    # 先输出原网络
    print('The original graph has',g.number_of_nodes(), 'nodes and', g.number_of_edges(), 'edges.')
    with open(base_dir+'results/original_graph.npy', 'wb') as f:
        np.save(f, nx.to_numpy_matrix(g))
    for level in range(len(den)):
        partition = den[level]
        subgraph = community_louvain.induced_graph(partition, g)
        g = subgraph
        print('subgraph at level', level, 'has', subgraph.number_of_nodes(), 'nodes and', subgraph.number_of_edges(), 'edges.')
        with open(base_dir+'results/subgraph_of_'+str(level)+'.npy', 'wb') as f:
            np.save(f, nx.to_numpy_matrix(subgraph)) '''
    
    # 转换为txt用于modet
    # 先输出原网络
    print('The original graph has ',g.number_of_nodes(), '/', g.number_of_edges(), ' nodes / edges.')
    with open(base_dir+'results/original_graph.txt', 'w') as f:
        edges = graph_to_txt(g)
        f.write(edges)
    for level in range(len(den)):
        partition = den[level]
        subgraph = community_louvain.induced_graph(partition, g)
        g = subgraph
        print('subgraph at level', level, 'has', subgraph.number_of_nodes(), '/', subgraph.number_of_edges(), ' nodes / edges.')
        with open(base_dir+'results/subgraph_of_'+str(level)+'.txt', 'w') as f:
            edges = graph_to_txt(subgraph)
            f.write(edges)

def get_communities(G, base_dir, dendrogram=False, resolution=1.0):
    '''
    Generate dendrogram of the graph.
    --------------------------
    Input: networkx graph
    Output: files of 3 classes. (1) 'louvain_go_level': nodes and their corresponding clusters in different levels.
                                (2) 'cluster_nodes_level': nodes' ensembl ids in different clusters of different levels.
                                (3) 'subgraph_of_level': subgraph numpy matrix of different levels. (changed)
    --------------------------
    Example:
        gofilter = 'biological_process'
        graph = gengraph.gen_noweight_from_txt(file_dir, ens_dir, gofilter)
        dendrogram(graph, base_dir)
        Output:
            ------ For GO Annotations -----
            Partition at level 0 has 3725 nodes; There are 721 clusters at level 0.
            Partition at level 1 has 721 nodes; There are 248 clusters at level 1.
            Partition at level 2 has 248 nodes; There are 169 clusters at level 2.
            Partition at level 3 has 169 nodes; There are 158 clusters at level 3.
            ----- For Topology Analysis -----
            subgraph at level 0 has 721 nodes and 3354 edges.
            subgraph at level 1 has 248 nodes and 1050 edges.
            subgraph at level 2 has 169 nodes and 408 edges.
            subgraph at level 3 has 158 nodes and 269 edges.
    '''
    print('Performing louvain on the graph.')
    dendro = community_louvain.generate_dendrogram(G, resolution=resolution)
    partition = community_louvain.partition_at_level(dendro, len(dendro)-1) # best partition

    if dendrogram:  # 判断是否要输出dendrogram
        # for GO annotations
        print('------ Saving Results For GO Annotations -----')
        for level in range(len(dendro)):
            # get original nodes and their corresponding communities in different levels
            print('Partition at level', level, 'has', len(dendro[level]), 'nodes', end='; ')
            with open(base_dir+'results/louvain_go_'+str(level)+'.json','w') as f:
                json.dump(community_louvain.partition_at_level(dendro,level),f)

            # 从louvain_go文件中，得到不同层次cluster中相应的节点
            cluster = defaultdict(list)
            with open(base_dir+'results/louvain_go_'+str(level)+'.json') as f:
                nodes = json.load(f)
                for n in nodes.keys():
                    cluster[nodes[n]].append(n)
            
            print('There are', len(cluster), 'clusters at level '+str(level)+'.')
            with open(base_dir+'results/cluster_nodes_'+str(level)+'.json', 'w') as f:
                json.dump(cluster, f)

        # for topology analysis
        # networkx graph in different levels
        print('----- Saving Results For Topology Analysis -----')
        prepare_for_topology(G, dendro, base_dir)

    else:   # 不输出dendrogram就输出best partition
        print('Outputting files for later analysis.')
        # 后续用于modularity, T, entropy计算
        with open(base_dir+'results/best_partition.json', 'w') as f:
            json.dump(partition, f)
        # 后续用于富集分析
        prepare_for_enrichment(partition, base_dir)
        # 后续用于拓扑分析（最后分析，不然partition可能会改变）
        prepare_for_topology(G, dendro, base_dir)

    return partition