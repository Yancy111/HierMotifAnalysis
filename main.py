import extract
import gengraph
import louvain
import metric

'''
TODO louvain中选择最后一层还是倒数第二层作为best_partition
'''

def main(base_dir, interaction_dir, weight=True, extraction=False, gofilter='biological_process', method='multi-cosine', ic_dir=None, sim_dir=None, ens_ent_dir=None, prosim_dir=None, dendrogram=False, resolution=1.0):
    print('[Parameters] WEIGHT:', weight, 'FILTER:', gofilter, 'METHOD:', method)

    # extract information with different gofilters from 'protein GO annotation file' or 'GOSemSim result file (protein functional similarity)' and produce protein features
    # generate the networkx graph
    if 'GOSemSim' in method:    # 调用GOSemSim相关程序
        if gofilter != None:
            note = gengraph.shortgofilter(gofilter)
            ens_feature_dir = base_dir+'ens_feature'+note+'.json'   # 从非GOSemSim中得到的ens_feature
        else:
            prosim_dir = None

        if extraction:
            extract.extract_from_sim(sim_dir, ens_ent_dir)
            prosim_dir = sim_dir[:-4] + '.json' # 避免未输入prosim_dir的情况或输入的prosim_dir未更新的情况；假如不extract，则利用现成的prosim_dir，假如extract，则现成提供的prosim_dir会更新。
        if weight:
            G = gengraph.gen_weightedgraph_gosim(interaction_dir, prosim_dir, gofilter, method)
        else:
            G = gengraph.gen_noweightgraph_gosim(interaction_dir, prosim_dir, gofilter)
    else:       # 调用one-hot code相关程序
        if gofilter != None:
            note = gengraph.shortgofilter(gofilter)
            ens_dir = base_dir+'ens_go'+note+'.json'
            ens_feature_dir = base_dir+'ens_feature'+note+'.json'
        else:
            ens_dir = None

        if extraction:
            if gofilter == None:
                print('Cannot perform extraction without filters.')
                return
            if 'GOIC' in method:
                extract.extract_from_csv(gofile_dir, base_dir, gofilter, ic=True, ic_dir=ic_dir)
            else:
                extract.extract_from_csv(gofile_dir, base_dir, gofilter)
        if weight:
            G = gengraph.gen_weightedgraph_from_txt(interaction_dir, ens_dir, ens_feature_dir, gofilter, method)
        else:
            G = gengraph.gen_noweightgraph_from_txt(interaction_dir, ens_dir, gofilter)

    # perform louvain on the graph
    best_partition = louvain.get_communities(G, base_dir, dendrogram=dendrogram, resolution=resolution)
    
    # evaluate louvain results: modularity, T value, entropy; prepare for enrichment analysis and topological analysis as well
    print('Modularity of the best partition:', metric.getModularity(best_partition, G))
    print('T of the best partition:', metric.getT(best_partition, G))
    print('E of the best partition:', metric.getEntropy(best_partition, G, ens_feature_dir))

if __name__ == '__main__':
    dataset = 'HI-II-14'
    base_dir = './data/'
    interaction_dir = base_dir + dataset + '.txt'
    gofile_dir = base_dir + 'GO_annotations_' + dataset + '.csv'
    goic_dir = base_dir + 'GO_IC.csv'
    similarity_dir = base_dir + 'similarity_Rel_IEA_BMA_bp.csv'
    ensembl_entrez_dir = base_dir + 'ensembl-entrez.json'
    prosim_dir = base_dir + 'similarity_Rel_IEA_BMA_bp.json'
    
    weight = True
    extraction = False
    gofilter = 'biological_process'    # 'biological_process', 'molecular_function'
    method = 'jaccard-GOIC+rank' # 'multi-cosine', 'jaccard', 'multi-cosine-GOIC', 'jaccard-GOIC', 'GOSemSim'
                            # 'multi-cosine+rank', 'jaccard+rank', 'multi-cosine-GOIC+rank', 'jaccard-GOIC+rank', 'GOSemSim+rank'
                            # 'GOSemSim'在'multi-cosine'而非'GOIC'之后进行
    dendrogram = False
    resolution = 1.0

    main(base_dir, interaction_dir, weight=weight, extraction=extraction, gofilter=gofilter, method=method, 
        ic_dir=goic_dir, sim_dir=similarity_dir, ens_ent_dir=ensembl_entrez_dir, prosim_dir=prosim_dir, 
        dendrogram=dendrogram, resolution=resolution)
    
    # topology analysis commond line:
    # cd C:\Users\ZYX\OneDrive - zju.edu.cn\Graduation_Project\program_Integration\src\modet
    # run -m modet -s 3 -g 

