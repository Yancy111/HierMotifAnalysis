from collections import defaultdict
import json
'''
备注都还没看
'''

def shortfilter(gofilter):
    if gofilter == 'biological_process':
        return '_bp'
    elif gofilter == 'molecular_function':
        return '_mf'
    elif gofilter == 'cellular_component':
        return '_cc'
    else:
        return '_others'

def extract_from_csv(fdir, outdir, gofilter, ic=False, ic_dir=None):
    '''
    Extract information from GO annotation csv file.
    ------------------
    Parameters:
        ic: elements of feature vectors are GO information content or not.
    TODO add multi-filter
    '''
    print('Extracting information from GO annotation csv file.')
    ens_go = defaultdict(set)     # dict: ensemble id and its corresponding go id
    go_id_def = {}  # dict: go id and its corresponding definition
    note = shortfilter(gofilter)

    # single filter
    with open(fdir) as f:
        for line in f.readlines():
            if  gofilter in line:
                go = [i.strip() for i in line[:line.find(','+gofilter)].split(',')] 
                ens_go[go[0]].add(go[1])
                go_id_def[go[1]] = ', '.join(go[2:])

    print('Number of proteins that have ' + gofilter + ' GO annotations:', len(ens_go))
    print('Number of GO ids:', len(go_id_def))

    # 确认ens_go中的go id在go def中都存在
    for key in ens_go.keys():
        ens_go[key] = list(ens_go[key])
    l = []
    for eg in ens_go.values():
        l.extend(eg)
    if len(set(l)) != len(go_id_def):
        print('[Warning] Not all the GO ids have definitions.')

    # 得到ensembl id 对应的特征向量
    ens_feature_vector = {}
    goids = sorted(list(go_id_def.keys()))  # 有顺序
    if not ic:
        for k in ens_go.keys():
            ens_feature_vector[k] = [0]*len(go_id_def)
            for g in ens_go[k]:
                ens_feature_vector[k][goids.index(g)] = 1
        if len(ens_feature_vector) != len(ens_go):
            print('[Warning] Not all the proteins have feature vectors.')
    else:
        # 没有ic的go id是否要删去（无影响）；评价标准是有无该特征，与元素的值无关（注意metric）
        goics = {}  # 存储GO ID及对应的ic
        exc_list = set()

        with open(ic_dir) as f:
            for line in f.readlines():
                g, c = line[:-1].split(',')
                if (g[:3] == 'GO:') and (c != 'Inf'):
                    goics[g] = float(c)
        
        for k in ens_go.keys():
            ens_feature_vector[k] = [0]*len(go_id_def)
            for g in ens_go[k]:
                if g in goics:
                    ens_feature_vector[k][goids.index(g)] = goics[g]
                else:
                    exc_list.add(g)
        print('There are', len(exc_list), 'GO IDs do not have valid information content.')        

        invalid_protein = 0
        keys = list(ens_feature_vector.keys())
        for efv in keys:  # 删除ens_go和ens_feature_vector中无有效go ic的蛋白质
            if not any(ens_feature_vector[efv]): # 字典的值为全零向量
                invalid_protein += 1
                del ens_feature_vector[efv]
                del ens_go[efv]
        if len(ens_feature_vector) != len(ens_go):
            print('[Warning] Ensembl-GO and Ensembl-Feature do not have the same proteins.')
        if invalid_protein:
            print('[Warning] Not all the proteins have feature vectors.', invalid_protein, 'proteins do not have IC feature.')

    # Output files
    with open(outdir+'ens_go'+note+'.json','w') as fe:
        json.dump(ens_go, fe)
    with open(outdir+'go_def'+note+'.json','w') as fg:
        json.dump(go_id_def, fg)
    with open(outdir+'ens_feature'+note+'.json','w') as ff:
        json.dump(ens_feature_vector, ff)
    print('Files have been saved in ' + outdir)

def extract_from_sim(sim_dir, ens_ent_dir):
    '''
    Extract information from GO similarity file which is produced by R package GOSemSim.
    Output: A dictionary whose value is dictionary too. 只有一半的结果（节省存储空间）.
    TODO figure out how the exception occur.
    '''
    prosim = defaultdict(dict)

    with open(ens_ent_dir) as f:
        fs = json.load(f)
        for s in fs:
            fs[s] = fs[s][0]    # 列表转为字符，且对应多个entrez id时只保留第一个；和getSimilarity中的输入id一致
        ent = dict(zip(fs.values(), fs.keys())) # 变成entrez-ensembl对应文件
        print('There are', len(ent), 'entrez ids have corresponding ensembl ids.')

    count = 0
    with open(sim_dir) as f:
        for line in f.readlines():
            if count == 0:  # title行
                count += 1
                pros = line[:-1].split(',') # pros[0] = ''
                print('There are', len(pros)-1, 'valid entrez ids in this GO similarity file.')
            else:
                values = line[:-1].split(',')
                for p in range(count+1, len(pros)):   # 只保留上三角矩阵（不包括对角线），减小存储需求
                    try:
                        prosim[ent[values[0]]][ent[pros[p]]] = float(values[p])
                        # prosim[(ent[values[0]], ent[pros[p]])] = float(values[p])
                    except KeyError:    # 会出现这个情况吗？
                        print('entrez', ent[pros[p]], 'does not have corresponding ensembl id.')
                count += 1

    outdir = sim_dir[:-4] + '.json'
    with open(outdir,'w') as f:
        json.dump(prosim, f)

