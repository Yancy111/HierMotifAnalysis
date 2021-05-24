from collections import defaultdict
import json

'''
Convert ensembl ids to entrez ids. Delete duplicated and invalid (NA) ones.
'''


with open('./data/entrez_ids_original_HI-II-14.txt') as f:
    ids = defaultdict(list)
    for line in f.readlines():
        num, ensembl, entrez = line[:-1].split()
        ids[ensembl].append(entrez)

keys = list(ids.keys())
for k in keys:
    if 'NA' in ids[k]:
        del ids[k]

keys = sorted(list(ids.keys()))
ent = []
for i in range(len(keys)):
    # ent.append("\t".join([str(i), keys[i], ids[keys[i]][0]]))
    ent.append(ids[keys[i]][0])
with open('./data/entrez_ids_HI-II-14.txt', 'w') as f:   # 该文件中无重复
    f.write("\n".join(ent))
with open('./data/ensembl-entrez_HI-II-14.json','w') as f:   # 该文件中包含重复，extract时会进行处理
    json.dump(ids, f)

'''
# 确认无重复
for k in ids:
    if len(set(ids[k])) != len(ids[k]):
        print(k)
'''
'''
# 找有多个entrez id 的 ensembl id
for k in ids:
    if len(ids[k])>1:
        print(k, len(ids[k]))
# 结果
ENSG00000076928 2
ENSG00000105793 2
ENSG00000120341 2
ENSG00000120709 2
ENSG00000124713 2
ENSG00000137843 2
ENSG00000137936 2
ENSG00000143702 2
ENSG00000145979 2
ENSG00000146112 2
ENSG00000158301 2
ENSG00000158747 2
ENSG00000163156 2
ENSG00000166272 2
ENSG00000178882 2
ENSG00000181135 2
ENSG00000188629 2
ENSG00000196605 2
ENSG00000204209 2
ENSG00000205571 2
ENSG00000214026 2
ENSG00000215269 5
ENSG00000224659 2
ENSG00000236362 7
ENSG00000268606 2
ENSG00000276070 2
'''
'''
# 找无entrez id 的 ensembl id
for k in ids:
    if 'NA' in ids[k]:
        print(k, ids[k])
# 结果
ENSG00000064489 ['NA']
ENSG00000158483 ['NA']
ENSG00000171570 ['NA']
ENSG00000183889 ['NA']
ENSG00000213204 ['NA']
ENSG00000233024 ['NA']
ENSG00000255104 ['NA']
ENSG00000257390 ['NA']
ENSG00000259288 ['NA']
ENSG00000259529 ['NA']
ENSG00000268173 ['NA']
ENSG00000268500 ['NA']
ENSG00000270136 ['NA']
ENSG00000272617 ['NA']
ENSG00000284341 ['NA']
'''