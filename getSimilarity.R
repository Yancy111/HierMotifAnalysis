options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")

'''
# ensembl id -> entrez id
proteins = read.table('C:/Users/ZYX/OneDrive - zju.edu.cn/Graduation_Project/program_Louvain/data/proteins_ensembl_ids.txt', header=F)
library(biomaRt)
mart <- useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl')
# entrezids <- getBM(attributes = 'entrezgene_id', filters = 'ensembl_gene_id', values = proteins, mart = mart)
ids <- getBM(attributes = c('ensembl_gene_id','entrezgene_id'), filters = 'ensembl_gene_id', values = proteins, mart = mart)
write.table(ids, "C:/Users/ZYX/OneDrive - zju.edu.cn/Graduation_Project/program_Louvain/data/entrez_ids.txt", row.names=F, col.names=T, quote=F)
'''

# go similarity
library(BiocGenerics)
library(AnnotationDbi)
library(GOSemSim)
if (!requireNamespace("BiocManager", quietly = TRUE))
+     install.packages("BiocManager")
# BiocManager::install('org.Hs.eg.db')
library('org.Hs.eg.db')
hsGO <- godata('org.Hs.eg.db', ont='BP')
ICs <- hsGO@IC
write.csv(ICs, 'C:/Users/ZYX/OneDrive - zju.edu.cn/Graduation_Project/program_Louvain/data/go_IC.csv', row.names=TRUE, col.names=TRUE, quote=FALSE)
entrez = scan('C:/Users/ZYX/OneDrive - zju.edu.cn/Graduation_Project/program_Louvain/data/entrez_ids.txt', '')  # 注意数据类型！

# Rel_IEA_BMA
sim <- mgeneSim(entrez, hsGO, measure='Rel', drop='IEA', combine='BMA', verbose=TRUE)
write.csv(sim, 'C:/Users/ZYX/OneDrive - zju.edu.cn/Graduation_Project/program_Louvain/data/similarity_Rel_IEA_BMA.csv', row.names=TRUE, col.names=TRUE, quote=FALSE)

# Wang_IEA_BMA
# Resnik_IEA_BMA
# Jiang_IEA_BMA
# Lin_IEA_BMA

# Rel_IEA_avg
# Wang_IEA_avg
# Resnik_IEA_avg
# Jiang_IEA_avg
# Lin_IEA_avg