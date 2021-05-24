options(BioC_mirror='http://mirrors.ustc.edu.cn/bioc/')

# ensembl id -> go id
library('biomaRt')
listMarts()
mart <- useMart('ensembl','hsapiens_gene_ensembl')
ensembl_ids <- read.table('C:/Users/ZYX/OneDrive - zju.edu.cn/Graduation_Project/program_Louvain/data/proteins_ensembl_ids_HI-II-14.txt')
goes <- getBM(attributes=c('ensembl_gene_id', 'go_id', 'definition_1006', 'namespace_1003'), filters='ensembl_gene_id', values=ensembl_ids, mart=mart)
write.csv(goes, 'C:/Users/ZYX/OneDrive - zju.edu.cn/Graduation_Project/program_Louvain/data/GO_annotations_HI-II-14.csv', row.names=FALSE, col.names=TRUE, quote=FALSE)