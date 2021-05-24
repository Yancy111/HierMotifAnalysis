# options('repos' = c(CRAN='https://mirrors.bfsu.edu.cn/CRAN'))
# options('repos' = c(CRAN='https://mirrors.tuna.tsinghua.edu.cn/CRAN/'))
options(BioC_mirror='http://mirrors.ustc.edu.cn/bioc/')
library(clusterProfiler)
library(BiocGenerics)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)

# 批量获取文件，读取到data列表中
path <- 'C:/Users/ZYX/OneDrive - zju.edu.cn/Graduation_Project/program_Integration/data/results/enrichment/'
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){paste(path,x,sep='/')})
data <- lapply(filePath, function(x){read.table(x)})
setwd('C:/Users/ZYX/OneDrive - zju.edu.cn/Graduation_Project/program_Integration/data/results/enrichment-result/')
file.remove(dir(".", pattern="(.jpeg)$") )  # 清空存储result的文件夹

# 存在的问题：1）有的entrez无对应kegg（if next）；2）有的无entrez（warning）；3）有的无富集结果（error）
for (i in 1:length(fileNames)){
    temp <- data[[i]]$V1
    temp_dir = gsub('.txt', '.jpeg', fileNames[i])   # 图片地址  TODO 这里需要两个[[]]吗
  
    tryCatch({  
        # 转换id: from ensembl to entrez
        entrez_lst <- bitr(temp, fromType='ENSEMBL', toType='ENTREZID', OrgDb='org.Hs.eg.db')
        entrez_lst = na.omit(entrez_lst)
        # GO analysis
        # go <- enrichGO(entrez_lst$ENTREZID, OrgDb='org.Hs.eg.db', ont='BP', pAdjustMethod='BH', pvalueCutoff=0.01, qvalueCutoff=0.05, keyType='ENTREZID', readable=TRUE) 
        # KEGG analysis
        kegg <- enrichKEGG(entrez_lst$ENTREZID, organism='hsa', pAdjustMethod='BH', pvalueCutoff=0.01, qvalueCutoff=0.05, keyType='kegg')
    }, warning=function(w){
        print('1st part warning')
        print(temp_dir)
    })
  
    tryCatch({
        # result display and output
        # goplot <- dotplot(go)
        # gg <- ggdraw() + draw_plot(goplot)
        com = gsub('.txt', '', fileNames[i])
        com = gsub('ensembls_in_', '', com)
        keplot <- dotplot(kegg, title=paste('The KEGG enrichment analysis result of community', com))
        gg <- ggdraw() + draw_plot(keplot)
        ggsave2(temp_dir, width=15, height=10)
    }, warning=function(w){
        print('2nd part warning')
        print(temp_dir)
    }, error=function(e){
        print('2nd part error')
        print(temp_dir)
    })
}


'''
# Reactome pathway analysis
# library(org.Hs.eg.db)
EG2Ensembl = toTable(org.Hs.egENSEMBL)
query = read.table('C:/Users/ZYX/OneDrive - zju.edu.cn/Graduation_Project/program_Integration/data/results/enrichment/ensembls_in_20.txt')
query = query$V1
geneLists = data.frame(ensembl_id=query)
results = merge(geneLists,EG2Ensembl,by='ensembl_id',all.x=T)
id = na.omit(results$gene_id)  #提取出非NA的ENTREZID
gene = id

library(ReactomePA)
reactome <- enrichPathway(gene=gene,pvalueCutoff=0.05, readable=T,organism = 'human')

barplot(reactome)
'''