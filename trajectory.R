## huaqiang 2021.1.4
## Prepare data to run scVelo (RNA velocity) on sc-RNA dataset

#-----------------------------------------parameter-----------------------------------------
#npcs: 选定PCA降维的维数
#dims: 根据ElbowPlot(seurat_obj)  这一步的结果设定的参数
#resolution: 设置聚类的分辨率
#-----------------------------------------parameter-----------------------------------------





#--------------------------------设置工作路径---------------------------------------------------
setwd("/cluster/huanglab/hhuang/project/SingleCell-pipline/Trajectory/results/trajectory")
getwd()
rm(list = ls())
#--------------------------------设置工作路径---------------------------------------------------





#----------------------------------载入库的路径，里面有一些必须的包----------------------------------------
.libPaths(c('/cluster/huanglab/hhuang/app/envs/R3.6.2/lib/R/library/'))
#----------------------------------载入库的路径，里面有一些必须的包----------------------------------------




#-------------------------------------library essential package-----------------------------------------------------------------
library(Matrix)
library(Matrix.utils)
library(plyr)
library(dplyr)
library(Seurat)
library(sctransform)
library(igraph)
library(factoextra)
library(ComplexHeatmap)
library(circlize)
require(Hmisc)
require(dplyr)
require(openxlsx)
require(ggplot2)
library(ggpubr)
require(cowplot)
library(data.table)
library(RColorBrewer)
library(rowr)
library(SingleR)
library(scater)
library(pheatmap)
library(nichenetr)
library(tidyverse)
library(FlexDotPlot)
#-------------------------------------library essential package-----------------------------------------------------------------



# load data ---------------------------------------------------------------
# covid_combined.nc <- readRDS(url("https://hosted-matrices-prod.s3-us-west-2.amazonaws.com/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection-25/blish_covid.seu.rds")) ## 读入作者已经处理好的数据，我已经保存成Rdata了，跑下面的代码就行
# save(covid_combined.nc, file = 'covid_combined.nc.Rdata')


#-----------------------------------------------载入数据--------------------------------------------------------------
load("/cluster/huanglab/hhuang/project/SingleCell-pipline/Trajectory/results/trajectory/covid_combined.nc.Rdata")
#-----------------------------------------------载入数据--------------------------------------------------------------




#--------------------------------paper focus on Granulocyte and PB, so extract Granulocyte and PB cells----------------------
# Granulocyte and PB Analysis #

#--------------------------------------------一个处理数据的function-----------------------------------------
processNewSeurat <- function(parent.object, idents = NULL, cells = NULL, npcs = 50, dims = 50, resolution = 1) {
  if (!is.null(idents)) {
    seu <- subset(parent.object, idents = idents)
  }
  if (!is.null(cells)) {
    seu <- subset(parent.object, cells = cells)
  }
  message("Running PCA")
  seu <- RunPCA(seu, verbose = FALSE, npcs = npcs)
  message("Running UMAP")
  seu <- RunUMAP(seu, dims = 1:dims, verbose = FALSE)
  message("Clustering")
  seu <- FindNeighbors(seu, dims = 1:dims, verbose = FALSE)
  seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)
  return(seu)
}

covid_pb <- subset(covid_combined.nc, idents = c("9", "16", "18", "21", "27", "29"))
#--------------------------------paper focus on Granulocyte and PB, so extract Granulocyte and PB cells----------------------



#---------------------------------------------处理数据，PCA----UMAP----Cluster-------------------------------------------------
covid_pb <- processNewSeurat(covid_combined.nc, idents = c("9", "16", "18", "21", "27", "29"))
save(covid_pb, file = 'covid_pb.Rdata') ## 保存数据
#---------------------------------------------处理数据，PCA----UMAP----Cluster-------------------------------------------------


#---------------------------------------------在UMAP上可视化一下提取出来的两群细胞--------------------------------------------------
# DimPlot(covid_pb, label = TRUE) + NoLegend()
DimPlot(covid_pb, group.by = "cell.type.fine", label = TRUE) + NoLegend() + labs(x = "UMAP1", y = "UMAP2") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_color_manual(values = c("steelblue", "orange2", "green3", "red2"))
ggsave("p-celltype-fine-GPB.pdf", path = "./", width = 5, height = 5)
#---------------------------------------------在UMAP上可视化一下提取出来的两群细胞--------------------------------------------------








####---------------------------------RNA velocity, need spliced and unspliced counts to run RNA velocity
##---------------------------------------------extract data for scvelo
load("/cluster/huanglab/hhuang/project/SingleCell-pipline/Trajectory/results/covid_combined.emat.Rdata") ## spliced counts
load("/cluster/huanglab/hhuang/project/SingleCell-pipline/Trajectory/results/covid_combined.nmat.Rdata") ## unspliced counts
load("/cluster/huanglab/hhuang/project/SingleCell-pipline/Trajectory/results/trajectory/covid_pb.Rdata") ## cells we interest


#---------------------------------------这个数据的列名有点问题，需要修改一下--------------------------------
colnames(covid_combined.emat)
colnames(covid_pb)
length(colnames(covid_pb))
a = sub('.............$','',colnames(covid_combined.emat))
head(a)
colnames(covid_combined.emat) <- a
length(colnames(covid_combined.nmat))
colnames(covid_combined.nmat) <- a
length(intersect(a, colnames(covid_pb)))

save(covid_combined.emat, file = 'covid_combined.emat.Rdata')
save(covid_combined.nmat, file = 'covid_combined.nmat.Rdata')
#---------------------------------------这个数据的列名有点问题，需要修改一下--------------------------------



#-----------------------------------------------------保留spliced和unsplicedcounts都测到的基因----------------
gene.retain = intersect(rownames(covid_combined.emat), rownames(covid_combined.nmat))
#-----------------------------------------------------保留spliced和unsplicedcounts都测到的基因----------------



#---------------------------------------------写出一些必要的文件，用于scVelo的输入 ----------------------------------
spliced = covid_combined.emat[gene.retain, colnames(covid_pb)]
dim(spliced)
write.csv(t(spliced),file = 'spliced.csv',row.names = T)

unspliced = covid_combined.nmat[gene.retain, colnames(covid_pb)]
write.csv(t(unspliced),file = 'unspliced.csv',row.names = T)

metadata = covid_pb@meta.data
write.csv(metadata,file = 'metadata.csv',row.names = T)

write.csv(Embeddings(covid_pb, reduction = "umap"), file = "cell_embeddings.csv")
#---------------------------------------------写出一些必要的文件，用于scVelo的输入 ----------------------------------


### 之后转到python中进行RNA velocity的分析



