library(factoextra)
library(abind)
library(viridis)
library(ggplot2)
library(NMF)
library(dbscan)
library(umap)
library(data.table)
library(pheatmap)
library(clusterSim)
###################################################################################################################
# Allele Clustering: Indentifying Populations of Alleles via Distance Matrices from ORCA experiments.
# 09-15-23
# Noah Burget
###################################################################################################################
setwd("/mnt/data0/noah/analysis/ORCA_analysis/allele_clustering")
source("processChrTracer.R")
source("plotContactFrequency.R")
###################################################################################################################
#                                         Current dataset being processed:
#   230902_Granta519cl27_24hdTAG_MYC5p_30mHyb_4phBl_30step_allfits.csv
str <- '230902_Granta519cl27_24hdTAG_MYC5p_30mHyb_4phBl_30step_allfits'
raw <- "ChTracer3_output/230902_Granta519cl27_24hdTAG_MYC5p_30mHyb_4phBl_30step_allfits.csv"
Forward_all_distances_adj.DTag <- process.chrtracer(raw, 10, T)
full_dist.DTag <- make.distTable(Forward_all_distances_adj.DTag)
dist.upper.lens <- unlist(lapply(Forward_all_distances_adj.DTag, function(x) length(which(!is.na(x)))))
###################################################################################################################
# Form of input data:
# TRACES/ROWS X DISTANCES/COLS

full_dist.noDTag <- as.data.frame(full_dist.noDTag)
##### 0: KMeans on distance table
fviz_nbclust(full_dist.noDTag, kmeans, "wss")
kmean.res <- kmeans(full_dist.noDTag, 3)
full_dist.noDTag$kmeans.cluster <- kmean.res$cluster

##### 1.1: Nonnegative Matrix Factorization (NMF)
set.seed(123)
# Estimate rank for input
nmf.res <- nmf(x = as.matrix(full_dist.DTag), 
               rank = 2, 
               method = 'lee',
               seed = 'random')
# Extract Basis matrix 
nmf.basis <- as.data.frame(basis(nmf.res))
# Plot basis matrix
ggplot(data=nmf.basis, aes(x = V1, y = V2))+
  geom_point()
#pdf(paste0(str,'snmf-l_basisHeatmap.pdf'))
#basismap(nmf.res[["snmf/l"]])
#dev.off()

##### 1.1.a: PCA
set.seed(123)
pca.res <- as.data.frame(prcomp(full_dist.DTag, center=F, scale=T)$x[,1:2])
# Plot first 2 PCs
ggplot(data = pca.res, aes(x = PC1, y = PC2))+
  geom_point()

##### 1.2: KMeans clustering of NMF basis matrix
set.seed(123)
# Estimate optimal k
fviz_nbclust(nmf.basis, kmeans, method='silhouette')
# KMeans clustering
kmeans.res <- kmeans(nmf.basis,3)
nmf.basis$kmeans.cluster <- kmeans.res$cluster
# Plot basis matrix by cluster
ggplot(data = nmf.basis, aes(x = V1, y = V2))+
  geom_point(aes(color = factor(kmeans.cluster)))
# Use factoextra package for better cluster viz
fviz_cluster(kmeans.res,
             data = nmf.basis[,-3],
             nmf.basispalette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw())

##### 1.3: KMeans clustering of first 2 PCs
set.seed(123)
# Estimate optimal k
fviz_nbclust(pca.res, kmeans, method='wss')
kmeans.res <- kmeans(as.matrix(pca.res), 3)
fviz_cluster(kmeans.res,
             data = pca.res,
             nmf.basispalette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw())
pca.res$kmeans.cluster <- kmeans.res$cluster

##### 2: UMAP
umap.res <- as.data.frame(umap(full_dist.noDTag)$layout)
# Plot UMAP embedding
ggplot(data = umap.res, aes(x = V1, y = V2))+
  geom_point()
##### 2.1: KMeans clustering of UMAP embedding
# Estimate optimal k
fviz_nbclust(umap.res, kmeans, method='silhouette')
# Run
kmeans.res <- kmeans(umap.res, 3)
fviz_cluster(kmeans.res,
             data = umap.res,
             nmf.basispalette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw())
umap.res$kmeans.cluster <- kmeans.res$cluster

#### Convert distances in each cluster to contact frequency
group1.bind <- plotContactFrequency(Forward_all_distances_adj.DTag[which(dist.upper.lens==100)[which(pca.res$kmeans.cluster == 1)]],
                           full_dist.DTag,
                           dist.upper.lens,
                           10)
group2.bind <- plotContactFrequency(Forward_all_distances_adj.DTag[which(dist.upper.lens==100)[which(pca.res$kmeans.cluster == 2)]],
                                    full_dist.DTag,
                                    dist.upper.lens,
                                    10)
group3.bind <- plotContactFrequency(Forward_all_distances_adj.DTag[which(dist.upper.lens==100)[which(pca.res$kmeans.cluster == 3)]],
                                    full_dist.DTag,
                                    dist.upper.lens,
                                    10)
pheatmap(group1.bind, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026","#4a0200"))(1000))
pheatmap(group2.bind, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026","#4a0200"))(1000))
pheatmap(group3.bind, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026","#4a0200"))(1000))

##### 3: Assessing the distinctness of the clusters
# Data to be assessed: reduced or full?


#### Put all heatmaps from each cluster into a dir
clusters <- 1:3
for (i in clusters){
  trace <- 1
  for (j in Forward_all_distances_adj.DTag[which(dist.upper.lens==100)[which(nmf.basis$kmeans.cluster == i)]]){
    pdf(paste0('nmf_kmeans_matrices/cluster',i,'/',str,'cluster',i,'_NMF_KMeans-3clus_trace',trace,'.pdf'))
    pheatmap(j, cluster_rows = F, cluster_cols = F, color = rev(colorRampPalette(c("#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026","#4a0200"))(1000)))
    dev.off()
    trace<-trace+1
  }
}



