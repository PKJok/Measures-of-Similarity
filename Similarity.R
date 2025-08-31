getwd()
setwd("C:/Users/ASUS/OneDrive/Desktop/R for Quantitative Genetics")

# load the libraries
library(tidyverse)
library(factoextra)
library(ggfortify)
library(vegan)

# simulating genotypic data
set.seed(123)
n_genotypes<- 10
n_markers<- 100 # 100 SNP markers

# Simulate marker data: 0, 1, 2
# if one of the alleles is the same, score the locus = 1;

# 0 = homozygous alternative alleles, 1= , and 2 to represent homozygous reference, heterozygous, and .

genotype_data <- matrix(sample(0:2, n_genotypes * n_markers, replace = TRUE),
                        nrow = n_genotypes,
                        ncol = n_markers)

rownames(genotype_data) <- paste0("Genotype_", 1:n_genotypes)
colnames(genotype_data) <- paste0("SNP_", 1:n_markers)


# lets do pca analysis and cluster the PC1 and 2 for initial exploration
pca<- prcomp(genotype_data, scale. = TRUE)

pca_per<- round(pca$sdev^2/ sum(pca$sdev^2)*100 ,2)

pca$x%>%
  as.data.frame()%>%
  select(PC1, PC2)%>%
  rownames_to_column(var = "Genotype")%>%
  ggplot(., aes(PC1,PC2, colour = Genotype))+
  geom_point(size=3)+
  theme_bw()+
  ylab(paste("PC1 (",pca_per[1],"%)"))+
  xlab(paste("PC2 (", pca_per[2],"%)"))

pca_df<- data.frame(pc1= pca$x[,1], pc2= pca$x[,2])

# distance calculation
dist_matrix<- dist(pca_df, method = "euclidean")
clust<- hclust(dist_matrix, method = "ward.D2")
plot(clust, hang = -0.2, main = "PCA Cluster")


# similarity index
# for genotype x and y
similarity_index<- function(x,y){
  sum(x * y)/ (2* length(x))
}


# computing pairwise similarity matrix
similarity_matrix<- matrix(0, nrow = n_genotypes, ncol = n_genotypes)

for (i in 1:n_genotypes){
  for (j in 1: n_genotypes){
    similarity_matrix[i, j] <- similarity_index(genotype_data[i, ], genotype_data[j, ])
  }
}

# setting the rownames and colnames of the similarity_matrix

rownames(similarity_matrix)<- rownames(genotype_data)
colnames(similarity_matrix)<- rownames(genotype_data)

# convert to distance matrix

distance_matrix= 1- similarity_matrix

heatmap(similarity_matrix)
heatmap(similarity_matrix, main = "Similarity Matrix", symm = TRUE)


# Genetic Distance Measures
# Euclidean distance
dist_euclidean <- dist(genotype_data, method = "euclidean")
clust1<- hclust(dist_euclidean, method = "ward.D2")
plot(clust1, main = "Euclidean_Cluster")

# cut tree into 3 clusters
clusters <- cutree(clust1, k = 3)
clusters


# manhattan distance
dist_manhattan<- dist(genotype_data, method = "manhattan")
clust2<- hclust(dist_manhattan, method="ward.D2")
plot(clust2, hang=-0.5, main = "Manhattan_Cluster")


# Gower distance, for SNP data
dist_gower <- vegdist(genotype_data, method = "gower")
clust3<- hclust(dist_gower, method="ward.D2")
plot(clust3, hang=-0.5, main = "Gower_Cluster")


# Nei's distance
#library(ape)
#dist_nei <- dist.dna(as.DNAbin(genotype_data), model = "N")


# Kmeans clustering
set.seed(123)
kmeans_result <- kmeans(genotype_data, centers = 3)
kmeans_result$cluster

# Add cluster assignment to PCA plot
plot(pca$x[, 1], pca$x[, 2],
     col = kmeans_result$cluster, pch = 16,
     xlab = "PC1", ylab = "PC2", main = "K-means Clustering (K=3)", 
     ylim = c(-10,7), xlim = c(-8,8))
text(pca$x[, 1], pca$x[, 2], labels = rownames(genotype_data), pos = 2)


# PCA with clusters
fviz_pca_ind(pca, geom = "point", col.ind = as.factor(kmeans_result$cluster)) +
  ggtitle("PCA with K-means Clusters")

# Dendrogram with clusters
fviz_dend(clust1, k = 3, cex = 0.6, k_colors = c("royalblue", "red3", "green3"))
