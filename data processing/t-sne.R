install.packages('Mercator')
library(Mercator)
kk<-10
#输出图片的目录
setwd("D:/Workspace/毕设/3.10/Paper/Graph/MEBF01")
#########################################################################
load("origin/01/GSE84465_scRNA.Rdata")
my.binmat <- BinaryMatrix(rec)
my.binmat <- t(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)

jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="GSE84465_scRNA-hclust-Histogram of Distances")
###################################
png(file = "GSE84465_scRNA_hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()
#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "GSE84465_scRNA_tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;GSE84465_scRNA-tsne")
dev.off()


#########################################################################
load("origin/01/MCA_Lung3.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="MCA_Lung3-hclust-Histogram of Distances")
###################################
png(file = "origin_01_MCA_Lung3-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "origin_01_MCA_Lung3-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;MCA_Lung3-tsne")
dev.off()
#########################################################################
load("origin/01/MCA_Pancreas.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="MCA_Pancreas-hclust-Histogram of Distances")
###################################
png(file = "origin_01_MCA_Pancreas-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "origin_01_MCA_Pancreas-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;MCA_Pancreas-tsne")
dev.off()
#########################################################################
load("origin/01/GSE81861_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="GSE81861_scRNA-hclust-Histogram of Distances")
###################################
png(file = "origin_01_GSE81861_scRNA-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "origin_01_GSE81861_scRNA-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;GSE81861_scRNA-tsne")
dev.off()
#########################################################################
load("origin/01/GSE67835_data_c.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="GSE67835_data_c-hclust-Histogram of Distances")
###################################
png(file = "origin_01_GSE67835_data_c-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "origin_01_GSE67835_data_c-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;GSE67835_data_c-tsne")
dev.off()
#########################################################################
load("origin/01/GSE75688_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="GSE75688_scRNA-hclust-Histogram of Distances")
###################################
png(file = "origin_01_GSE75688_scRNA-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "GSE75688_scRNA-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;GSE75688_scRNA-tsne")
dev.off()






#########################################################################
load("origin/cutoff/MCA_Brain1.Rdata")
my.binmat <- BinaryMatrix(result)
my.binmat <- t(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)

jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="MCA_Brain1-hclust-Histogram of Distances")
###################################
png(file = "origin_cutoff_MCA_Brain1_hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()
#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "origin_cutoff_MCA_Brain1_tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;MCA_Brain1-tsne")
dev.off()


#########################################################################
load("origin/cutoff/MCA_Lung3.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="MCA_Lung3-hclust-Histogram of Distances")
###################################
png(file = "origin_cutoff_MCA_Lung3-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "origin_cutoff_MCA_Lung3-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;MCA_Lung3-tsne")
dev.off()
#########################################################################
load("origin/cutoff/MCA_Pancreas.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="MCA_Pancreas-hclust-Histogram of Distances")
###################################
png(file = "origin_cutoff_MCA_Pancreas-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "origin_cutoff_MCA_Pancreas-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;MCA_Pancreas-tsne")
dev.off()
#########################################################################
load("origin/cutoff/GSE81861_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="GSE81861_scRNA-hclust-Histogram of Distances")
###################################
png(file = "origin_cutoff_GSE81861_scRNA-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "origin_cutoff_GSE81861_scRNA-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;GSE81861_scRNA-tsne")
dev.off()
#########################################################################
load("origin/cutoff/GSE67835_data_c.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="GSE67835_data_c-hclust-Histogram of Distances")
###################################
png(file = "origin_cutoff_GSE67835_data_c-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "origin_cutoff_GSE67835_data_c-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;GSE67835_data_c-tsne")
dev.off()
#########################################################################
load("origin/cutoff/GSE75688_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="GSE75688_scRNA-hclust-Histogram of Distances")
###################################
png(file = "origin_cutoff_GSE75688_scRNA-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "GSE75688_scRNA-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;GSE75688_scRNA-tsne")
dev.off()









#######################################################################################################################################

#########################################################################
load("MEBF/01/MCA_Brain1.Rdata")
my.binmat <- BinaryMatrix(result)
my.binmat <- t(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)

jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="MCA_Brain1-hclust-Histogram of Distances")
###################################
png(file = "MEBF_01_MCA_Brain1_hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()
#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "MEBF_01_MCA_Brain1_tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;MCA_Brain1-tsne")
dev.off()


#########################################################################
load("MEBF/01/MCA_Lung3.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="MCA_Lung3-hclust-Histogram of Distances")
###################################
png(file = "MEBF_01_MCA_Lung3-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "MEBF_01_MCA_Lung3-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;MCA_Lung3-tsne")
dev.off()
#########################################################################
load("MEBF/01/MCA_Pancreas.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="MCA_Pancreas-hclust-Histogram of Distances")
###################################
png(file = "MEBF_01_MCA_Pancreas-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "MEBF_01_MCA_Pancreas-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;MCA_Pancreas-tsne")
dev.off()
#########################################################################
load("MEBF/01/GSE81861_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="GSE81861_scRNA-hclust-Histogram of Distances")
###################################
png(file = "MEBF_01_GSE81861_scRNA-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "MEBF_01_GSE81861_scRNA-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;GSE81861_scRNA-tsne")
dev.off()
#########################################################################
load("MEBF/01/GSE67835_data_c.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="GSE67835_data_c-hclust-Histogram of Distances")
###################################
png(file = "MEBF_01_GSE67835_data_c-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "MEBF_01_GSE67835_data_c-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;GSE67835_data_c-tsne")
dev.off()
#########################################################################
load("MEBF/01/GSE75688_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="GSE75688_scRNA-hclust-Histogram of Distances")
###################################
png(file = "MEBF_01_GSE75688_scRNA-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "GSE75688_scRNA-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;GSE75688_scRNA-tsne")
dev.off()






#########################################################################
load("origin/cutoff/MCA_Brain1.Rdata")
my.binmat <- BinaryMatrix(result)
my.binmat <- t(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)

jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="MCA_Brain1-hclust-Histogram of Distances")
###################################
png(file = "MEBF_cutoff_MCA_Brain1_hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()
#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "MEBF_cutoff_MCA_Brain1_tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;MCA_Brain1-tsne")
dev.off()


#########################################################################
load("origin/cutoff/MCA_Lung3.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="MCA_Lung3-hclust-Histogram of Distances")
###################################
png(file = "MEBF_cutoff_MCA_Lung3-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "MEBF_cutoff_MCA_Lung3-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;MCA_Lung3-tsne")
dev.off()
#########################################################################
load("origin/cutoff/MCA_Pancreas.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="MCA_Pancreas-hclust-Histogram of Distances")
###################################
png(file = "MEBF_cutoff_MCA_Pancreas-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "MEBF_cutoff_MCA_Pancreas-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;MCA_Pancreas-tsne")
dev.off()
#########################################################################
load("origin/cutoff/GSE81861_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="GSE81861_scRNA-hclust-Histogram of Distances")
###################################
png(file = "MEBF_cutoff_GSE81861_scRNA-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "MEBF_cutoff_GSE81861_scRNA-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;GSE81861_scRNA-tsne")
dev.off()
#########################################################################
load("origin/cutoff/GSE67835_data_c.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="GSE67835_data_c-hclust-Histogram of Distances")
###################################
png(file = "MEBF_cutoff_GSE67835_data_c-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "MEBF_cutoff_GSE67835_data_c-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;GSE67835_data_c-tsne")
dev.off()
#########################################################################
load("origin/cutoff/GSE75688_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)


jacc.Vis <- Mercator(my.binmat, "jaccard", "hclust", K=kk)
#hclust
hist(jacc.Vis, 
     xlab="Jaccard Distance", main="GSE75688_scRNA-hclust-Histogram of Distances")
###################################
png(file = "MEBF_cutoff_GSE75688_scRNA-hclust.png")
plot(jacc.Vis, view = "hclust")
dev.off()

#t-sne
par(pty="s")
jacc.Vis <- addVisualization(jacc.Vis, "tsne", 
            perplexity=5, 
            xlab="T1", ylab="T2")
names(jacc.Vis@view)
###################################
png(file = "GSE75688_scRNA-tsne.png")
plot(jacc.Vis, view = "tsne", main="t-SNE; Jaccard Distance;GSE75688_scRNA-tsne")
dev.off()