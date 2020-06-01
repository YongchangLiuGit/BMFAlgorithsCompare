
library(Mercator)

#数据目录
setwd("D:/Workspace/毕设/3.10/Paper/Data/MEBF/01")

#########################################################################
load("MCA_Lung3.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("MCA_Pancreas.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("MCA_Brain1.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("MCA_Stomach.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
#########################################################################


load("GSE81861_scRNA.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("GSE67835_data_c.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("GSE75688_scRNA.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("GSE84465_scRNA.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
#########################################################################
#########################################################################

setwd("D:/Workspace/毕设/3.10/Paper/Data/MEBF/cutoff")


#########################################################################
load("MCA_Lung3.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("MCA_Pancreas.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("MCA_Brain1.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("MCA_Stomach.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
#########################################################################


load("GSE81861_scRNA.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("GSE67835_data_c.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("GSE75688_scRNA.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("GSE84465_scRNA.Rdata")
my.binmat <- BinaryMatrix(rec)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
setwd("D:/Workspace/毕设/3.10/Paper/Data/origin/01")


#########################################################################
load("MCA_Lung3.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("MCA_Pancreas.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("MCA_Brain1.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("MCA_Stomach.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
#########################################################################


load("GSE81861_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("GSE67835_data_c.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("GSE75688_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("GSE84465_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
#########################################################################
#########################################################################

setwd("D:/Workspace/毕设/3.10/Paper/Data/origin/cutoff")


#########################################################################
load("MCA_Lung3.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("MCA_Pancreas.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("MCA_Brain1.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("MCA_Stomach.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
#########################################################################


load("GSE81861_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("GSE67835_data_c.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("GSE75688_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
load("GSE84465_scRNA.Rdata")
my.binmat <- BinaryMatrix(result)
summary(my.binmat)
my.binmat <- t(my.binmat)
summary(my.binmat)
my.binmat <- removeDuplicateFeatures(my.binmat)
summary(my.binmat)

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################



