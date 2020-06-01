# source('D:/Workspace/毕设/3.10/Test/R/BMFAlgorithmsCompare.R')
# res<-BMFAlgorithmsCompare()



library(ggplot2)
library(reshape2)
library(ggbeeswarm)
library(gcookbook)
library(ggthemes)
library(caret)
BMFAlgorithmsCompare <- function(n=300,k=5,p1=0.3,p2=0.01,DIM=10,Thres=0.6,Times=11){

	output<-list()
	
###########################################################################################
	# Run 
	omASSO <-OutputMatrix()
	omPANDA <-OutputMatrix()
	omMEBF <-OutputMatrix()
	for(i in 1:Times){
	
	# Generate Data with Noise.
	input<-AddNoise(DataGenerator(n,k,p1),p2)
	
	# Each Time
	omASSO<-cbind(ASSOPerformance(input,DIM,Thres),omASSO)
	omPANDA<-cbind(PANDAPerformance(input,DIM),omPANDA)
	omMEBF<-cbind(MEBFPerformance(input,DIM,Thres),omMEBF)
	}
	
	# The Average Performance Of Each Algorithm
	omASSO<-AveragePerformance(omASSO,Times)
	omPANDA<-AveragePerformance(omPANDA,Times)
	omMEBF<-AveragePerformance(omMEBF,Times)
	
	output[[1]]<-omASSO
	output[[2]]<-omPANDA
	output[[3]]<-omMEBF
	 
  return(output)
}

# RowCompare
# This function has the same level with BMFAlgorithmsCompare
RowCompare<-function(input,rec,name='row1'){
input<-as.matrix(input[name,],byrow=T)
rec<-as.matrix(rec[name,],byrow=T)
om <-OutputMatrixrc()
# PN
tp<-0#1-1
tn<-0#0-0
fp<-0#0-1
fn<-0#1-0
p<-MatL1(input)
n<-((nrow(input)*ncol(input))-p)
p1<-MatL1(rec)
n1<-((nrow(input)*ncol(input))-p1)

  for(x in 1:nrow(input)){
  for(y in 1:ncol(input)){
	# Cover 1
    if(input[x,y]==1){
	if(rec[x,y]==1){tp<-tp+1}
	if(rec[x,y]==0){fn<-fn+1}
	}
	if(input[x,y]==0){
	if(rec[x,y]==0){tn<-tn+1}
	if(rec[x,y]==1){fp<-fp+1}
	}
  }
 
}
res<-c(tp,tn,fp,fn)
 return (res)
}

# ColCompare
ColCompare<-function(input,rec,name='col1'){
input<-as.matrix(input[,name])
rec<-as.matrix(rec[,name])
om <-OutputMatrixrc()
# PN
tp<-0#1-1
tn<-0#0-0
fp<-0#0-1
fn<-0#1-0
p<-MatL1(input)
n<-((nrow(input)*ncol(input))-p)
p1<-MatL1(rec)
n1<-((nrow(input)*ncol(input))-p1)

  for(x in 1:nrow(input)){
  for(y in 1:ncol(input)){
	# Cover 1
    if(input[x,y]==1){
	if(rec[x,y]==1){tp<-tp+1}
	if(rec[x,y]==0){fn<-fn+1}
	}
	if(input[x,y]==0){
	if(rec[x,y]==0){tn<-tn+1}
	if(rec[x,y]==1){fp<-fp+1}
	}
  }
   
}
res<-c(tp,tn,fp,fn)
 return (res)
}


# Generate Algorithm Performance Data Matrix
OutputMatrix <- function(){
	# om <-matrix(ncol=3,nrow=6)
	om <-matrix(nrow=16)
	rownames(om) <- c('Time','Correct','Accuracy','ErrorRate','Precision','Recall','F1Score','Reconstruction Error','Density','Sensitive','Specificity','Jaccard Distance','TP','TN','FP','FN')
	
  return(om)
}


###################################################################################################################################################################
# Generate Algorithm Performance Data Of Each Algorithm
ASSOPerformance<- function(input,DIM,Thres){

# Decompose And Timing
t1=proc.time()
dec<-ASSO(input,DIM,Thres)
t2=proc.time()
t<-t2-t1
t<-t[[1]]

###########################################################################################
# Recompose the matrix.
rec<-Recompose(dec)

# Difference Matrix in -1/0/1
difm<-DifferenceMatrix(input,rec)

# Difference Matrix in 0/1. Use this matrix to calculate the L1 norm of a matrix.
difm01<-DifferenceMatrix01(input,rec)
###########################################################################################
# PN
tp<-0#1-1
tn<-0#0-0
fp<-0#0-1
fn<-0#1-0
p<-MatL1(input)
n<-((nrow(input)*ncol(input))-p)
p1<-MatL1(rec)
n1<-((nrow(input)*ncol(input))-p1)

  for(x in 1:nrow(input)){
  for(y in 1:ncol(input)){
	# Cover 1
    if(input[x,y]==1){
	if(rec[x,y]==1){tp<-tp+1}
	if(rec[x,y]==0){fn<-fn+1}
	}
	if(input[x,y]==0){
	if(rec[x,y]==0){tn<-tn+1}
	if(rec[x,y]==1){fp<-fp+1}
	}
  }
}
###########################################################################################
# 'Time','Correct','Accuracy','ErrorRate','Precision','Recall','F1 Score','ReconstructionError','Density','Sensitive','Specificity','TP','TN','FP','FN')

# Correct
correct<-tp+tn

# Accuracy
accuracy<-(tp+tn)/(tp+tn+fp+fn)

# ErrorRate
errorrate<-(fp+fn)/(tp+tn+fp+fn)

# Precision
precision<-tp/(tp+fp)

# Recall
recall<-tp/p

#F1 Score
F1<-(2*precision*recall)/(precision+recall)

# Reconstruction Error.
re<-ReconstructionError(input,difm01)

# Density.
d<-Density(dec)

# 1 Coverage Rate
#

#JaccardDistance
jd<-JaccardDistance(fp,fn,tp)

# Sensitive
sensitive<-tp/p

# Specificity
specificity<-tn/n

###########################################################################################
# Generate the image of Difference Matrix. Red means wrong. White means correct.
#library(raster)
#library(latticeExtra)
#png(filename="D:/Workspace/毕设/3.10/Test/output/ASSO.png")
#spplot(raster(difm01), colorkey=FALSE, col.regions=c('white', 'red')) + layer(panel.grid(h=nrow(difm)-1, v=ncol(difm01)-1, col=1))
#dev.off()

###########################################################################################
# Generate Output
# '1Time','2Correct','3Accuracy','4ErrorRate','5Precision','6Recall','7F1 Score','8ReconstructionError','9Density','10Sensitive','11Specificity',16overlaperror','12TP','13TN','14FP','15FN'
	count<-c(t,correct,accuracy,errorrate,precision,recall,F1,re,d,sensitive,specificity,jd,tp,tn,fp,fn)
    return(count)
}

###########################################################################################
# Generate Algorithm Performance Data Of PANDA Algorithm
PANDAPerformance<- function(input,DIM){

# Decompose And Timing
t1=proc.time()
dec<-PANDA(input,DIM)
t2=proc.time()
t<-t2-t1
t<-t[[1]]

###########################################################################################
# Recompose the matrix.
rec<-Recompose(dec)

# Difference Matrix in -1/0/1
difm<-DifferenceMatrix(input,rec)

# Difference Matrix in 0/1. Use this matrix to calculate the L1 norm of a matrix.
difm01<-DifferenceMatrix01(input,rec)
###########################################################################################
# PN
tp<-0#1-1
tn<-0#0-0
fp<-0#0-1
fn<-0#1-0
p<-MatL1(input)
n<-((nrow(input)*ncol(input))-p)
p1<-MatL1(rec)
n1<-((nrow(input)*ncol(input))-p1)

  for(x in 1:nrow(input)){
  for(y in 1:ncol(input)){
	# Cover 1
    if(input[x,y]==1){
	if(rec[x,y]==1){tp<-tp+1}
	if(rec[x,y]==0){fn<-fn+1}
	}
	if(input[x,y]==0){
	if(rec[x,y]==0){tn<-tn+1}
	if(rec[x,y]==1){fp<-fp+1}
	}
  }
}
###########################################################################################
# 'Time','Correct','Accuracy','ErrorRate','Precision','Recall','F1 Score','ReconstructionError','Density','Sensitive','Specificity','TP','TN','FP','FN')

# Correct
correct<-tp+tn

# Accuracy
accuracy<-(tp+tn)/(tp+tn+fp+fn)

# ErrorRate
errorrate<-(fp+fn)/(tp+tn+fp+fn)

# Precision
precision<-tp/(tp+fp)

# Recall
recall<-tp/p

#F1 Score
F1<-(2*precision*recall)/(precision+recall)

# Reconstruction Error.
re<-ReconstructionError(input,difm01)

# Density.
d<-Density(dec)

# 1 Coverage Rate



#JaccardDistance
jd<-JaccardDistance(fp,fn,tp)

# Sensitive
sensitive<-tp/p

# Specificity
specificity<-tn/n

###########################################################################################
# Generate the image of Difference Matrix. Red means wrong. White means correct.
#library(raster)
#library(latticeExtra)
#png(filename="D:/Workspace/毕设/3.10/Test/output/ASSO.png")
#spplot(raster(difm01), colorkey=FALSE, col.regions=c('white', 'red')) + layer(panel.grid(h=nrow(difm)-1, v=ncol(difm01)-1, col=1))
#dev.off()

###########################################################################################
# Generate Output
# '1Time','2Correct','3Accuracy','4ErrorRate','5Precision','6Recall','7F1 Score','8ReconstructionError','9Density','10Sensitive','11Specificity','12TP','13TN','14FP','15FN'
	count<-c(t,correct,accuracy,errorrate,precision,recall,F1,re,d,sensitive,specificity,jd,tp,tn,fp,fn)
    return(count)
}

###########################################################################################
# Generate Algorithm Performance Data Of MEBF Algorithm
MEBFPerformance<- function(input,DIM,Thres){

# Decompose And Timing
t1=proc.time()
dec<-MEBF(input,DIM,Thres)
t2=proc.time()
t<-t2-t1
t<-t[[1]]

###########################################################################################
# Recompose the matrix.
rec<-Recompose(dec)

# Difference Matrix in -1/0/1
difm<-DifferenceMatrix(input,rec)

# Difference Matrix in 0/1. Use this matrix to calculate the L1 norm of a matrix.
difm01<-DifferenceMatrix01(input,rec)
###########################################################################################
# PN
tp<-0#1-1
tn<-0#0-0
fp<-0#0-1
fn<-0#1-0
p<-MatL1(input)
n<-((nrow(input)*ncol(input))-p)
p1<-MatL1(rec)
n1<-((nrow(input)*ncol(input))-p1)

  for(x in 1:nrow(input)){
  for(y in 1:ncol(input)){
	# Cover 1
    if(input[x,y]==1){
	if(rec[x,y]==1){tp<-tp+1}
	if(rec[x,y]==0){fn<-fn+1}
	}
	if(input[x,y]==0){
	if(rec[x,y]==0){tn<-tn+1}
	if(rec[x,y]==1){fp<-fp+1}
	}
  }
}
###########################################################################################
# 'Time','Correct','Accuracy','ErrorRate','Precision','Recall','F1 Score','ReconstructionError','Density','Sensitive','Specificity','TP','TN','FP','FN')

# Correct
correct<-tp+tn

# Accuracy
accuracy<-(tp+tn)/(tp+tn+fp+fn)

# ErrorRate
errorrate<-(fp+fn)/(tp+tn+fp+fn)

# Precision
precision<-tp/(tp+fp)

# Recall
recall<-tp/p

#F1 Score
F1<-(2*precision*recall)/(precision+recall)

# Reconstruction Error.
re<-ReconstructionError(input,difm01)

# Density.
d<-Density(dec)

# 1 Coverage Rate



#JaccardDistance
jd<-JaccardDistance(fp,fn,tp)

# Sensitive
sensitive<-tp/p

# Specificity
specificity<-tn/n

###########################################################################################
# Generate the image of Difference Matrix. Red means wrong. White means correct.
#library(raster)
#library(latticeExtra)
#png(filename="D:/Workspace/毕设/3.10/Test/output/ASSO.png")
#spplot(raster(difm01), colorkey=FALSE, col.regions=c('white', 'red')) + layer(panel.grid(h=nrow(difm)-1, v=ncol(difm01)-1, col=1))
#dev.off()

###########################################################################################
# Generate Output
# '1Time','2Correct','3Accuracy','4ErrorRate','5Precision','6Recall','7F1 Score','8ReconstructionError','9Density','10Sensitive','11Specificity','12TP','13TN','14FP','15FN'
	count<-c(t,correct,accuracy,errorrate,precision,recall,F1,re,d,sensitive,specificity,jd,tp,tn,fp,fn)
    return(count)
}
####################################################################################################################################################################

# Tools

# Data Generator 
# (X^n*k)*(X^k*m) Bernoulli(p)
DataGenerator<-function(n,k,p){
list<-list()
u<-matrix(rbinom((n*k),1,p),nrow=n,ncol=k)
v<-matrix(rbinom((k*n),1,p),nrow=k,ncol=n)
list[[1]]<-u
list[[2]]<-v
test<-Recompose(list)
return (test)
}

# Add Noise
# Flip test Bernoulli(p)
AddNoise<-function(test,p){
x<-test
noise<-matrix(rbinom((nrow(test)*ncol(test)),1,p),nrow=nrow(test),ncol=ncol(test))
for(i in 1:nrow(test)){
	for(j in 1:ncol(test)){
	if(noise[i,j]==1){
	if(test[i,j]==0){x[i,j]<-1}
	if(test[i,j]==1){x[i,j]<-0}
	}
	}
}
return (x)
}

# the product of two boolean matrics
Recompose<- function(list){
	library(Matrix)
	rec<-as.matrix(list[[1]]%&%list[[2]])
	rec<-rec*1
    return(rec)
}

# -1/0/1
DifferenceMatrix <- function(mat1, mat2){
  result<-matrix(nrow=nrow(mat1),ncol=ncol(mat1))
	if(ncol(mat1)!=ncol(mat2)){return (-1)}
	if(nrow(mat1)!=nrow(mat2)){return (-1)}
for(x in 1:nrow(mat1)){
  for(y in 1:ncol(mat1)){
	# correct
    if(mat1[x,y]==mat2[x,y]){result[x,y]<-0}
	# 0 to 1
    if(mat1[x,y]<mat2[x,y]){result[x,y]<--1}
	# 1 to 0
    if(mat1[x,y]>mat2[x,y]){result[x,y]<-1}
	
  }
}
  return(result)
}

# 0/1
DifferenceMatrix01 <- function(mat1, mat2){
  result<-matrix(nrow=nrow(mat1),ncol=ncol(mat1))
for(x in 1:nrow(mat1)){
  for(y in 1:ncol(mat1)){
	# correct=0
    if(mat1[x,y]==mat2[x,y]){result[x,y]<-0}
	# wrong=1
    else{result[x,y]<-1}
  }
}
  return(result)
}

# L1 norm of a matrix
MatL1<-function(mat){
    rs<-rowSums(mat)
    L1<-sum(rs)
    return (L1)
}

# Reconstruction Error
ReconstructionError <- function(input, difm01){
re<-(MatL1(difm01)/MatL1(input))
  return(re)
}

# Density
Density <- function(list){
  d<- (MatL1(list[[1]])+MatL1(list[[2]])) / (nrow(list[[1]])+ncol(list[[2]])*nrow(list[[1]]))
  return(d)
}

# Average Performance 
# Rename Columns
AveragePerformance<- function(om,Times){
	colnames(om) <- colnames(om, do.NULL = FALSE, prefix = "Time_")
	om<-om[,-(Times+1)]
	ave<-as.matrix(rowMeans(om))
	colnames(ave) <- c('Average')
	om<-cbind(ave,om)
	return (om)
}

#Overlap Error
OverlapError<- function(list){
	oe<-0
	#row sum
	rs<-rowSums(list[[1]])
	#col sum
	cs<-colSums(list[[2]])
	
	#Overlap Error
	for(i in 1:nrow(list[[1]])){
	if(rs[i]>1){oe<-oe+(rs[i]-1)}
	}
	for(j in 1:ncol(list[[2]])){
	if(cs[j]>1){oe<-oe+(cs[i]-1)}
	}
  return(oe)
}

#JaccardDistance
JaccardDistance<-function(fp,fn,tp){
jd<-((fp+fn)/(fp+fn+tp))
return (jd)
}

###########################################################################################################################################################
# ASSO Algorithm
ASSO<-function(MAT,DIM=10,Thres){
  
  if(min(ncol(MAT),nrow(MAT))<DIM){
    DIM<-min(ncol(MAT),nrow(MAT))-1
  }
  
  
  
  X<-MAT %*% t(MAT)
  VEC<-diag(X)
  X<-sweep(X,2,VEC,FUN = "/")
  X<-apply(X>Thres, 2, as.numeric)
  X[is.na(X)]<-0
  
  A<-MAT
  B<-NULL
  C<-NULL
  D<-A*0
  Thres_use<-0.05*sum(A)
  COL<-0
  NOW<-0
  for (i in 1:ncol(X)) {
    B_TEMP<-X[,i]
    C_TEMP<-rep(1,ncol(A))
    A_TEMP<-B_TEMP%o%C_TEMP
    TEMP<-colSums(A_TEMP*A*2-A_TEMP)
    C_TEMP<-as.numeric(TEMP>0)
    if(NOW<sum(TEMP[TEMP>0])){
      B<-B_TEMP
      C<-C_TEMP
      NOW<-sum(TEMP[TEMP>0])
      COL<-i
    }
  }
  X<-X[,-COL]
  A_stop<-sum(apply((A-B%o%C)>0,2,as.numeric))
  
  
  for (m in 2:DIM) {
    
    if(A_stop>Thres_use){
      B_use<-NULL
      C_use<-NULL
      COL<-0
      for (i in 1:ncol(X)) {
        B_TEMP<-cbind(B,X[,i])
        C_TEMP<-rbind(C,rep(1,ncol(A)))
        A_TEMP<-apply((B_TEMP%*%C_TEMP)>0,2,as.numeric)
        TEMP<-colSums(A_TEMP*A*2-A_TEMP)
        if(NOW<sum(TEMP[TEMP>0])){
          NOW<-sum(TEMP[TEMP>0])
          B_use<-X[,i]
          C_use<-as.numeric(TEMP>0)
          COL<-i
        }
      }
      
      
      if(COL==0){
        break
      }
      B<-cbind(B,B_use)
      C<-rbind(C,C_use)
      X<-X[,-COL]
      A_stop<-sum(apply((A-B%*%C)>0,2,as.numeric))
    }
  }
  
  result<-list()
  result[[1]]<-B
  result[[2]]<-C
  return(result)
}

# MEBF Algorithm
MEBF<-function(MAT,DIM=10,Thres){
  
  if(min(ncol(MAT),nrow(MAT))<DIM){
    DIM<-min(ncol(MAT),nrow(MAT))-1
  }
  
  m1<-MAT
  SUM<-sum(MAT)
  MAT_B<-NULL
  MAT_C<-NULL
  for (i in 1:DIM) {
    
    if(sum(m1)<=0.05*SUM){
      result<-list()
      result[[1]]<-MAT_B
      result[[2]]<-MAT_C
      return(result)
    }
    
    C1<-0
    B1<-rep(0,nrow(m1))
    B1_use<-B1
    B2<-rep(0,ncol(m1))
    B2_use<-B2
    COL<-colSums(m1)
    ROW<-rowSums(m1)
    
    # START with column
    if(!is.na(median(COL[COL>1]))&median(COL[COL>1])>1){
      TEMP<-which(colSums(m1)==min(COL[which(COL>=median(COL[COL>1]))]))
      if(length(TEMP)==1){
        B1<-m1[,TEMP]
        B2[which(colSums(m1[which(B1==1),])>=min(Thres*sum(B1)+1,sum(B1)))]<-1
        C2<-(sum(B1)-1)*(sum(B2)-1)-sum(m1[which(B1==1),which(B2==1)]==0)
        if(C2>C1){
          C1<-C2
          B1_use<-B1
          B2_use<-B2
        }
        B1<-rep(0,nrow(m1))
        B2<-rep(0,ncol(m1))
      }else{
        for(j in 1:length(TEMP)){
          B1<-m1[,TEMP[j]]
          B2[which(colSums(m1[which(B1==1),])>=min(Thres*sum(B1)+1,sum(B1)))]<-1
          C2<-(sum(B1)-1)*(sum(B2)-1)-sum(m1[which(B1==1),which(B2==1)]==0)
          if(C2>C1){
            C1<-C2
            B1_use<-B1
            B2_use<-B2
          }
          B1<-rep(0,nrow(m1))
          B2<-rep(0,ncol(m1))
        }
      }
    }
    
    # START with ROW
    if(!is.na(median(ROW[ROW>1]))&median(ROW[ROW>1])>1){
      TEMP<-which(rowSums(m1)==min(ROW[ROW>=median(ROW[ROW>1])]))
      if(length(TEMP)==1){
        B2<-m1[TEMP,]
        B1[which(rowSums(m1[,which(B2==1)])>=min(Thres*sum(B2)+1,sum(B2)))]<-1
        C2<-(sum(B1)-1)*(sum(B2)-1)-sum(m1[which(B1==1),which(B2==1)]==0)
        if(C2>C1){
          C1<-C2
          B1_use<-B1
          B2_use<-B2
        }
        B1<-rep(0,nrow(m1))
        B2<-rep(0,ncol(m1))
      }else{
        for(j in 1:length(TEMP)){
          B2<-m1[TEMP[j],]
          B1[which(rowSums(m1[,which(B2==1)])>=min(Thres*sum(B2)+1,sum(B2)))]<-1
          C2<-(sum(B1)-1)*(sum(B2)-1)-sum(m1[which(B1==1),which(B2==1)]==0)
          if(C2>C1){
            C1<-C2
            B1_use<-B1
            B2_use<-B2
          }
          B1<-rep(0,nrow(m1))
          B2<-rep(0,ncol(m1))
        }
      }
    }
    
    if(C1==0){
      ROW<-order(rowSums(m1),decreasing = T)
      COL<-order(colSums(m1),decreasing = T)
      # start from ROW
      B1_1<-rep(0,nrow(m1))
      if(length(which(rowSums(m1[,COL[1:2]])==2))>1){
        B1_1[which(rowSums(m1[,COL[1:2]])==2)]<-1
        B1_2<-rep(0,ncol(m1))
        B1_2[colSums(m1[which(B1_1==1),])>=min((Thres*sum(m1[which(B1_1==1),COL[1]])+1),sum(m1[which(B1_1==1),COL[1]]))]<-1
        C1<-(sum(B1_1)-1)*(sum(B1_2)-1)-sum(m1[which(B1_1==1),which(B1_2==1)]==0)
      }else{
        C1<-(-Inf)
      }
      
      # start from COL
      B2_2<-rep(0,ncol(m1))
      if(length(which(colSums(m1[ROW[1:2],])==2))>1){
        B2_2[which(colSums(m1[ROW[1:2],])==2)]<-1
        B2_1<-rep(0,nrow(m1))
        B2_1[rowSums(m1[,which(B2_2==1)])>=min(Thres*sum(m1[ROW[1],which(B2_2==1)])+1,sum(m1[ROW[1],which(B2_2==1)]))]<-1
        C2<-(sum(B2_1)-1)*(sum(B2_2)-1)-sum(m1[which(B2_1==1),which(B2_2==1)]==0)
      }else{
        C2<-(-Inf)
      }
      
      if((C1==(-Inf))&(C2==(-Inf))){
        break
      }else{
        if(C1>C2){
          B1_use<-B1_1
          B2_use<-B1_2
          m1[which(B1_1==1),which(B1_2==1)]<-0
        }else{
          B1_use<-B2_1
          B2_use<-B2_2
          m1[which(B2_1==1),which(B2_2==1)]<-0
        }
      }
    }
    MAT_B<-cbind(MAT_B,B1_use)
    MAT_C<-rbind(MAT_C,B2_use)
    m1[which(MAT_B[,i]==1),which(MAT_C[i,]==1)]<-0
  }
  result<-list()
  result[[1]]<-MAT_B
  result[[2]]<-MAT_C
  return(result)
}

# PANDA Algorithm
GAMMA<-function(MAT_A,MAT_AC,Core){
  TEMP<-Core[[1]]%o%Core[[2]]
  TEMP<-TEMP+MAT_AC
  TEMP<-apply(TEMP>0, 2, as.numeric)
  TEMP<-abs(MAT_A-TEMP)
  return(sum(TEMP))
}

Find_Core<-function(MAT_A,MAT_AR,MAT_AC){
  Core<-list()
  E<-NULL
  S<-order(rowSums(MAT_AR),decreasing = T)
  Cl<-rep(0,nrow(MAT_AR))
  Cr<-rep(0,ncol(MAT_AR))
  Cl[S[1]]<-1
  Cr<-MAT_AR[S[1],]
  Core[[1]]<-Cl
  Core[[2]]<-Cr
  G<-GAMMA(MAT_A,MAT_AC,Core)+sum(Core[[1]])+sum(Core[[2]])
  for (i in 2:length(S)) {
    Core1<-list()
    Cl1<-Cl
    Cr1<-Cr
    Cl1[S[i]]<-1
    Cr1[which(MAT_AR[S[i],]==0)]<-0
    Core1[[1]]<-Cl1
    Core1[[2]]<-Cr1
    G1<-GAMMA(MAT_A,MAT_AC,Core1)+sum(Core1[[1]])+sum(Core1[[2]])
    if(G1<G){
      Core<-Core1
      G<-G1
    }else{
      E<-append(E,S[i])
    }
  }
  
  Core[[3]]<-E
  names(Core)<-c("left","right","Extension")
  return(Core)
}

Extend_Core<-function(Core,MAT_A,MAT_AC){
  
  Cl<-Core[[1]]
  Cr<-Core[[2]]
  
  G<-GAMMA(MAT_A,MAT_AC,Core)+sum(Core[[1]])+sum(Core[[2]])
  for (i in Core[[3]]) {
    Core1<-Core[1:2]
    Core1[[1]]<-Cl
    Core1[[1]][i]<-1
    G1<-GAMMA(MAT_A,MAT_AC,Core1)+sum(Core1[[1]])+sum(Core1[[2]])
    if(G1<G){
      Core[1:2]<-Core1
      G<-G1
    }
  }
  for (i in which(Core[[2]]==0)) {
    Core1<-Core[1:2]
    Core1[[2]]<-Cr
    Core1[[2]][i]<-1
    G1<-GAMMA(MAT_A,MAT_AC,Core1)+sum(Core1[[1]])+sum(Core1[[2]])
    if(G1<G){
      Core[1:2]<-Core1
      G<-G1
    }
  }
  return(Core)
} 

PANDA<-function(MAT,DIM=10){
  
  if(min(ncol(MAT),nrow(MAT))<DIM){
    DIM<-min(ncol(MAT),nrow(MAT))-1
  }
  Thres_use<-0.05*sum(MAT)
  MAT_A<-MAT
  MAT_B<-NULL
  MAT_C<-NULL
  MAT_AR<-MAT_A
  MAT_AC<-0*MAT
  G<-sum(MAT_A)
  G1<-0
  for (i in 1:DIM) {
    if(sum(MAT_AR)<=Thres_use){
      result<-list()
      result[[1]]<-MAT_B
      result[[2]]<-MAT_C
      return(result)
      break()
    }else{
      Core<-Find_Core(MAT_A,MAT_AR,MAT_AC)
      Core<-Extend_Core(Core,MAT_A,MAT_AC)
      MAT_B<-cbind(MAT_B,Core[[1]])
      MAT_C<-rbind(MAT_C,Core[[2]])
      
      TEMP<-Core[[1]] %o% Core[[2]]
      MAT_AR<-MAT_AR-TEMP
      MAT_AR<-apply(MAT_AR>0,2,as.numeric)
      MAT_AC<-MAT_AC+TEMP
      MAT_AC<-apply(MAT_AC>0,2,as.numeric)
      
      G1<-sum(abs(MAT_A-MAT_AC))+sum(MAT_B)+sum(MAT_C)
      if(G1<G){
        G<-G1
      }else{
        result<-list()
        result[[1]]<-MAT_B
        result[[2]]<-MAT_C
        return(result)
        break()
      }
    }
  }
  result<-list()
  result[[1]]<-MAT_B
  result[[2]]<-MAT_C
  return(result)
}
