library(ggplot2)
library(reshape2)
library(ggbeeswarm)

setwd("D:/Workspace/毕设/3.10/Test/Graph")
png(file = "Violin.png")

dev.off()

term<-as.matrix(res[[1]][c('Time'),])
term<-cbind(as.matrix(res[[2]][c('Time'),]),term)
term<-cbind(as.matrix(res[[3]][c('Time'),]),term)
term<-term[-1,]
term
colnames(term)<-c('ASSO','PANDA','ASSO')
#term<-t(term)
#data_m <- melt(term)
data_m <- melt(term,id.vars="Algorithm",variable.name="Term",value.name="Value")
data_m<-data_m[,-1]
colnames(data_m)<-c('Algorithm','Value')
data_m



#四舍五入
#data_m[,2]<-round(data_m[,2]) 
#data_m

#############################################################################
# 小提琴图
{
name<-c('Accuracy')

term<-as.matrix(res[[1]][name,])
term<-cbind(as.matrix(res[[2]][name,]),term)
term<-cbind(as.matrix(res[[3]][name,]),term)
term<-term[-1,]
colnames(term)<-c('MEBF','PANDA','ASSO')
data_m <- melt(term,id.vars="Algorithm",variable.name="Term",value.name=name)
data_m<-data_m[,-1]
colnames(data_m)<-c('Algorithm',name)
setwd("D:/Workspace/毕设/3.10/Test/Graph")
png(file = paste(name,'png',sep='.'))
# violin

#这里还要再改一次
p <- ggplot(data_m, aes(x=Algorithm, y=Accuracy),color=Value) +
    geom_violin(aes(fill=factor(Algorithm))) +
    geom_quasirandom() +
    theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) +
    theme(legend.position="none")+labs(x="算法",y="准确率",title="准确率（Accuracy）")+theme(plot.title=element_text(hjust=0.5))
p
dev.off()

# Box
p <- ggplot(data_m, aes(x=Algorithm, y=Accuracy),color=Algorithm) +
geom_boxplot(aes(fill=factor(Algorithm))) +
theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) +
theme(legend.position="none")
p



library(ggbeeswarm)


theme_bw() + theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), legend.key=element_blank()) +
theme(legend.position="none")








}
#############################################################################

#堆叠柱形图
{
term<-as.matrix(res[[1]][c('TP','TN','FP','FN'),1])
term<-cbind(as.matrix(res[[2]][c('TP','TN','FP','FN'),1]),term)
term<-cbind(as.matrix(res[[3]][c('TP','TN','FP','FN'),1]),term)
colnames(term)<-c('MEBF','PANDA','ASSO')
term<-t(term)
data_m <- melt(term,id.vars="Algorithm",variable.name="Term",value.name="Value")
colnames(data_m)<-c('Algorithm','Term','Value')
data_m[,3]<-round(data_m[,3])
setwd("D:/Workspace/毕设/3.10/Test/Graph")
png(file = "TP TN FP FN of Three Algorithms.png")
ggplot(data_m,aes(Algorithm,Value,fill=Term))+
geom_bar(stat="identity",position="stack")+
geom_text(aes(label=Value), position=position_stack(vjust=0.5))+
ggtitle("TP TN FP FN of Three Algorithms")+
theme_wsj()+
scale_fill_wsj("rgby", "")+
theme(axis.ticks.length=unit(0.5,'cm'))+
guides(fill=guide_legend(title=NULL))
dev.off()
}

#############################################################################

# 散点图
{
name<-c('Time')

term<-as.matrix(res[[1]][name,])
term<-cbind(as.matrix(res[[2]][name,]),term)
term<-cbind(as.matrix(res[[3]][name,]),term)
term<-term[-1,]
colnames(term)<-c('MEBF','PANDA','ASSO')

#term<-t(term)
#data_m <- melt(term)
data_m <- melt(term,id.vars="Algorithm",variable.name="Term",value.name=name)
colnames(data_m)<-c('Times','Algorithm',name)
#data_m
setwd("D:/Workspace/毕设/3.10/Test/Graph")
png(file = paste(name,'png',sep='.'))
p<-ggplot(data_m,aes(x=Times,y=Time,colour=Algorithm))+geom_point()
p
dev.off()
}

# Jaccard Distance
{
term<-as.matrix(resnoaverage[[1]][c('TP','TN','FP','FN'),1])
term<-cbind(as.matrix(resnoaverage[[2]][c('TP','TN','FP','FN'),1]),term)
term<-cbind(as.matrix(resnoaverage[[3]][c('TP','TN','FP','FN'),1]),term)
term
jaccard1<-(term[3,1]+term[4,1])/(term[4,1]+term[3,1]+term[1,1])
jaccard2<-(term[3,2]+term[4,2])/(term[4,2]+term[3,2]+term[1,2])
jaccard3<-(term[3,3]+term[4,3])/(term[4,3]+term[3,3]+term[1,3])
term<-t(as.matrix(c(jaccard1,jaccard2,jaccard3),byrow=true))
colnames(term)<-c('MEBF','PANDA','ASSO')
rownames(term)<-c('Jaccard Distance')
data_m <- melt(term)
data_m<-data_m[,-1]
colnames(data_m)<-c('Algorithm','Jaccard Distance')
data_m
}

caret
d1f<-factor(d1)
d2f<-factor(d2)
confusionMatrix(d1f,d2f)