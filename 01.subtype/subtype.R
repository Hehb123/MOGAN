library(ConsensusClusterPlus)
library(ggplot2)
library(GSVA)
library(GSEABase)
library(scatterplot3d)
library(tidyverse)
library(Rtsne)
library(ggsci)
load('./exp.rda')
exp =as.data.frame(t(exp))
dc = as.matrix(exp)

#################test
gset <- read.csv("CellReports.CSV",header = F)
gset =gset[,-2]
gset=t(gset)
colnames(gset)=gset[1,]
gset=gset[-1,]
a <- gset
a= a[1:nrow(a),]
set <- colnames(a)
l <- list()

for (i in set) { 
  x <- as.character(a[,i]) 
  x <- x[nchar(x)!=0] 
  x <- as.character(x) 
  l[[i]] <-x }

gsva_data = gsva(dc,l, method = "ssgsea")

b = gsva_data %>% t() %>% as.data.frame()

write.table(a,"ssGSEA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

library(NMF)
a = read.delim('ssGSEA.txt',row.names = 1,check.names = F)
ranks = 2:8
normalize <- function(x){
  return((x-min(x))/(max(x)-min(x)))
  }
dc = normalize(t(a))
result = nmf(dc,ranks, nrun=50,seed=123,.options=list(keep.all=TRUE))

pdf('./plot/01.rankSelection.pdf',width=7,height = 6)
plot(result)
dev.off()

pdf('./plot/02.consensusmap.pdf',width = 15,height = 15,onefile = F)
consensusmap(result,
             
             annRow = NA,
             
             annCol = NA,
             
             main = "Consensus matrix",
             
             info = FALSE)
dev.off() 



estimate2 <- nmf(dc,rank = 2,seed = 123,method = "brunet")
group <- predict(estimate2)
group <- as.data.frame(group) 



write.csv(group,file = 'subtype.csv')
group = read.csv('subtype.csv',row.names = 1)

#作图PCA
group = read.csv('subtype.csv',row.names = 1)
group$group = ifelse(group$group==1,'subtype1','subtype2')
pcadata <- as.data.frame(t(exp))
com1 <- prcomp(pcadata, center = TRUE,scale. = TRUE)
pcas <- data.frame(com1$x, Group = group$group) 

pdf('./plot/03.PCA_2D.pdf',width = 7,height = 6)
ggplot(pcas,aes(x=PC1,y=PC2,color=Group))+ geom_point()+scale_color_manual(values=c('#BC3C29FF','#0072B5FF'))+theme_bw()
dev.off()

library(FactoMineR)
library(factoextra)
library(scatterplot3d)

fviz_screeplot(com1, addlabels = TRUE, ylim = c(0, 30))
dat.pca <- com1$x
my_color = pal_nejm()(2)
colors =  my_color[as.numeric(group$group)]
pdf('./plot/03.PCA.pdf',width = 6,height = 6)
p3 = scatterplot3d(dat.pca[,1:3],color = "black",pch = 21,bg=colors)
group$group = ifelse(group$group==1,'subtype1','subtype2')
group$group <- as.factor(group$group)
legend("top",col = "black", legend = levels(group$group),pt.bg =  my_color, pch = 21,
       inset = -0.2, xpd = TRUE, horiz = TRUE)
dev.off()
