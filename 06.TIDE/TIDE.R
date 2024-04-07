library(ggplot2)
library(ggpubr)
exp  = read.csv('../01.subtype/omics1_4.csv',row.names = 1)
info = read.csv('../01.subtype/subtype.csv')
info$group = ifelse(info$group==1,'subtype1','subtype2')
subtype1 = info[info$group=='subtype1',]$X
subtype2 = info[info$group=='subtype2',]$X

data = as.data.frame(t(scale(exp)))
subtype1_expr = data[,colnames(data) %in% subtype1]
subtype2_expr = data[,colnames(data) %in% subtype2]


write.table(subtype1_expr, file='subtype1_expr.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(subtype2_expr, file='subtype2_expr.tsv', quote=FALSE, sep='\t', col.names = NA)



subtype1 <- read.csv('./result/subtype1_result.csv')
subtype2 <- read.csv('./result/subtype2_result.csv')

subtype2TIDE <-subtype2$TIDE
subtype1TIDE <- subtype1$TIDE

subtype2$group <- 'subtype2'
subtype1$group <- 'subtype1'

all <- rbind(subtype2,subtype1)
all = all[order(all$group),]
all$group = ifelse(all$group=='subtype1','S1','S2')
pdf('../01.subtype/plot/06.TIDE.pdf',width = 7,height = 6)
ggboxplot(all,x='group',y='TIDE',add = 'jitter',color = 'group',size = 0.7,fill = '#CACACA',alpha=0.5)+
  stat_boxplot(geom = "errorbar",width=0.3,aes(color=group,ymin=NULL,ymax=NULL),size=0.7,position=position_dodge(0))+
  stat_compare_means()+
  scale_color_manual(values = c('#BC3C29','#0072B5'))+
  theme(panel.background = element_rect(fill = "#EBEBEB"))+
  theme_bw()

dev.off()