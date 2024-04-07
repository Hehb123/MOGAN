library(ggpubr)

MSI_info = read.delim('luad_tcga_pan_can_atlas_2018_clinical_data.tsv',check.names = F)
MSI_info = MSI_info[,c(2,6)]
subtype = read.csv('../01.subtype/subtype.csv')
subtype$group = ifelse(subtype$group==1,'subtype1','subtype2')
colnames(subtype)[1] = 'Sample ID'

rt = merge(MSI_info,subtype,by='Sample ID')

pdf('./plot/01.MSI_diff.pdf',width = 6,height = 6)
ggboxplot(rt, x = "group", y = "MSI MANTIS Score",
          fill = "group") +
  stat_compare_means(label.y = 0.3,label.x=1.42,method = 't.test')
dev.off()

MSI_info_1 = rt[rt$group=='subtype1',]
MSI_info_2 = rt[rt$group=='subtype2',]
ggdensity(MSI_info_2$`MSI MANTIS Score`, 
          main = "Density plot of sepal length",
          xlab = "sepal length")
ggqqplot(MSI_info_1$`MSI MANTIS Score`)
shapiro.test(MSI_info_2$`MSI MANTIS Score`)
