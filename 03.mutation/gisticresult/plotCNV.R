BiocManager::install("PoisonAlien/maftools")
library(maftools)
subtype1.gistic = readGistic(gisticAllLesionsFile = 'subtype1_conf_99.txt',
                         gisticAmpGenesFile = 'subtype1_amp_conf_99.txt', 
                         gisticDelGenesFile = 'subtype1_del_conf_99.txt', 
                         gisticScoresFile = 'subtype1_scores.gistic', 
                         isTCGA = T)

gisticChromPlot(gistic = subtype1.gistic, ref.build = 'hg38')
S1_info = subtype1.gistic@cnv.summary
S1_info = S1_info[,c(1,4)]
S1_info$group = 'subtype1'

subtype2.gistic = readGistic(gisticAllLesionsFile = 'subtype2_conf_99.txt',
                         gisticAmpGenesFile = 'subtype2_amp_conf_99.txt', 
                         gisticDelGenesFile = 'subtype2_del_conf_99.txt', 
                         gisticScoresFile = 'subtype2_scores.gistic', 
                         isTCGA = T)
S2_info = subtype2.gistic@cnv.summary
S2_info = S2_info[,c(1,4)]
S2_info$group = 'subtype2'
all_info = rbind(S1_info,S2_info)
all_info$`log10(CNV+1)` = log10(all_info$total+1)

library(ggpubr)
pdf('./CNVplot/CNV.pdf',height = 6,width = 7)
ggboxplot(all_info,x='group',y="`log10(CNV+1)`",add = 'jitter',color = 'group',size = 0.7,fill = '#CACACA',alpha=0.5)+
  stat_boxplot(geom = "errorbar",width=0.3,aes(color=group,ymin=NULL,ymax=NULL),size=0.7,position=position_dodge(0))+
  stat_compare_means()+
  scale_color_manual(values = c('#BC3C29','#0072B5'))+
  theme(panel.background = element_rect(fill = "#EBEBEB"))+
  theme_bw()
dev.off()


pdf('./CNVplot/subtype2_CNV.pdf',width = 8,height = 3)
gisticChromPlot(gistic = subtype2.gistic, ref.build = 'hg38',fdrCutOff = 0.05,mutGenes = T)
dev.off()
pdf('./CNVplot/subtype1_CNV.pdf',width = 8,height = 3)
gisticChromPlot(gistic = subtype1.gistic, ref.build = 'hg38',fdrCutOff = 0.05,mutGenes = T)
dev.off()

#基因查询
# 安装和加载所需的包
library(biomaRt)  # 加载biomaRt包

# 连接到Ensembl数据库
ensembl = useMart('ensembl',dataset = "hsapiens_gene_ensembl",host = "https://asia.ensembl.org/")

# 指定染色质区域
chromosome = "1"
arm = "p"
region = "13.2"

# 获取染色质区域上的基因信息
genes = getBM(attributes = c("ensembl_gene_id", "external_gene_name", "start_position", "end_position"),
              filters = c("chromosome_name", "chromosomal_region"),
              values = list(chromosome_name = chromosome, chromosomal_region = paste(arm, region, sep = "")),
              mart = ensembl)

# 打印基因信息
print(genes)

