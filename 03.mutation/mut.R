#突变分析
library(tidyverse)
library(maftools)
library(ggplotify)
library(ggpubr)

result = read.csv('../基线资料/clini.csv',row.names = 1)
subtype1 = result[result$group=='subtype1',]$sampleID
subtype2 = result[result$group=='subtype2',]$sampleID


load('TCGA-LUAD_SNP.Rdata')
data = data[substring(data$Tumor_Sample_Barcode,14,16)=='01A',]
subtype_maf = data[substring(data$Tumor_Sample_Barcode,1,15) %in% c(subtype1,subtype2),]
maf = read.maf(subtype_maf,)

vc_cols <- RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) <- c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
a = maf@clinical.data
a$sampleID = substring(a$Tumor_Sample_Barcode,1,15)
a = merge(a,result,by='sampleID')
a= a[,-1]
a= a[order(a$group),]
maf@clinical.data=a
sample = maf@clinical.data$Tumor_Sample_Barcode

groupcol = c('#FF4500','#87CEFA')
names(groupcol)=c('subtype1','subtype2')

gendercol = c('#82B29A','#F2CC8E')
names(gendercol)=c('MALE','FEMALE')

stagecol = c('#000000','#3A0964','#7A1B6D','#BD3752','#ED6825')
names(stagecol) = c('Unknow','IV','III','II','I')

smokecol = c('#DAA520','#AFEEEE')
names(smokecol) = c('smoked','non-smoked')

agecol =c('#B0E0E6','#4682B4','#E6E6FA')
names(agecol) = c('<65','>=65','Unknown')

anno_col = list(group=groupcol,gender=gendercol,pathologic_stage=stagecol,smoke_status=smokecol,age=agecol)
pdf('../01.subtype/plot/07.oncoplot.pdf',width = 9,height = 7.5)
oncoplot(maf = maf,top = 20,clinicalFeatures = c('group'),colors = vc_cols,annotationColor = anno_col)
dev.off()

pdf('../01.subtype/plot/07.oncoplot_clini.pdf',width = 9,height = 7.5)
oncoplot(maf = maf,top = 20,clinicalFeatures = c('group','pathologic_stage','gender','smoke_status','age'),colors = vc_cols,annotationColor = anno_col,sortByAnnotation =T)
dev.off()

data_1 = subtype_maf[substring(subtype_maf$Tumor_Sample_Barcode,1,15)%in%subtype1,]
data_2 = subtype_maf[substring(subtype_maf$Tumor_Sample_Barcode,1,15)%in%subtype2,]
maf_1 = read.maf(data_1)
maf_2 = read.maf(data_2)

#################TMB分析
alltmb = tmb(maf)
alltmb$sampleID = substring(alltmb$Tumor_Sample_Barcode,1,15)
alltmb = merge(alltmb,result,by='sampleID')

pdf('../01.subtype/plot/08.TMB.pdf',height = 6,width = 7)
ggboxplot(alltmb,x='group',y='total_perMB',add = 'jitter',color = 'group',size = 0.7,fill = '#CACACA',alpha=0.5)+
  stat_boxplot(geom = "errorbar",width=0.3,aes(color=group,ymin=NULL,ymax=NULL),size=0.7,position=position_dodge(0))+
  stat_compare_means()+
  scale_color_manual(values = c('#BC3C29','#0072B5'))+
  theme(panel.background = element_rect(fill = "#EBEBEB"))+
  theme_bw()
dev.off()

surv = read.delim('../02.prog/survival_LUAD_survival.txt')
surv = surv[,c(1,3,4)]
colnames(surv)[1] = 'sampleID'
alltmb = merge(alltmb,surv,by='sampleID')
alltmb = alltmb[,c(1,4,11,12)]
alltmb$group = ifelse(alltmb$total_perMB>=median(alltmb$total_perMB),'high','low')
table(alltmb$group)
library(survival)
library(survminer)
fit<- survfit(Surv(OS.time, OS) ~ group, data = alltmb)
pdf('../01.subtype/plot/08.TMBsurv.pdf',width = 8,height = 6)
ggsurvplot(fit, pval = TRUE,risk.table = TRUE, risk.table.y.text.col = T,palette='nejm')
dev.off()



library(reshape2)
Chitest = function(maf1,maf2,gene){
  total_sample_1 = nrow(maf1@clinical.data)
  total_sample_2 = nrow(maf2@clinical.data)
  mut_info1 = maf1@gene.summary
  mut_info2 = maf2@gene.summary
  mut_num1 = mut_info1[mut_info1$Hugo_Symbol==gene,]$MutatedSamples
  mut_num2 = mut_info2[mut_info2$Hugo_Symbol==gene,]$MutatedSamples
  non_mut_num1 = total_sample_1-mut_num1
  non_mut_num2 = total_sample_2-mut_num2
  df = data.frame(subtype1=c(mut_num1,non_mut_num1),
                  subtype2=c(mut_num2,non_mut_num2),row.names = c('mut','non_mut'))
  pvalue = chisq.test(df,correct = F)
  return(pvalue)
}
genes = c('TP53','TTN','CSMD3','MUC16','RYR2','USH2A','ZFHX4','LRP1B','FLG',
          'KRAS','XIRP2','SPTA1','COL11A1','MUC17','ADAMTS12','ANK2','PTPRD','TNR','PAPPA2','PCLO')
df = data.frame(gene=genes,pvalue=NA)
for(i in genes){
  pvalue = Chitest(maf_1,maf_2,i)[["p.value"]]
  df[df$gene==i,]$pvalue=pvalue
}
df = df[df$pvalue<=0.001,]


stackplot = function(mafdata){
  N = nrow(mafdata@clinical.data)
  mut_info = mafdata@gene.summary
  mut_info = mut_info[,c(1,12)]
  mut_info$Frequency = mut_info$MutatedSamples/N
  colnames(mut_info)[1] = 'Gene'
  mut_info = mut_info[c(1:20),]
  mut_info$Gene = factor(mut_info$Gene,levels = mut_info$Gene)
  ggplot(data=mut_info,aes(Gene,Frequency,fill='#2272A9'))+
    geom_col()+
    theme_classic()+
    guides(fill=F)+
    scale_fill_manual(values=c("#2272A9"))+
    theme(panel.background=element_rect(fill="white",colour="black",linewidth=0.5))+
    ylim(0,0.6)
    
}
stackplot_2 = function(mafdata){
  mut_info = mafdata@gene.summary
  mut_info = mut_info[order(mut_info$total,decreasing = T),]
  mut_info = mut_info[,c(1:10)]
  mut_info = column_to_rownames(mut_info,'Hugo_Symbol')
  mut_info = as.data.frame(t(mut_info))
  mut_info = mut_info[,c(1:20)]
  mut_info = rownames_to_column(mut_info,'mutation_type')
  data =melt(mut_info,id.vars = 'mutation_type')
  colnames(data)[2:3] = c('Gene','Frequency')
  ggplot(data=data,aes(Gene,Frequency,fill=mutation_type))+
    geom_bar(stat="identity",position="stack", color="black", width=0.7,size=0.25)+
    scale_fill_manual(values=c("#A6CEE3","#FB9A99", "#FF7F00", "#E31A1C","#1F78B4" ,"#B2DF8A" ,"#33A02C" ,"#FDBF6F", "#CAB2D6"))+
    theme_classic()+
    theme(panel.background=element_rect(fill="white",colour="black",linewidth=0.5))
  
}
pdf('../01.subtype/plot/11.mut_frequency_subtype1.pdf',width=12,height = 6)
stackplot(maf_1)
dev.off()

pdf('../01.subtype/plot/11.mut_frequency_subtype2.pdf',width=12,height = 6)
stackplot(maf_2)
dev.off()
