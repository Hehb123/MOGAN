library(SummarizedExperiment)
library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-LUAD", 
  data.category = "Copy Number Variation",
  data.type = "Masked Copy Number Segment",
)

GDCdownload(query)

GDCprepare(query, save = T,save.filename = "TCGA-LUAD_CNV.Rdata")
a = load(file = "./TCGA-LUAD_CNV.Rdata")
tumorCNV  = eval(parse(text = a))
tumorCNV = tumorCNV[,2:7]
tumorCNV = tumorCNV[,c('Sample','Chromosome','Start','End','Num_Probes',"Segment_Mean")]

rt  = read.csv('../01.subtype/subtype.csv')
subtype_1 = rt[rt$group==1,]$X
subtype_2 = rt[rt$group==2,]$X
subtype1_CNV = tumorCNV[substring(tumorCNV$Sample,1,15) %in% subtype_1,]
subtype2_CNV = tumorCNV[substring(tumorCNV$Sample,1,15) %in% subtype_2,]
write.table(subtype1_CNV,file = './gisticdata/subtype1_CNV.txt',sep = '\t',quote = F,row.names = F)
write.table(subtype2_CNV,file = './gisticdata/subtype2_CNV.txt',sep = '\t',quote = F,row.names = F)

marker_file = read.delim('./gisticdata/snp6.na35.remap.hg38.subset.txt')
marker_file = marker_file[marker_file$freqcnv==FALSE,]
marker_file = marker_file[,1:3]
colnames(marker_file) = c('Marker Name','Chromosome','Marker Position')
write.table(marker_file,file = './gisticdata/marker_file.txt',sep = '\t',quote = F,row.names = F)
