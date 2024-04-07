library(survival)
library(survminer)
###生存分析
group = read.csv('subtype.csv',row.names = 1)
group$group = ifelse(group$group==1,'subtype1','subtype2')
table(group$group)
surv = read.delim('survival_LUAD_survival.txt')
surv = surv[,c('sample','OS','OS.time')]

surv = surv[surv$sample %in% rownames(group),]
identical(rownames(group),surv$sample)
surv = cbind(surv,group)

fit<- survfit(Surv(OS.time, OS) ~ group, data = surv)
pdf('../01.subtype/plot/04.surv.pdf',width = 8,height = 6)
ggsurvplot(fit, pval = TRUE,risk.table = TRUE, risk.table.y.text.col = T,palette='nejm')
dev.off()
