
#鐢熷瓨鏇茬嚎鍏嶇柅寰幆澧?
a = read.delim('../01.subtype/ssGSEA.txt',check.names = F)

surv = read.delim('../01.subtype/survival_LUAD_survival.txt')
surv = surv[,c('sample','OS','OS.time')]
surv = surv[surv$sample %in% rownames(a),]
identical(surv$sample,rownames(a))
surv = cbind(surv,a)

for(i in colnames(surv)[4:ncol(surv)]){
  df = data.frame(OS=surv$OS,OS.time=surv$OS.time,i=surv[i])
  df$group = ifelse(df[,3]>=median(df[,3]),'high','low')
  fit = survfit(Surv(OS.time,OS) ~ group,data=df)
  a = ggsurvplot(fit, pval = TRUE,
                 palette='nejm',
                 ggtheme = theme_bw(),
                 legend.title = i,
                 legend.labs = c("High Fraction", "Low Fraction"))
  ggsave(paste0('./plot/',i,'.pdf'),a$plot,width = 6,height = 6)
}