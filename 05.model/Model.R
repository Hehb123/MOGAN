library(tibble)
library(data.table)
library(caret)
fivegene = read.delim('./fivegene.txt')
fivegene = fivegene$x
exp = fread('../07.DEGs/exp.csv')
exp = column_to_rownames(exp,'V1')

exp = exp[fivegene,]
exp = as.data.frame(t(exp))
subtype = read.csv('../01.subtype/subtype.csv')
exp = exp[subtype$X,]
identical(rownames(exp),subtype$X)
exp = as.data.frame(scale(exp))
exp$subtype = subtype$group
exp$subtype = ifelse(exp$subtype==1,'subtype1','subtype2')
exp$subtype = as.factor(exp$subtype)

set.seed(123)
trainIndex = createDataPartition(exp$subtype, p = 0.6, list = FALSE)
trainData = exp[trainIndex, ]
testData = exp[-trainIndex, ]

ctrl = trainControl(method="repeatedcv", number=5, repeats=10, search="random")
fit = train(subtype ~ ., data = trainData, method = "glm", family = "binomial",trControl=ctrl)
predictions = predict(fit, newdata = testData)

save(fit,file= 'fit.rda')





library(caret)
library(ROCR)

# 假设已经拟合好的模型存储在fit对象�?

# 预测概率
prob <- predict(fit, testData, type = "prob")[,'subtype1']

# 真实类别标签
actual <- testData$subtype
actual = as.character(actual)

# 设置新的分类阈�?
new_threshold <- 0.4

# 根据新的分类阈值生成预测结�?
new_pred <- ifelse(prob > new_threshold, "subtype1", "subtype2")
table(actual,new_pred)[1]+table(actual,new_pred)[4]

tp = sum(new_pred == 'subtype1' &testData$subtype == 'subtype1')
fp = sum(new_pred == 'subtype1' &testData$subtype == 'subtype2')
fn = sum(new_pred == 'subtype2' &testData$subtype == 'subtype1')
tn = sum(new_pred == 'subtype2' &testData$subtype == 'subtype2')
precision = tp / (tp + fp)
recall = tp / (tp + fn)
f1_score = 2 * precision * recall / (precision + recall)
accuracy = sum(new_pred == testData$subtype) / nrow(testData)

predictions = predict(fit, testData, type = "prob")[,'subtype1']
labels = ifelse(testData$subtype=='subtype1',1,0)
pred <- prediction(predictions = predictions, 
                   labels = labels)
perf <- performance(prediction.obj = pred,
                    measure = "tpr",
                    x.measure = "fpr")
plot(perf,colorize=TRUE,
     main="ROCR fingerpainting toolkit",
     xlab="False Positive Rate", ylab="True Positive Rate", 
     box.lty=7, box.lwd=2, box.col="gray",linewidth=10)
auc=performance(pred,"auc")
auc_area<-slot(auc,"y.values")[[1]]
auc_area<-round(auc_area,4)
#添加文本注释
text_auc<-paste("AUC=", auc_area,sep="")
pdf('AUC.pdf',width = 6,height = 6)
plot(perf,colorize=TRUE,
     main="ROCR fingerpainting toolkit",
     xlab="False Positive Rate", ylab="True Positive Rate", 
     box.lty=7, box.lwd=2, box.col="gray",linewidth=10)
text(0.25,0.9,text_auc)
dev.off()