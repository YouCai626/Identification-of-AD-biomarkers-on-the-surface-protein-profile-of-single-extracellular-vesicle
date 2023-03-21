# Set the maximum memory
if(.Platform$OS.type == "windows") withAutoprint({
  memory.size()
  memory.size(TRUE)
  memory.limit()
})

library(tidyverse) 
library(caret)
library(pROC)
library(pheatmap)

# -- ------------------------------------------------------------------------------------------------------------------
# 1.Pseudo processing of samples
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
sampleID <- read.csv('Human_SampleID_SampleName.csv')
head(sampleID)
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('RawData/',sampleID[i,1],'.total_ev_protein.csv')
  data <- read.csv(fn1, header = T)
  wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
  rownames(wide_data) <- wide_data$ev
  wide_data <- wide_data[, -c(1:2)]
  pseudoSize10 <- floor(nrow(wide_data)/10) # 一个raw转换为10个sample
  set.seed(2022)
  random.row <- sample(1:nrow(wide_data), pseudoSize10*10, replace = F)
  raw.data <- wide_data[random.row,] %>% as.data.frame()
  raw.data$Group <- rep(c(1:10), each = pseudoSize10)
  raw.data <- relocate(raw.data, Group)
  pseudo.data <- aggregate(list(raw.data[2:ncol(raw.data)]), by = list(raw.data$Group), FUN = sum)
  rownames(pseudo.data) <- str_c(sampleID[i,2], 'sEV', 1:10, sep = '_')
  fn2 <- str_c('pseudoData/pseudo10/', sampleID[i,1], '.pseudo10.csv')
  write.csv(pseudo.data[-1], fn2)
  rm(data, wide_data,raw.data, pseudo.data)
}

# -- ------------------------------------------------------------------------------------------------------------------
# 2.Correlation of pseudo samples with raw
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read_csv('Results/Human_Raw.Counts.csv', show_col_types = F)
names(data)[1] <- 'Pro'
data.anno <- read_csv('Human_SampleID_SampleName.csv', show_col_types = F)
for (i in 1:nrow(data.anno)) {
  dt <- data %>% select(Pro, data.anno$SampleName[i])
  fn1 <- str_c('pseudoData/pseudo10/', data.anno$SampleID[i], '.pseudo10.csv')
  temp <- read.csv(fn1, row.names = 1)
  temp <- t(temp) %>% as.data.frame()
  temp$Pro <- rownames(temp)
  dt <- inner_join(dt, temp, by = 'Pro')
  dt.cor <- cor(dt[-1])
  fn2 <- str_c('Results/FigureS6/', data.anno$SampleName[i], '_Raw.pseudo10.png')
  pheatmap(dt.cor, color = pal_material(palette = c("orange"), n = 20)(20), cellwidth = 10, cellheight = 6,
           display_numbers = T, number_color = 'black', number_format = '%.3f', 
           fontsize_number = 3, fontsize_row = 6, fontsize_col = 8,
           cluster_rows = F, cluster_cols = F,
           filename = fn2) 
  fn2 <- str_c('Results/FigureS6/', data.anno$SampleName[i], '_Raw.pseudo10.pdf')
  pheatmap(dt.cor, color = pal_material(palette = c("orange"), n = 20)(20), cellwidth = 10, cellheight = 6,
           display_numbers = T, number_color = 'black', number_format = '%.3f', 
           fontsize_number = 3, fontsize_row = 6, fontsize_col = 8,
           cluster_rows = F, cluster_cols = F,
           filename = fn2) 
}

# -- ------------------------------------------------------------------------------------------------------------------
# 3.Machine learning
# -- ------------------------------------------------------------------------------------------------------------------
library(edgeR)
library(preprocessCore)
# -- ------------------------------------------------------------------------------------------------------------------
# HumanUrine
# -- ------------------------------------------------------------------------------------------------------------------
# total expression matrix
rm(list=ls())
gc()
data.anno <- read_csv('Results/HumanUrine_Samples.Annotation.csv', show_col_types = F)
fn1 <- str_c('pseudoData/pseudo10/', data.anno$SampleID[1], '.pseudo10.csv')
dt <- read.csv(fn1, row.names = 1)
dt <- dt %>% t() %>% as.data.frame()
dt$Pro <- rownames(dt)
dt <- relocate(dt, Pro)
for (j in 2:nrow(data.anno)) {
  fn2 <- str_c('pseudoData/pseudo10/', data.anno$SampleID[j], '.pseudo10.csv')
  dt2 <- read.csv(fn2, row.names = 1)
  dt2 <- dt2 %>% t() %>% as.data.frame()
  dt2$mean <- rowMeans(dt2)
  dt2 <- dt2 %>% filter(mean >= 1) %>% dplyr::select(!mean)
  dt2$Pro <- rownames(dt2)
  dt <- inner_join(dt, dt2, by = 'Pro')
}
write_csv(dt, 'pseudoData/HumanUrine.AllSample.pseudo10.csv')
# Standardization of data and annotation of sample groups
rm(list=ls())
gc()
data <- read.csv('pseudoData/HumanUrine.AllSample.pseudo10.csv', row.names = 1)
data.nor <- data %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(data.nor) <- rownames(data)
names(data.nor) <- names(data)
data.nor$mean <- rowMeans(data.nor)
data.nor <- data.nor %>% filter(mean >= 10) %>% dplyr::select(!mean) # Removal of Pro with too low expression
data <- data.nor %>% t() %>% as.data.frame() 
ID <- rownames(data) %>% str_split('_s') %>% as.data.frame() %>% t() %>% as.data.frame()
data$ID <- ID$V1 %>% factor
data$Group <- ifelse(str_detect(data$ID, 'D'), 1, 0)
data$Group <- factor(data$Group, levels = c(1,0), labels = c('AD', 'Con'))
data <- relocate(data, ID, Group)
write.csv(data, 'pseudoData/HumanUrine.AllSample.pseudo10.CPM.quantiles.Group.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# Machine learning with the caret package
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('pseudoData/HumanUrine.AllSample.pseudo10.CPM.quantiles.Group.csv', row.names = 1)
data$ID <- data$ID %>% factor()
data$Group <- data$Group %>% factor()
# -- ------------------------------------------------------------------------------------------------------------------
# a.Dividing the training set and the test set:
# -- ------------------------------------------------------------------------------------------------------------------
data.ID <- levels(data$ID) 
data.ID <- data.frame(ID = data.ID)
data.ID$Group <- ifelse(str_detect(data.ID$ID, 'D'), 'AD', 'Con')
set.seed(2022)
train.ID <- createDataPartition(data.ID$Group, p = 0.6, list = FALSE)
train.ID <- data.ID$ID[train.ID]
data.train <- data %>% filter(ID %in% train.ID) %>% dplyr::select(!ID)
data.test <- data %>% filter(!ID %in% train.ID) %>% dplyr::select(!ID)
# -- ------------------------------------------------------------------------------------------------------------------
# b.Filtering variables with LASSO
# -- ------------------------------------------------------------------------------------------------------------------
X <- as.matrix(data.train[2:ncol(data.train)])
Y <- as.matrix(data.train[1])
myLasso <- glmnet(X, Y, alpha = 1, family = 'binomial', nlambda = 200)
NROW(myLasso$lambda)
min(myLasso$lambda)
pdf('Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO.pdf', width = 6, height = 8)
plot(myLasso, xvar = 'lambda', label = T)
dev.off()
lasso.coef <- coef(myLasso, s = min(myLasso$lambda))
my.lasso.coef <- as.matrix(lasso.coef) %>% as.data.frame()
my.lasso.coef$abs <- abs(my.lasso.coef$s1)
my.lasso.coef <- arrange(my.lasso.coef, desc(abs)) %>% filter(abs > 0)
write.csv(my.lasso.coef, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO.csv')
my.var <- rownames(my.lasso.coef)
my.var <- my.var[-1]
res <- c()
for (i in 1:length(my.var)) {
  fit <- glm(substitute(Group ~ x, list(x = as.name(my.var[i]))), family = binomial(), data = data.train)
  coef.data <- summary(fit)$coef
  or.data <- exp(cbind('OR' = coef(fit), confint(fit)))
  res1 <- cbind(coef.data, or.data) %>% as.data.frame()
  rownames(res1)[1] <- str_c('Intercept', my.var[i], sep = '.')
  res <- rbind(res, res1[-1,]) %>% as.data.frame()
}
res1 <- res %>% arrange(`Pr(>|z|)`)
write.csv(res1, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.csv')
res2 <- res1 %>% filter(`Pr(>|z|)` < 0.05)
write.csv(res2, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# c.Construction of Train and Test datasets with positive variables
# -- ------------------------------------------------------------------------------------------------------------------
SigPro <- read.csv('Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos.csv', row.names = 1)
myPro <- str_replace_all(rownames(SigPro)[1:20], '/', '.') %>% str_replace_all('-', '.')
myPro <- intersect(myPro, names(data.train))
data.train <- data.train %>% dplyr::select(Group, all_of(myPro))
data.test <- data.test %>% dplyr::select(Group, all_of(myPro))
# Normalize the data and fill in the missing values, this can be done with the preProcess function
preProcValues <- preProcess(data.train, method = c('center','scale'))
train.transformed <- predict(preProcValues, data.train)
test.transformed <- predict(preProcValues, data.test)
str(train.transformed)
# Save data
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')
# -- ------------------------------------------------------------------------------------------------------------------
# d.Model Training and Tuning
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')
myMethod <- c('glmnet','kknn','C5.0','rf','AdaBag','gbm','nnet','svmLinear','svmPoly','svmRadial','nb')
fitControl = trainControl(method = 'cv', number = 5,
                          classProbs = T, summaryFunction = twoClassSummary,
                          search = 'random')

for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed, method = myMethod[i],
               trControl = fitControl, 
               tuneLength = 30, 
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
# -- ------------------------------------------------------------------------------------------------------------------
# e.Model prediction and evaluation
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for ( a in 1:11) {
  Fit <- readRDS(filenames[a])
  predict.train <- predict(Fit, newdata = train.transformed)
  predict.test <- predict(Fit, newdata = test.transformed)
  # The ROC analysis should distinguish the direction
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test) %>% mutate(method = myMethod[a])
  rownames(Fit.roc) <- myMethod[a]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}

# -- ------------------------------------------------------------------------------------------------------------------
# 4.The AUC predicted by the 3 randomized experiments were presented in a heat map
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data1 <- read.csv('Results/Figure5/Files/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_1_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data2 <- read.csv('Results/Figure5/Files/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data3 <- read.csv('Results/Figure5/Files/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_3_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data <- left_join(data1, data2, by = 'method')
data <- left_join(data, data3, by = 'method')
rownames(data) <- data$method
data.Urine <- data %>% dplyr::select(method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data.Urine$Train.Mean <- rowMeans(data.Urine[2:4])
data.Urine$Train.sd <- apply(data.Urine[2:4], 1, sd)
data.Urine$Test.Mean <- rowMeans(data.Urine[5:7])
data.Urine$Test.sd <- apply(data.Urine[5:7], 1, sd)
data.Urine.heatmap <- data.Urine %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)

data1 <- read.csv('Results/Figure5/Files/HumanPlasma.AllSample.pseudo10.CPM.quantiles.train_1_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data2 <- read.csv('Results/Figure5/Files/HumanPlasma.AllSample.pseudo10.CPM.quantiles.train_2_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data3 <- read.csv('Results/Figure5/Files/HumanPlasma.AllSample.pseudo10.CPM.quantiles.train_3_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data <- left_join(data1, data2, by = 'method')
data <- left_join(data, data3, by = 'method')
rownames(data) <- data$method
data.Plasma <- data %>% dplyr::select(method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data.Plasma$Train.Mean <- rowMeans(data.Plasma[2:4])
data.Plasma$Train.sd <- apply(data.Plasma[2:4], 1, sd)
data.Plasma$Test.Mean <- rowMeans(data.Plasma[5:7])
data.Plasma$Test.sd <- apply(data.Plasma[5:7], 1, sd)
data.Plasma.heatmap <- data.Plasma %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)

data1 <- read.csv('Results/Figure5/Files/HumanSaliva.AllSample.pseudo10.CPM.quantiles.train_1_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data2 <- read.csv('Results/Figure5/Files/HumanSaliva.AllSample.pseudo10.CPM.quantiles.train_2_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data3 <- read.csv('Results/Figure5/Files/HumanSaliva.AllSample.pseudo10.CPM.quantiles.train_3_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data <- left_join(data1, data2, by = 'method')
data <- left_join(data, data3, by = 'method')
rownames(data) <- data$method
data.Saliva <- data %>% dplyr::select(method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data.Saliva$Train.Mean <- rowMeans(data.Saliva[2:4])
data.Saliva$Train.sd <- apply(data.Saliva[2:4], 1, sd)
data.Saliva$Test.Mean <- rowMeans(data.Saliva[5:7])
data.Saliva$Test.sd <- apply(data.Saliva[5:7], 1, sd)
data.Saliva.heatmap <- data.Saliva %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)

data1 <- read.csv('Results/Figure5/Files/HumanTear.AllSample.pseudo10.CPM.quantiles.train_1_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data2 <- read.csv('Results/Figure5/Files/HumanTear.AllSample.pseudo10.CPM.quantiles.train_2_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data3 <- read.csv('Results/Figure5/Files/HumanTear.AllSample.pseudo10.CPM.quantiles.train_3_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data <- left_join(data1, data2, by = 'method')
data <- left_join(data, data3, by = 'method')
rownames(data) <- data$method
data.Tear <- data %>% dplyr::select(method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data.Tear$Train.Mean <- rowMeans(data.Tear[2:4])
data.Tear$Train.sd <- apply(data.Tear[2:4], 1, sd)
data.Tear$Test.Mean <- rowMeans(data.Tear[5:7])
data.Tear$Test.sd <- apply(data.Tear[5:7], 1, sd)
data.Tear.heatmap <- data.Tear %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)
# Plot
data.heatmap <- rbind(data.Urine.heatmap, data.Saliva.heatmap, data.Tear.heatmap, data.Plasma.heatmap)
data.heatmap <- data.heatmap %>% round(., 3)
data.heatmap$Train.Show <- str_c(data.heatmap$Train.Mean, '±', data.heatmap$Train.sd, sep = ' ')
data.heatmap$Test.Show <- str_c(data.heatmap$Test.Mean, '±', data.heatmap$Test.sd, sep = ' ')
data.heatmap[1:2] %>% range()
breakList <- seq(0.45, 1, 0.05)
show_col(pal_material('blue')(12))
colorRampPalette(c('#90CAF8FF', '#FFF2DFFF', '#FF9800FF'))(12) %>% show_col()
pheatmap(data.heatmap[1:2], colorRampPalette(c('#90CAF8FF', '#FFF2DFFF', '#FF9800FF'))(12), 
         cellwidth = 40, cellheight = 10,
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6], 
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(11, 22, 33),
         filename = 'Results/Figure5/Human.pseudo10_caret.ROC 2.png') 
pheatmap(data.heatmap[1:2], colorRampPalette(c('#90CAF8FF', '#FFF2DFFF', '#FF9800FF'))(12), 
         cellwidth = 40, cellheight = 10,
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6],
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(11, 22, 33),
         filename = 'Results/Figure5/Human.pseudo10_caret.ROC 2.pdf') 
pheatmap(data.heatmap[1:2], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6],
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(11, 22, 33),
         filename = 'Results/Figure5/Human.pseudo10_caret.ROC.png')
pheatmap(data.heatmap[1:2], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6],
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(11, 22, 33),
         filename = 'Results/Figure5/Human.pseudo10_caret.ROC.pdf')

# -- ------------------------------------------------------------------------------------------------------------------
# 5.Showing ROC plot and Feature ranking of the best model
# -- ------------------------------------------------------------------------------------------------------------------
# HumanUrine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure5/Files/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2_LASSO_Logistic.Pos_Top20_caret.RData')
Fit <- readRDS('Results/Figure5/Files/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2_LASSO_Logistic.Pos_Top20.nb.Fit.rds')
var.importance <- varImp(Fit, scale = T)
png('Results/Figure5/HumanUrine.pseudo10.CPM.quantiles.nb.Fit.2022.VarImportance.png', width = 300, height = 400)
plot(var.importance, xlab = 'Importance', main = 'nb')
dev.off()
pdf('Results/Figure5/HumanUrine.pseudo10.CPM.quantiles.nb.Fit.2022.VarImportance.pdf', width = 3, height = 4)
plot(var.importance, xlab = 'Importance', main = 'nb')
dev.off()
# ROC
predict.test <- predict(Fit, newdata = test.transformed, type = 'prob')
roc.test <- roc(test.transformed$Group, predict.test$AD)$auc
rocPlot <- roc(test.transformed$Group, predict.test$AD, levels = c('AD', 'Con'), direction = '>')
png('Results/Figure5/HumanUrine.pseudo10.CPM.quantiles.nb.Fit.ROC.png', width = 400, height = 400)
plot.roc(rocPlot, col="red", lwd = 4, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.5, auc.polygon = TRUE, auc.polygon.col = "#F0E685CC",
         print.auc.x=0.7, print.auc.y=0.1,
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
pdf('Results/Figure5/HumanUrine.pseudo10.CPM.quantiles.nb.Fit.ROC.pdf', width = 4, height = 4)
plot.roc(rocPlot, col="red", lwd = 4, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.5, auc.polygon = TRUE, auc.polygon.col = "#F0E685CC",
         print.auc.x=0.7, print.auc.y=0.1,
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()

# The analysis for other groups is similar to the urine group
