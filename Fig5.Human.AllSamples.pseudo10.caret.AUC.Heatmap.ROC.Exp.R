# -- ------------------------------------------------------------------------------------------------------------------
getwd()
# 设置最大的内存
if(.Platform$OS.type == "windows") withAutoprint({
  memory.size()
  memory.size(TRUE)
  memory.limit()
})
# options(future.globals.maxSize = 30*1024^3) #设置为30G

# -- ------------------------------------------------------------------------------------------------------------------
# 用caret包对压缩后的数据进行机器学习
library(tidyverse) # 处理数据神包
library(caret) # 机器学习集成包
library(pROC)
library(pheatmap)

# -- ------------------------------------------------------------------------------------------------------------------
# 1.对人的样本原始数据进行pseudo处理:1-10,1-20
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
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('RawData/',sampleID[i,1],'.total_ev_protein.csv')
  data <- read.csv(fn1, header = T)
  wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
  rownames(wide_data) <- wide_data$ev
  wide_data <- wide_data[, -c(1:2)]
  pseudoSize20 <- floor(nrow(wide_data)/20) # 一个raw转换为20个sample
  set.seed(2022)
  random.row <- sample(1:nrow(wide_data), pseudoSize20*20, replace = F)
  raw.data <- wide_data[random.row,] %>% as.data.frame()
  raw.data$Group <- rep(c(1:20), each = pseudoSize20)
  raw.data <- relocate(raw.data, Group)
  pseudo.data <- aggregate(list(raw.data[2:ncol(raw.data)]), by = list(raw.data$Group), FUN = sum)
  rownames(pseudo.data) <- str_c(sampleID[i,2], 'sEV', 1:20, sep = '_')
  fn2 <- str_c('pseudoData/pseudo20/', sampleID[i,1], '.pseudo20.csv')
  write.csv(pseudo.data[-1], fn2)
  rm(data, wide_data,raw.data, pseudo.data)
}
# 由于之前已经跑过这部分的内容了，所以直接把之前的结果COPY过来
rm(list=ls())
gc()
setwd("D:/2.AD.PBA.CAM.Article.V3.Proj")
sampleID <- read.csv('Human_SampleID_SampleName.csv')
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('D:/2.AD.PBA.CAM.Article.V1.Proj/pseudoData/pseudo10/',sampleID[i,1],'.pseudo10.csv')
  fn2 <- str_c('D:/2.AD.PBA.CAM.Article.V3.Proj/pseudoData/pseudo10/',sampleID[i,1],'.pseudo10.csv')
  file.copy(from = fn1, to = fn2)
}
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('D:/2.AD.PBA.CAM.Article.V1.Proj/pseudoData/pseudo20/',sampleID[i,1],'.pseudo20.csv')
  fn2 <- str_c('D:/2.AD.PBA.CAM.Article.V3.Proj/pseudoData/pseudo20/',sampleID[i,1],'.pseudo20.csv')
  file.copy(from = fn1, to = fn2)
}

# -- ------------------------------------------------------------------------------------------------------------------
# 2.pseudo数据与raw的correlation
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read_csv('Results/Human_Raw.Counts.csv', show_col_types = F)
names(data)[1] <- 'Pro'
data.anno <- read_csv('Human_SampleID_SampleName.csv', show_col_types = F)
psedoSize <- c(10,20)
for (i in 1:nrow(data.anno)) {
  for (j in psedoSize) {
    dt <- data %>% select(Pro, data.anno$SampleName[i])
    fn1 <- str_c('pseudoData/pseudo', j, '/', data.anno$SampleID[i], '.pseudo', j ,'.csv')
    temp <- read.csv(fn1, row.names = 1)
    temp <- t(temp) %>% as.data.frame()
    temp$Pro <- rownames(temp)
    dt <- inner_join(dt, temp, by = 'Pro')
    dt.cor <- cor(dt[-1])
    fn2 <- str_c('Results/FigureS6/', data.anno$SampleName[i], '_Raw.pseudo', j, '.png')
    pheatmap(dt.cor, color = pal_material(palette = c("orange"), n = 20)(20), cellwidth = 10, cellheight = 6,
             display_numbers = T, number_color = 'black', number_format = '%.3f', 
             fontsize_number = 3, fontsize_row = 6, fontsize_col = 8,
             cluster_rows = F, cluster_cols = F,
             filename = fn2) # 保存，自动调整纸张大小
    fn2 <- str_c('Results/FigureS6/', data.anno$SampleName[i], '_Raw.pseudo', j, '.pdf')
    pheatmap(dt.cor, color = pal_material(palette = c("orange"), n = 20)(20), cellwidth = 10, cellheight = 6,
             display_numbers = T, number_color = 'black', number_format = '%.3f', 
             fontsize_number = 3, fontsize_row = 6, fontsize_col = 8,
             cluster_rows = F, cluster_cols = F,
             filename = fn2) # 保存，自动调整纸张大小
  }
}

# -- ------------------------------------------------------------------------------------------------------------------
# 3.pseudo10数据进行机器学习
# -- ------------------------------------------------------------------------------------------------------------------
library(edgeR) # 用于cpm标准化
library(preprocessCore) # CPM.quantiles函数进行标准化
# -- ------------------------------------------------------------------------------------------------------------------
# HumanUrine
# -- ------------------------------------------------------------------------------------------------------------------
# 读取pseudo10文件合成总的表达矩阵
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
# 对数据进行标准化，并进行样本分组注释
rm(list=ls())
gc()
data <- read.csv('pseudoData/HumanUrine.AllSample.pseudo10.csv', row.names = 1)
data.nor <- data %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(data.nor) <- rownames(data)
names(data.nor) <- names(data)
data.nor$mean <- rowMeans(data.nor)
data.nor <- data.nor %>% filter(mean >= 10) %>% dplyr::select(!mean) # 将表达量太低的Pro去除
data <- data.nor %>% t() %>% as.data.frame() # 行为样本，列为基因的数据框
ID <- rownames(data) %>% str_split('_s') %>% as.data.frame() %>% t() %>% as.data.frame()
data$ID <- ID$V1 %>% factor
data$Group <- ifelse(str_detect(data$ID, 'D'), 1, 0)
data$Group <- factor(data$Group, levels = c(1,0), labels = c('AD', 'Con'))
data <- relocate(data, ID, Group)
write.csv(data, 'pseudoData/HumanUrine.AllSample.pseudo10.CPM.quantiles.Group.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# 用caret包进行机器学习
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('pseudoData/HumanUrine.AllSample.pseudo10.CPM.quantiles.Group.csv', row.names = 1)
data$ID <- data$ID %>% factor()
data$Group <- data$Group %>% factor()
# -- ------------------------------------------------------------------------------------------------------------------
# a.划分训练集与测试集:
# -- ------------------------------------------------------------------------------------------------------------------
data.ID <- levels(data$ID) # 只有factor才能用levels这个函数
data.ID <- data.frame(ID = data.ID)
data.ID$Group <- ifelse(str_detect(data.ID$ID, 'D'), 'AD', 'Con')
set.seed(1234)
set.seed(5678)
set.seed(2020)
set.seed(2021)
set.seed(2022)
train.ID <- createDataPartition(data.ID$Group, p = 0.6, list = FALSE)
train.ID <- data.ID$ID[train.ID]
data.train <- data %>% filter(ID %in% train.ID) %>% dplyr::select(!ID)
data.test <- data %>% filter(!ID %in% train.ID) %>% dplyr::select(!ID)
# -- ------------------------------------------------------------------------------------------------------------------
# b.用LASSO筛选变量
# -- ------------------------------------------------------------------------------------------------------------------
X <- as.matrix(data.train[2:ncol(data.train)])
Y <- as.matrix(data.train[1])
# glmnet()语法中alpha=0表示岭回归，1表示LASSO回归
myLasso <- glmnet(X, Y, alpha = 1, family = 'binomial', nlambda = 200) # glmnet默认运行100次
NROW(myLasso$lambda)
min(myLasso$lambda) # 选最优模型的lambda值
# 绘制图形
pdf('Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO.pdf', width = 6, height = 8)
pdf('Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO.pdf', width = 6, height = 8)
pdf('Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO.pdf', width = 6, height = 8)
pdf('Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO.pdf', width = 6, height = 8)
pdf('Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO.pdf', width = 6, height = 8)
plot(myLasso, xvar = 'lambda', label = T)
dev.off()
lasso.coef <- coef(myLasso, s = min(myLasso$lambda)) # 最优解下的回归系数
# 只有筛选出的自变量才有回归系数
# 将变量根据回归系数进行排序，回归系数越大的变量越重要
my.lasso.coef <- as.matrix(lasso.coef) %>% as.data.frame()
my.lasso.coef$abs <- abs(my.lasso.coef$s1)
my.lasso.coef <- arrange(my.lasso.coef, desc(abs)) %>% filter(abs > 0)
write.csv(my.lasso.coef, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO.csv')
# 用LASSO筛选的变量进行单因素Logistic回归，获得每一个因素的P值，结合LASSO筛选的相关系数获得变量的权重排序
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
write.csv(res1, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.csv')
res2 <- res1 %>% filter(`Pr(>|z|)` < 0.05)
write.csv(res2, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# c.用阳性变量构建Train和Test数据集
# -- ------------------------------------------------------------------------------------------------------------------
SigPro <- read.csv('Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos.csv', row.names = 1)
myPro <- str_replace_all(rownames(SigPro)[1:20], '/', '.') %>% str_replace_all('-', '.')
myPro <- intersect(myPro, names(data.train))
data.train <- data.train %>% dplyr::select(Group, all_of(myPro))
data.test <- data.test %>% dplyr::select(Group, all_of(myPro))
# 将数据进行标准化并补足缺失值，这时可以用preProcess命令，缺省参数是标准化数据。?preProcess
preProcValues <- preProcess(data.train, method = c('center','scale')) # 用于进行其他验证数据集的标准化
train.transformed <- predict(preProcValues, data.train)
test.transformed <- predict(preProcValues, data.test)
str(train.transformed)
# 保存数据
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')
# -- ------------------------------------------------------------------------------------------------------------------
# d.模型训练及调参
# -- ------------------------------------------------------------------------------------------------------------------
# 常用的分类算法包括：NBC（Naive Bayesian Classifier，朴素贝叶斯分类）算法、
# LR（Logistic Regress，逻辑回归）算法、C5.0 决策树算法、SVM（Support Vector Machine，支持向量机）算法、
# KNN(K-Nearest Neighbor，K 最近邻近)算法、ANN（Artificial Neural Network，人工神经网络）算法
# Stochastic Gradient Boosting(GBM)
rm(list=ls())
gc()
load('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')

myMethod <- c('glmnet','kknn','C5.0','rf','AdaBag','gbm','nnet','svmLinear','svmPoly','svmRadial','nb')
fitControl = trainControl(method = 'cv', number = 5,
                          classProbs = T, summaryFunction = twoClassSummary,
                          search = 'random') # 随机调参

for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed, method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed, method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed, method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed, method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed, method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
# -- ------------------------------------------------------------------------------------------------------------------
# e.模型预测及评价
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for ( a in 1:11) {
  Fit <- readRDS(filenames[a])
  predict.train <- predict(Fit, newdata = train.transformed)
  predict.test <- predict(Fit, newdata = test.transformed)
  # ROC分析的时候要区分方向，否则容易得到相反的AUC
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test) %>% mutate(method = myMethod[a])
  rownames(Fit.roc) <- myMethod[a]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for ( a in 1:11) {
  Fit <- readRDS(filenames[a])
  predict.train <- predict(Fit, newdata = train.transformed)
  predict.test <- predict(Fit, newdata = test.transformed)
  # ROC分析的时候要区分方向，否则容易得到相反的AUC
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test) %>% mutate(method = myMethod[a])
  rownames(Fit.roc) <- myMethod[a]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for ( a in 1:11) {
  Fit <- readRDS(filenames[a])
  predict.train <- predict(Fit, newdata = train.transformed)
  predict.test <- predict(Fit, newdata = test.transformed)
  # ROC分析的时候要区分方向，否则容易得到相反的AUC
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test) %>% mutate(method = myMethod[a])
  rownames(Fit.roc) <- myMethod[a]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for ( a in 1:11) {
  Fit <- readRDS(filenames[a])
  predict.train <- predict(Fit, newdata = train.transformed)
  predict.test <- predict(Fit, newdata = test.transformed)
  # ROC分析的时候要区分方向，否则容易得到相反的AUC
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test) %>% mutate(method = myMethod[a])
  rownames(Fit.roc) <- myMethod[a]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
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
  # ROC分析的时候要区分方向，否则容易得到相反的AUC
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test) %>% mutate(method = myMethod[a])
  rownames(Fit.roc) <- myMethod[a]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/HumanUrine/caret/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}

# -- ------------------------------------------------------------------------------------------------------------------
# HumanSaliva
# -- ------------------------------------------------------------------------------------------------------------------
# -- ------------------------------------------------------------------------------------------------------------------
# HumanPlasma
# -- ------------------------------------------------------------------------------------------------------------------
# -- ------------------------------------------------------------------------------------------------------------------
# HumanTear
# -- ------------------------------------------------------------------------------------------------------------------
# -- ------------------------------------------------------------------------------------------------------------------
# 由于之前已经跑过这部分的内容了，所以直接把之前的结果COPY过来
filnames <- Sys.glob('D:/2.AD.PBA.CAM.Article.V1.Proj/Results/*/*/*.AllSample.pseudo10.CPM.quantiles.train_*_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
file.copy(from = filnames, to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/Results/Figure5/Files/')
# 将这些AUC的值整理成一个EXCEL

# -- ------------------------------------------------------------------------------------------------------------------
# 4.将caret建模的AUC热图展示，从中挑选3次效果好的
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data1 <- read.csv('Results/Figure5/Files/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data2 <- read.csv('Results/Figure5/Files/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data3 <- read.csv('Results/Figure5/Files/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data <- left_join(data1, data2, by = 'method')
data <- left_join(data, data3, by = 'method')
rownames(data) <- data$method
data.Urine <- data %>% dplyr::select(method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data.Urine$Train.Mean <- rowMeans(data.Urine[2:4])
data.Urine$Train.sd <- apply(data.Urine[2:4], 1, sd)
data.Urine$Test.Mean <- rowMeans(data.Urine[5:7])
data.Urine$Test.sd <- apply(data.Urine[5:7], 1, sd)
data.Urine.heatmap <- data.Urine %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)

data1 <- read.csv('Results/Figure5/Files/HumanPlasma.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data2 <- read.csv('Results/Figure5/Files/HumanPlasma.AllSample.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data3 <- read.csv('Results/Figure5/Files/HumanPlasma.AllSample.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data <- left_join(data1, data2, by = 'method')
data <- left_join(data, data3, by = 'method')
rownames(data) <- data$method
data.Plasma <- data %>% dplyr::select(method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data.Plasma$Train.Mean <- rowMeans(data.Plasma[2:4])
data.Plasma$Train.sd <- apply(data.Plasma[2:4], 1, sd)
data.Plasma$Test.Mean <- rowMeans(data.Plasma[5:7])
data.Plasma$Test.sd <- apply(data.Plasma[5:7], 1, sd)
data.Plasma.heatmap <- data.Plasma %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)

data1 <- read.csv('Results/Figure5/Files/HumanSaliva.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data2 <- read.csv('Results/Figure5/Files/HumanSaliva.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data3 <- read.csv('Results/Figure5/Files/HumanSaliva.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data <- left_join(data1, data2, by = 'method')
data <- left_join(data, data3, by = 'method')
rownames(data) <- data$method
data.Saliva <- data %>% dplyr::select(method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data.Saliva$Train.Mean <- rowMeans(data.Saliva[2:4])
data.Saliva$Train.sd <- apply(data.Saliva[2:4], 1, sd)
data.Saliva$Test.Mean <- rowMeans(data.Saliva[5:7])
data.Saliva$Test.sd <- apply(data.Saliva[5:7], 1, sd)
data.Saliva.heatmap <- data.Saliva %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)

data1 <- read.csv('Results/Figure5/Files/HumanTear.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data2 <- read.csv('Results/Figure5/Files/HumanTear.AllSample.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data3 <- read.csv('Results/Figure5/Files/HumanTear.AllSample.pseudo10.CPM.quantiles.train_12345_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data <- left_join(data1, data2, by = 'method')
data <- left_join(data, data3, by = 'method')
rownames(data) <- data$method
data.Tear <- data %>% dplyr::select(method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data.Tear$Train.Mean <- rowMeans(data.Tear[2:4])
data.Tear$Train.sd <- apply(data.Tear[2:4], 1, sd)
data.Tear$Test.Mean <- rowMeans(data.Tear[5:7])
data.Tear$Test.sd <- apply(data.Tear[5:7], 1, sd)
data.Tear.heatmap <- data.Tear %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)
# 作图
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
         display_numbers = data.heatmap[5:6], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(11, 22, 33),
         filename = 'Results/Figure5/Human.pseudo10_caret.ROC 2.png') # 保存，自动调整纸张大小
pheatmap(data.heatmap[1:2], colorRampPalette(c('#90CAF8FF', '#FFF2DFFF', '#FF9800FF'))(12), 
         cellwidth = 40, cellheight = 10,
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(11, 22, 33),
         filename = 'Results/Figure5/Human.pseudo10_caret.ROC 2.pdf') # 保存，自动调整纸张大小
pheatmap(data.heatmap[1:2], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(11, 22, 33),
         filename = 'Results/Figure5/Human.pseudo10_caret.ROC.png') # 保存，自动调整纸张大小
pheatmap(data.heatmap[1:2], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(11, 22, 33),
         filename = 'Results/Figure5/Human.pseudo10_caret.ROC.pdf') # 保存，自动调整纸张大小

# -- ------------------------------------------------------------------------------------------------------------------
# 5.展示每种体液建模效果最好的ROC图和Feature排序
# -- ------------------------------------------------------------------------------------------------------------------
# 将最好的Rdata copy到这里的文件
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData',
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/Results/Figure5/Files/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/Results/HumanUrine/caret/RDS/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.nb.Fit.rds',
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/Results/Figure5/Files/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/Results/HumanPlasma/caret/RDS/HumanPlasma.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData',
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/Results/Figure5/Files/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/Results/HumanPlasma/caret/RDS/HumanPlasma.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.svmRadial.Fit.rds',
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/Results/Figure5/Files/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/Results/HumanSaliva/caret/RDS/HumanSaliva.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData',
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/Results/Figure5/Files/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/Results/HumanSaliva/caret/RDS/HumanSaliva.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.AdaBag.Fit.rds',
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/Results/Figure5/Files/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/Results/HumanTear/caret/RDS/HumanTear.AllSample.pseudo10.CPM.quantiles.train_12345_LASSO_Logistic.Pos_Top20_caret.RData',
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/Results/Figure5/Files/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/Results/HumanTear/caret/RDS/HumanTear.AllSample.pseudo10.CPM.quantiles.train_12345_LASSO_Logistic.Pos_Top20.nb.Fit.rds',
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/Results/Figure5/Files/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/Results/HumanTear/caret/RDS/HumanTear.AllSample.pseudo10.CPM.quantiles.train_12345_LASSO_Logistic.Pos_Top20.C5.0.Fit.rds',
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/Results/Figure5/Files/')
# -- ------------------------------------------------------------------------------------------------------------------
# HumanUrine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure5/Files/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
Fit <- readRDS('Results/Figure5/Files/HumanUrine.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.nb.Fit.rds')
var.importance <- varImp(Fit, scale = T)
png('Results/Figure5/HumanUrine.pseudo10.CPM.quantiles.nb.Fit.2022.VarImportance.png', width = 300, height = 400)
plot(var.importance, xlab = 'Importance', main = 'nb')
dev.off()
pdf('Results/Figure5/HumanUrine.pseudo10.CPM.quantiles.nb.Fit.2022.VarImportance.pdf', width = 3, height = 4)
plot(var.importance, xlab = 'Importance', main = 'nb')
dev.off()
# ROC图
predict.test <- predict(Fit, newdata = test.transformed, type = 'prob')
roc.test <- roc(test.transformed$Group, predict.test$AD)$auc
rocPlot <- roc(test.transformed$Group, predict.test$AD, levels = c('AD', 'Con'), direction = '>')
png('Results/Figure5/HumanUrine.pseudo10.CPM.quantiles.nb.Fit.ROC.png', width = 400, height = 400)
plot.roc(rocPlot, col="red", lwd = 4, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.5, auc.polygon = TRUE, auc.polygon.col = "#F0E685CC",
         print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
pdf('Results/Figure5/HumanUrine.pseudo10.CPM.quantiles.nb.Fit.ROC.pdf', width = 4, height = 4)
plot.roc(rocPlot, col="red", lwd = 4, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.5, auc.polygon = TRUE, auc.polygon.col = "#F0E685CC",
         print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
# TOP3蛋白Train的表达值图
for (i in 2:ncol(data.train)) {
  myTitle <- gsub("/", "_", colnames(data.train)[i])
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(data.train, aes(x = Group, fill = Group)) + theme_bw() +
    geom_boxplot(aes_(y = as.name(colnames(data.train)[i])), width = 0.6, size = 0.8) +  
    geom_violin(aes_(y = as.name(colnames(data.train)[i])), scale = 'width', width = 0.8, alpha = 0.6) +
    geom_jitter(aes_(y = as.name(colnames(data.train)[i])), width = 0.3, alpha = 0.6, size = 2) + 
    geom_signif(aes_(y = as.name(colnames(data.train)[i])), 
                comparisons = list(c('AD', 'Con')), 
                map_signif_level = F, step_increase = 0.01, 
                tip_length = 0.01, test = "wilcox.test", textsize = 2) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure5/HumanUrine_pseudo10.CPM.quantiles_2021.nb_Train_', myTitle, '.png')
  ggsave(p2, filename = fn2, width = 4, height = 4)
  fn2 <- str_c('Results/Figure5/HumanUrine_pseudo10.CPM.quantiles_2021.nb_Train_', myTitle, '.pdf')
  ggsave(p2, filename = fn2, width = 4, height = 4)
}

# -- ------------------------------------------------------------------------------------------------------------------
# HumanPlasma
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure5/Files/HumanPlasma.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
Fit <- readRDS('Results/Figure5/Files/HumanPlasma.AllSample.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.svmRadial.Fit.rds')
var.importance <- varImp(Fit, scale = T)
png('Results/Figure5/HumanPlasma.pseudo10.CPM.quantiles.svmRadial.Fit.2022.VarImportance.png', width = 300, height = 400)
plot(var.importance, xlab = 'Importance', main = 'svmRadial')
dev.off()
pdf('Results/Figure5/HumanPlasma.pseudo10.CPM.quantiles.svmRadial.Fit.2022.VarImportance.pdf', width = 3, height = 4)
plot(var.importance, xlab = 'Importance', main = 'svmRadial')
dev.off()
# ROC图
predict.test <- predict(Fit, newdata = test.transformed, type = 'prob')
roc.test <- roc(test.transformed$Group, predict.test$AD)$auc
rocPlot <- roc(test.transformed$Group, predict.test$AD, levels = c('AD', 'Con'), direction = '>')
png('Results/Figure5/HumanPlasma.pseudo10.CPM.quantiles.svmRadial.Fit.ROC.png', width = 400, height = 400)
plot.roc(rocPlot, col="red", lwd = 4, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.5, auc.polygon = TRUE, auc.polygon.col = "#D595A7CC",
         print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
pdf('Results/Figure5/HumanPlasma.pseudo10.CPM.quantiles.svmRadial.Fit.ROC.pdf', width = 4, height = 4)
plot.roc(rocPlot, col="red", lwd = 4, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.5, auc.polygon = TRUE, auc.polygon.col = "#D595A7CC",
         print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
# TOP3蛋白Train的表达值图
for (i in 2:ncol(data.train)) {
  myTitle <- gsub("/", "_", colnames(data.train)[i])
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(data.train, aes(x = Group, fill = Group)) + theme_bw() +
    geom_boxplot(aes_(y = as.name(colnames(data.train)[i])), width = 0.6, size = 0.8) +  
    geom_violin(aes_(y = as.name(colnames(data.train)[i])), scale = 'width', width = 0.8, alpha = 0.6) +
    geom_jitter(aes_(y = as.name(colnames(data.train)[i])), width = 0.3, alpha = 0.6, size = 2) + 
    geom_signif(aes_(y = as.name(colnames(data.train)[i])), 
                comparisons = list(c('AD', 'Con')), 
                map_signif_level = F, step_increase = 0.01, 
                tip_length = 0.01, test = "wilcox.test", textsize = 2) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure5/HumanPlasma_pseudo10.CPM.quantiles_1234.svmRadial_Train_', myTitle, '.png')
  ggsave(p2, filename = fn2, width = 4, height = 4)
  fn2 <- str_c('Results/Figure5/HumanPlasma_pseudo10.CPM.quantiles_1234.svmRadial_Train_', myTitle, '.pdf')
  ggsave(p2, filename = fn2, width = 4, height = 4)
}

# -- ------------------------------------------------------------------------------------------------------------------
# HumanSaliva
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure5/Files/HumanSaliva.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
Fit <- readRDS('Results/Figure5/Files/HumanSaliva.AllSample.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.AdaBag.Fit.rds')
var.importance <- varImp(Fit, scale = T)
png('Results/Figure5/HumanSaliva.pseudo10.CPM.quantiles.AdaBag.Fit.2022.VarImportance.png', width = 300, height = 400)
plot(var.importance, xlab = 'Importance', main = 'AdaBag')
dev.off()
pdf('Results/Figure5/HumanSaliva.pseudo10.CPM.quantiles.AdaBag.Fit.2022.VarImportance.pdf', width = 3, height = 4)
plot(var.importance, xlab = 'Importance', main = 'AdaBag')
dev.off()
# ROC图
predict.test <- predict(Fit, newdata = test.transformed, type = 'prob')
roc.test <- roc(test.transformed$Group, predict.test$AD)$auc
rocPlot <- roc(test.transformed$Group, predict.test$AD, levels = c('AD', 'Con'), direction = '>')
png('Results/Figure5/HumanSaliva.pseudo10.CPM.quantiles.AdaBag.Fit.ROC.png', width = 400, height = 400)
plot.roc(rocPlot, col="red", lwd = 4, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.5, auc.polygon = TRUE, auc.polygon.col = "#925E9FAA",
         print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
pdf('Results/Figure5/HumanSaliva.pseudo10.CPM.quantiles.AdaBag.Fit.ROC.pdf', width = 4, height = 4)
plot.roc(rocPlot, col="red", lwd = 4, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.5, auc.polygon = TRUE, auc.polygon.col = "#925E9FAA",
         print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
# TOP3蛋白Train的表达值图
for (i in 2:ncol(data.train)) {
  myTitle <- gsub("/", "_", colnames(data.train)[i])
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(data.train, aes(x = Group, fill = Group)) + theme_bw() +
    geom_boxplot(aes_(y = as.name(colnames(data.train)[i])), width = 0.6, size = 0.8) +  
    geom_violin(aes_(y = as.name(colnames(data.train)[i])), scale = 'width', width = 0.8, alpha = 0.6) +
    geom_jitter(aes_(y = as.name(colnames(data.train)[i])), width = 0.3, alpha = 0.6, size = 2) + 
    geom_signif(aes_(y = as.name(colnames(data.train)[i])), 
                comparisons = list(c('AD', 'Con')), 
                map_signif_level = F, step_increase = 0.01, 
                tip_length = 0.01, test = "wilcox.test", textsize = 2) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure5/HumanSaliva_pseudo10.CPM.quantiles_2021.AdaBag_Train_', myTitle, '.png')
  ggsave(p2, filename = fn2, width = 4, height = 4)
  fn2 <- str_c('Results/Figure5/HumanSaliva_pseudo10.CPM.quantiles_2021.AdaBag_Train_', myTitle, '.pdf')
  ggsave(p2, filename = fn2, width = 4, height = 4)
}

# -- ------------------------------------------------------------------------------------------------------------------
# HumanTear
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure5/Files/HumanTear.AllSample.pseudo10.CPM.quantiles.train_12345_LASSO_Logistic.Pos_Top20_caret.RData')
Fit <- readRDS('Results/Figure5/Files/HumanTear.AllSample.pseudo10.CPM.quantiles.train_12345_LASSO_Logistic.Pos_Top20.nb.Fit.rds')
var.importance <- varImp(Fit, scale = T)
png('Results/Figure5/HumanTear.pseudo10.CPM.quantiles.nb.Fit.2022.VarImportance.png', width = 300, height = 400)
plot(var.importance, xlab = 'Importance', main = 'nb')
dev.off()
pdf('Results/Figure5/HumanTear.pseudo10.CPM.quantiles.nb.Fit.2022.VarImportance.pdf', width = 3, height = 4)
plot(var.importance, xlab = 'Importance', main = 'nb')
dev.off()
# ROC图
Fit <- readRDS('Results/Figure5/Files/HumanTear.AllSample.pseudo10.CPM.quantiles.train_12345_LASSO_Logistic.Pos_Top20.C5.0.Fit.rds')
predict.test <- predict(Fit, newdata = test.transformed, type = 'prob')
roc.test <- roc(test.transformed$Group, predict.test$AD)$auc
rocPlot <- roc(test.transformed$Group, predict.test$AD, levels = c('AD', 'Con'), direction = '>')
png('Results/Figure5/HumanTear.pseudo10.CPM.quantiles.nb.Fit.ROC.png', width = 400, height = 400)
plot.roc(rocPlot, col="red", lwd = 4, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.5, auc.polygon = TRUE, auc.polygon.col = "#6BD76BCC",
         print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
pdf('Results/Figure5/HumanTear.pseudo10.CPM.quantiles.nb.Fit.ROC.pdf', width = 4, height = 4)
plot.roc(rocPlot, col="red", lwd = 4, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.5, auc.polygon = TRUE, auc.polygon.col = "#6BD76BCC",
         print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
# TOP3蛋白Train的表达值图
for (i in 2:ncol(data.train)) {
  myTitle <- gsub("/", "_", colnames(data.train)[i])
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(data.train, aes(x = Group, fill = Group)) + theme_bw() +
    geom_boxplot(aes_(y = as.name(colnames(data.train)[i])), width = 0.6, size = 0.8) +  
    geom_violin(aes_(y = as.name(colnames(data.train)[i])), scale = 'width', width = 0.8, alpha = 0.6) +
    geom_jitter(aes_(y = as.name(colnames(data.train)[i])), width = 0.3, alpha = 0.6, size = 2) + 
    geom_signif(aes_(y = as.name(colnames(data.train)[i])), 
                comparisons = list(c('AD', 'Con')), 
                map_signif_level = F, step_increase = 0.01, 
                tip_length = 0.01, test = "wilcox.test", textsize = 2) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure5/HumanTear_pseudo10.CPM.quantiles_12345.nb_Train_', myTitle, '.png')
  ggsave(p2, filename = fn2, width = 4, height = 4)
  fn2 <- str_c('Results/Figure5/HumanTear_pseudo10.CPM.quantiles_12345.nb_Train_', myTitle, '.pdf')
  ggsave(p2, filename = fn2, width = 4, height = 4)
}

















