# -- ------------------------------------------------------------------------------------------------------------------
getwd()
setwd("D:/2.AD.PBA.CAM.Article.V3.Proj")
# 设置最大的内存
if(.Platform$OS.type == "windows") withAutoprint({
  memory.size()
  memory.size(TRUE)
  memory.limit()
})
# options(future.globals.maxSize = 30*1024^3) #设置为30G

# -- ------------------------------------------------------------------------------------------------------------------
# 1.对小鼠的样本原始数据进行pseudo处理:1-10,1-20
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
sampleID <- read.csv('Mouse.PT_SampleID_SampleName.csv')
head(sampleID)
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('RawData/Mouse/',sampleID[i,1],'.total_ev_protein.csv')
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
  fn1 <- str_c('RawData/Mouse/',sampleID[i,1],'.total_ev_protein.csv')
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
sampleID <- read.csv('Mouse.PT_SampleID_SampleName.csv')
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
data <- read_csv('Results/Mouse_Raw.Counts.csv', show_col_types = F)
names(data)[1] <- 'Pro'
data.anno <- read_csv('Mouse.PT_SampleID_SampleName.csv', show_col_types = F)
psedoSize <- c(10,20)
for (i in 1:nrow(data.anno)) {
  for (j in psedoSize) {
    dt <- data %>% dplyr::select(Pro, data.anno$SampleName[i])
    fn1 <- str_c('pseudoData/pseudo', j, '/', data.anno$SampleID[i], '.pseudo', j ,'.csv')
    temp <- read.csv(fn1, row.names = 1)
    temp <- t(temp) %>% as.data.frame()
    temp$Pro <- rownames(temp)
    dt <- inner_join(dt, temp, by = 'Pro')
    dt.cor <- cor(dt[-1])
    fn2 <- str_c('Results/FigureS4/', data.anno$SampleName[i], '_Raw.pseudo', j, '.png')
    pheatmap(dt.cor, color = pal_material(palette = c("orange"), n = 20)(20), cellwidth = 10, cellheight = 6,
             display_numbers = T, number_color = 'black', number_format = '%.3f', 
             fontsize_number = 3, fontsize_row = 6, fontsize_col = 8,
             cluster_rows = F, cluster_cols = F,
             filename = fn2) # 保存，自动调整纸张大小
    fn2 <- str_c('Results/FigureS4/', data.anno$SampleName[i], '_Raw.pseudo', j, '.pdf')
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
library(tidyverse) # 处理数据神包
library(glmnet) # 进行正则技术所需
library(edgeR) # 用于cpm标准化
library(preprocessCore) # normalize.quantiles函数进行标准化
library(caret)
library(pROC)
library(gbm)
# -- ------------------------------------------------------------------------------------------------------------------
# MouseUrine
# -- ------------------------------------------------------------------------------------------------------------------
# 读取pseudo10文件合成总的表达矩阵
rm(list=ls())
gc()
data.anno <- read_csv('Results/MouseUrine_Samples.Annotation.csv', show_col_types = F)
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
write_csv(dt, 'pseudoData/MouseUrine.AllSample.pseudo10.csv')
# 对数据进行标准化，并进行样本分组注释
rm(list=ls())
gc()
data <- read.csv('pseudoData/MouseUrine.AllSample.pseudo10.csv', row.names = 1)
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
data$Class <- ifelse(str_detect(data$ID, '15'), 15,
                        ifelse(str_detect(data$ID, '10'), 10, 3))
data <- relocate(data, ID, Group, Class)
write.csv(data, 'pseudoData/MouseUrine.AllSample.pseudo10.CPM.quantiles.Group.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# 用caret包进行机器学习
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('pseudoData/MouseUrine.AllSample.pseudo10.CPM.quantiles.Group.csv', row.names = 1)
data$ID <- data$ID %>% factor()
data$Group <- data$Group %>% factor()
data$Class <- data$Class %>% factor(levels = c(3,10,15))
# -- ------------------------------------------------------------------------------------------------------------------
# a.划分训练集与测试集:
# -- ------------------------------------------------------------------------------------------------------------------
data.ID <- levels(data$ID) # 只有factor才能用levels这个函数
data.ID <- data.frame(ID = data.ID)
data.ID$Group <- ifelse(str_detect(data.ID$ID, 'D'), 'AD', 'Con')
data.ID$Class <- ifelse(str_detect(data.ID$ID, '15'), 15,
                     ifelse(str_detect(data.ID$ID, '10'), 10, 3))
data.ID$Class <- factor(data.ID$Class, levels = c(3,10,15), labels = c('3M', '10M', '15M'))
data.ID$Split <- str_c(data.ID$Class, data.ID$Group, sep = '.')
data.ID$Split <- factor(data.ID$Split, levels = c('3M.AD', '3M.Con', '10M.AD', '10M.Con', '15M.AD', '15M.Con'))
set.seed(1234)
set.seed(5678)
set.seed(2020)
set.seed(2021)
set.seed(2022)
train.ID <- createDataPartition(data.ID$Split, p = 0.5, list = FALSE)
train.ID <- data.ID$ID[train.ID]
# load('D:/2.AD.PBA.CAM.Article.V2.Proj/Results/MouseUrine/caret_data.ID_1/RDS/MouseUrine.PT.pseudo10.Lasso.CPM.quantiles.caret.RData')
# ID <- rownames(data.train) %>% str_split('_s') %>% as.data.frame() %>% t() %>% as.data.frame()
# train.ID <- ID$V1 %>% factor %>% levels()
data.train <- data %>% filter(ID %in% train.ID) %>% dplyr::select(!ID)
data.test <- data %>% filter(!ID %in% train.ID) %>% dplyr::select(!ID)
# -- ------------------------------------------------------------------------------------------------------------------
# b.用LASSO筛选变量
# -- ------------------------------------------------------------------------------------------------------------------
X <- as.matrix(data.train[3:ncol(data.train)])
Y <- as.matrix(data.train[1])
# glmnet()语法中alpha=0表示岭回归，1表示LASSO回归
myLasso <- glmnet(X, Y, alpha = 1, family = 'binomial', nlambda = 200) # glmnet默认运行100次
NROW(myLasso$lambda)
min(myLasso$lambda) # 选最优模型的lambda值
# 绘制图形
pdf('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_5678_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2020_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2021_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2022_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO.pdf', width = 6, height = 8)
plot(myLasso, xvar = 'lambda', label = T)
dev.off()
lasso.coef <- coef(myLasso, s = min(myLasso$lambda)) # 最优解下的回归系数
# 只有筛选出的自变量才有回归系数
# 将变量根据回归系数进行排序，回归系数越大的变量越重要
my.lasso.coef <- as.matrix(lasso.coef) %>% as.data.frame()
my.lasso.coef$abs <- abs(my.lasso.coef$s1)
my.lasso.coef <- arrange(my.lasso.coef, desc(abs)) %>% filter(abs > 0)
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_5678_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2020_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2021_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2022_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO.csv')
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
write.csv(res1, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.csv')
res2 <- res1 %>% filter(`Pr(>|z|)` < 0.05)
write.csv(res2, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# c.用阳性变量构建Train和Test数据集
# -- ------------------------------------------------------------------------------------------------------------------
SigPro <- read.csv('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos.csv', row.names = 1)
myPro <- str_replace_all(rownames(SigPro)[1:20], '/', '.') %>% str_replace_all('-', '.')
myPro <- intersect(myPro, names(data.train))
data.train <- data.train %>% dplyr::select(Class, Group, all_of(myPro))
data.test <- data.test %>% dplyr::select(Class, Group, all_of(myPro))
# 将数据进行标准化并补足缺失值，这时可以用preProcess命令，缺省参数是标准化数据。?preProcess
preProcValues <- preProcess(data.train, method = c('center','scale')) # 用于进行其他验证数据集的标准化
train.transformed <- predict(preProcValues, data.train)
test.transformed <- predict(preProcValues, data.test)
str(train.transformed)
# 保存数据
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20_caret.RData')
# -- ------------------------------------------------------------------------------------------------------------------
# d.模型训练及调参
# -- ------------------------------------------------------------------------------------------------------------------
# 常用的分类算法包括：NBC（Naive Bayesian Classifier，朴素贝叶斯分类）算法、
# LR（Logistic Regress，逻辑回归）算法、C5.0 决策树算法、SVM（Support Vector Machine，支持向量机）算法、
# KNN(K-Nearest Neighbor，K 最近邻近)算法、ANN（Artificial Neural Network，人工神经网络）算法
# Stochastic Gradient Boosting(GBM)
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20_caret.RData')

myMethod <- c('glmnet','kknn','C5.0','rf','AdaBag','gbm','nnet','svmLinear','svmPoly','svmRadial','nb')
fitControl = trainControl(method = 'cv', number = 5,
                          classProbs = T, summaryFunction = twoClassSummary,
                          search = 'random') # 随机调参

for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}

# -- ------------------------------------------------------------------------------------------------------------------
# e.模型预测及评价
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}

# -- ------------------------------------------------------------------------------------------------------------------
# MouseNIS
# -- ------------------------------------------------------------------------------------------------------------------
# 读取pseudo10文件合成总的表达矩阵
rm(list=ls())
gc()
data.anno <- read_csv('Results/MouseNIS_Samples.Annotation.csv', show_col_types = F)
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
write_csv(dt, 'pseudoData/MouseNIS.AllSample.pseudo10.csv')
# 对数据进行标准化，并进行样本分组注释
rm(list=ls())
gc()
data <- read.csv('pseudoData/MouseNIS.AllSample.pseudo10.csv', row.names = 1)
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
data$Class <- ifelse(str_detect(data$ID, '15'), 15,
                     ifelse(str_detect(data$ID, '10'), 10, 3))
data <- relocate(data, ID, Group, Class)
write.csv(data, 'pseudoData/MouseNIS.AllSample.pseudo10.CPM.quantiles.Group.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# 用caret包进行机器学习
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('pseudoData/MouseNIS.AllSample.pseudo10.CPM.quantiles.Group.csv', row.names = 1)
data$ID <- data$ID %>% factor()
data$Group <- data$Group %>% factor()
data$Class <- data$Class %>% factor(levels = c(3,10,15))
# -- ------------------------------------------------------------------------------------------------------------------
# a.划分训练集与测试集:
# -- ------------------------------------------------------------------------------------------------------------------
data.ID <- levels(data$ID) # 只有factor才能用levels这个函数
data.ID <- data.frame(ID = data.ID)
data.ID$Group <- ifelse(str_detect(data.ID$ID, 'D'), 'AD', 'Con')
data.ID$Class <- ifelse(str_detect(data.ID$ID, '15'), 15,
                        ifelse(str_detect(data.ID$ID, '10'), 10, 3))
data.ID$Class <- factor(data.ID$Class, levels = c(3,10,15), labels = c('3M', '10M', '15M'))
data.ID$Split <- str_c(data.ID$Class, data.ID$Group, sep = '.')
data.ID$Split <- factor(data.ID$Split, levels = c('3M.AD', '3M.Con', '10M.AD', '10M.Con', '15M.AD', '15M.Con'))
set.seed(1234)
set.seed(5678)
set.seed(2020)
set.seed(2021)
set.seed(2022)
train.ID <- createDataPartition(data.ID$Split, p = 0.5, list = FALSE)
train.ID <- data.ID$ID[train.ID]
# load('D:/2.AD.PBA.CAM.Article.V2.Proj/Results/MouseNIS/caret_data.ID_2/RDS/MouseNIS.PT.pseudo10.Lasso.CPM.quantiles.caret.RData')
# ID <- rownames(data.train) %>% str_split('_s') %>% as.data.frame() %>% t() %>% as.data.frame()
# train.ID <- ID$V1 %>% factor %>% levels()
data.train <- data %>% filter(ID %in% train.ID) %>% dplyr::select(!ID)
data.test <- data %>% filter(!ID %in% train.ID) %>% dplyr::select(!ID)
# -- ------------------------------------------------------------------------------------------------------------------
# b.用LASSO筛选变量
# -- ------------------------------------------------------------------------------------------------------------------
X <- as.matrix(data.train[3:ncol(data.train)])
Y <- as.matrix(data.train[1])
# glmnet()语法中alpha=0表示岭回归，1表示LASSO回归
myLasso <- glmnet(X, Y, alpha = 1, family = 'binomial', nlambda = 200) # glmnet默认运行100次
NROW(myLasso$lambda)
min(myLasso$lambda) # 选最优模型的lambda值
# 绘制图形
pdf('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_1234_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_5678_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2020_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2022_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_old_LASSO.pdf', width = 6, height = 8)
plot(myLasso, xvar = 'lambda', label = T)
dev.off()
lasso.coef <- coef(myLasso, s = min(myLasso$lambda)) # 最优解下的回归系数
# 只有筛选出的自变量才有回归系数
# 将变量根据回归系数进行排序，回归系数越大的变量越重要
my.lasso.coef <- as.matrix(lasso.coef) %>% as.data.frame()
my.lasso.coef$abs <- abs(my.lasso.coef$s1)
my.lasso.coef <- arrange(my.lasso.coef, desc(abs)) %>% filter(abs > 0)
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_1234_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_5678_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2020_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2022_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_old_LASSO.csv')
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
write.csv(res1, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.csv')
res2 <- res1 %>% filter(`Pr(>|z|)` < 0.05)
write.csv(res2, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# c.用阳性变量构建Train和Test数据集
# -- ------------------------------------------------------------------------------------------------------------------
SigPro <- read.csv('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos.csv', row.names = 1)
myPro <- str_replace_all(rownames(SigPro)[1:20], '/', '.') %>% str_replace_all('-', '.')
myPro <- intersect(myPro, names(data.train))
data.train <- data.train %>% dplyr::select(Class, Group, all_of(myPro))
data.test <- data.test %>% dplyr::select(Class, Group, all_of(myPro))
# 将数据进行标准化并补足缺失值，这时可以用preProcess命令，缺省参数是标准化数据。?preProcess
preProcValues <- preProcess(data.train, method = c('center','scale')) # 用于进行其他验证数据集的标准化
train.transformed <- predict(preProcValues, data.train)
test.transformed <- predict(preProcValues, data.test)
str(train.transformed)
# 保存数据
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20_caret.RData')
# -- ------------------------------------------------------------------------------------------------------------------
# d.模型训练及调参
# -- ------------------------------------------------------------------------------------------------------------------
# 常用的分类算法包括：NBC（Naive Bayesian Classifier，朴素贝叶斯分类）算法、
# LR（Logistic Regress，逻辑回归）算法、C5.0 决策树算法、SVM（Support Vector Machine，支持向量机）算法、
# KNN(K-Nearest Neighbor，K 最近邻近)算法、ANN（Artificial Neural Network，人工神经网络）算法
# Stochastic Gradient Boosting(GBM)
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20_caret.RData')

myMethod <- c('glmnet','kknn','C5.0','rf','AdaBag','gbm','nnet','svmLinear','svmPoly','svmRadial','nb')
fitControl = trainControl(method = 'cv', number = 5,
                          classProbs = T, summaryFunction = twoClassSummary,
                          search = 'random') # 随机调参

for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
# -- ------------------------------------------------------------------------------------------------------------------
# e.模型预测及评价
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}


# -- ------------------------------------------------------------------------------------------------------------------
# MouseSerum
# -- ------------------------------------------------------------------------------------------------------------------
# 读取pseudo10文件合成总的表达矩阵
rm(list=ls())
gc()
data.anno <- read_csv('Results/MouseSerum_Samples.Annotation.csv', show_col_types = F)
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
write_csv(dt, 'pseudoData/MouseSerum.AllSample.pseudo10.csv')
# 对数据进行标准化，并进行样本分组注释
rm(list=ls())
gc()
data <- read.csv('pseudoData/MouseSerum.AllSample.pseudo10.csv', row.names = 1)
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
data$Class <- ifelse(str_detect(data$ID, '15'), 15,
                     ifelse(str_detect(data$ID, '10'), 10, 3))
data <- relocate(data, ID, Group, Class)
write.csv(data, 'pseudoData/MouseSerum.AllSample.pseudo10.CPM.quantiles.Group.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# 用caret包进行机器学习
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('pseudoData/MouseSerum.AllSample.pseudo10.CPM.quantiles.Group.csv', row.names = 1)
data$ID <- data$ID %>% factor()
data$Group <- data$Group %>% factor()
data$Class <- data$Class %>% factor(levels = c(3,10,15))
# -- ------------------------------------------------------------------------------------------------------------------
# a.划分训练集与测试集:
# -- ------------------------------------------------------------------------------------------------------------------
data.ID <- levels(data$ID) # 只有factor才能用levels这个函数
data.ID <- data.frame(ID = data.ID)
data.ID$Group <- ifelse(str_detect(data.ID$ID, 'D'), 'AD', 'Con')
data.ID$Class <- ifelse(str_detect(data.ID$ID, '15'), 15,
                        ifelse(str_detect(data.ID$ID, '10'), 10, 3))
data.ID$Class <- factor(data.ID$Class, levels = c(3,10,15), labels = c('3M', '10M', '15M'))
data.ID$Split <- str_c(data.ID$Class, data.ID$Group, sep = '.')
data.ID$Split <- factor(data.ID$Split, levels = c('3M.AD', '3M.Con', '10M.AD', '10M.Con', '15M.AD', '15M.Con'))
set.seed(1234)
set.seed(5678)
set.seed(2020)
set.seed(2021)
set.seed(2022)
train.ID <- createDataPartition(data.ID$Split, p = 0.5, list = FALSE)
train.ID <- data.ID$ID[train.ID]
data.train <- data %>% filter(ID %in% train.ID) %>% dplyr::select(!ID)
data.test <- data %>% filter(!ID %in% train.ID) %>% dplyr::select(!ID)
# -- ------------------------------------------------------------------------------------------------------------------
# b.用LASSO筛选变量
# -- ------------------------------------------------------------------------------------------------------------------
X <- as.matrix(data.train[3:ncol(data.train)])
Y <- as.matrix(data.train[1])
# glmnet()语法中alpha=0表示岭回归，1表示LASSO回归
myLasso <- glmnet(X, Y, alpha = 1, family = 'binomial', nlambda = 200) # glmnet默认运行100次
NROW(myLasso$lambda)
min(myLasso$lambda) # 选最优模型的lambda值
# 绘制图形
pdf('Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_5678_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2020_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2021_LASSO.pdf', width = 6, height = 8)
pdf('Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2022_LASSO.pdf', width = 6, height = 8)
plot(myLasso, xvar = 'lambda', label = T)
dev.off()
lasso.coef <- coef(myLasso, s = min(myLasso$lambda)) # 最优解下的回归系数
# 只有筛选出的自变量才有回归系数
# 将变量根据回归系数进行排序，回归系数越大的变量越重要
my.lasso.coef <- as.matrix(lasso.coef) %>% as.data.frame()
my.lasso.coef$abs <- abs(my.lasso.coef$s1)
my.lasso.coef <- arrange(my.lasso.coef, desc(abs)) %>% filter(abs > 0)
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_5678_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2020_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2021_LASSO.csv')
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2022_LASSO.csv')
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
write.csv(res1, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.csv')
write.csv(res1, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.csv')
res2 <- res1 %>% filter(`Pr(>|z|)` < 0.05)
write.csv(res2, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos.csv')
write.csv(res2, file = 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# c.用阳性变量构建Train和Test数据集
# -- ------------------------------------------------------------------------------------------------------------------
SigPro <- read.csv('Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos.csv', row.names = 1)
SigPro <- read.csv('Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos.csv', row.names = 1)
myPro <- str_replace_all(rownames(SigPro)[1:20], '/', '.') %>% str_replace_all('-', '.')
myPro <- intersect(myPro, names(data.train))
data.train <- data.train %>% dplyr::select(Class, Group, all_of(myPro))
data.test <- data.test %>% dplyr::select(Class, Group, all_of(myPro))
# 将数据进行标准化并补足缺失值，这时可以用preProcess命令，缺省参数是标准化数据。?preProcess
preProcValues <- preProcess(data.train, method = c('center','scale')) # 用于进行其他验证数据集的标准化
train.transformed <- predict(preProcValues, data.train)
test.transformed <- predict(preProcValues, data.test)
str(train.transformed)
# 保存数据
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')
# -- ------------------------------------------------------------------------------------------------------------------
# d.模型训练及调参
# -- ------------------------------------------------------------------------------------------------------------------
# 常用的分类算法包括：NBC（Naive Bayesian Classifier，朴素贝叶斯分类）算法、
# LR（Logistic Regress，逻辑回归）算法、C5.0 决策树算法、SVM（Support Vector Machine，支持向量机）算法、
# KNN(K-Nearest Neighbor，K 最近邻近)算法、ANN（Artificial Neural Network，人工神经网络）算法
# Stochastic Gradient Boosting(GBM)
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
load('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')

myMethod <- c('glmnet','kknn','C5.0','rf','AdaBag','gbm','nnet','svmLinear','svmPoly','svmRadial','nb')
fitControl = trainControl(method = 'cv', number = 5,
                          classProbs = T, summaryFunction = twoClassSummary,
                          search = 'random') # 随机调参

for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
               trControl = fitControl, 
               tuneLength = 30, # 随机调参的随机次数
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
# -- ------------------------------------------------------------------------------------------------------------------
# e.模型预测及评价
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.RData')
filenames <- Sys.glob('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20.*.Fit.rds')
myMethod <- str_remove(filenames, 'Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20.') %>%
  str_remove('.Fit.rds')
res <- data.frame()
for (j in 1:11) {
  Fit <- readRDS(filenames[j])
  predict.train <- predict(Fit, newdata = train.transformed[-1])
  predict.test <- predict(Fit, newdata = test.transformed[-1])
  test.transformed.3M <- test.transformed %>% filter(Class == 3)
  predict.test.3M <- predict(Fit, newdata = test.transformed.3M[-1])
  test.transformed.10M <- test.transformed %>% filter(Class == 10)
  predict.test.10M <- predict(Fit, newdata = test.transformed.10M[-1])
  test.transformed.15M <- test.transformed %>% filter(Class == 15)
  predict.test.15M <- predict(Fit, newdata = test.transformed.15M[-1])
  roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
  roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
  roc.test.3M <- roc(test.transformed.3M$Group, as.numeric(predict.test.3M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.10M <- roc(test.transformed.10M$Group, as.numeric(predict.test.10M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'), direction = '>')$auc
  roc.test.15M.2 <- roc(test.transformed.15M$Group, as.numeric(predict.test.15M), levels = c('Con','AD'))$auc
  Fit.roc <- data.frame(Train = roc.train, Test = roc.test, Test.3M = roc.test.3M, Test.10M = roc.test.10M,
                        Test.15M = roc.test.15M, Test.15M.2 = roc.test.15M.2) %>% mutate(method = myMethod[j])
  rownames(Fit.roc) <- myMethod[j]
  res <- rbind(Fit.roc, res)
  write.csv(res, 'Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.ROC.csv')
}

# -- ------------------------------------------------------------------------------------------------------------------
# 4.将caret建模的AUC热图展示，从中挑选3次效果好的
# -- ------------------------------------------------------------------------------------------------------------------
library(pheatmap)
rm(list=ls())
gc()
data1 <- read.csv('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_2020_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data2 <- read.csv('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data3 <- read.csv('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_5678_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data <- left_join(data1, data2, by = 'method')
data <- left_join(data, data3, by = 'method')
rownames(data) <- data$method
data.Urine <- data %>% dplyr::select(method, Train, Train.x, Train.y, Test, Test.x, Test.y,
                                     Test.3M, Test.3M.x, Test.3M.y, Test.10M, Test.10M.x, Test.10M.y,
                                     Test.15M, Test.15M.x, Test.15M.y)
data.Urine$Train.Mean <- rowMeans(data.Urine[2:4])
data.Urine$Train.sd <- apply(data.Urine[2:4], 1, sd)
data.Urine$Test.Mean <- rowMeans(data.Urine[5:7])
data.Urine$Test.sd <- apply(data.Urine[5:7], 1, sd)
data.Urine$Test.3M.Mean <- rowMeans(data.Urine[8:10])
data.Urine$Test.3M.sd <- apply(data.Urine[8:10], 1, sd)
data.Urine$Test.10M.Mean <- rowMeans(data.Urine[11:13])
data.Urine$Test.10M.sd <- apply(data.Urine[11:13], 1, sd)
data.Urine$Test.15M.Mean <- rowMeans(data.Urine[14:16])
data.Urine$Test.15M.sd <- apply(data.Urine[14:16], 1, sd)
data.Urine.heatmap <- data.Urine %>% dplyr::select(Train.Mean, Test.Mean, Test.3M.Mean, Test.10M.Mean, Test.15M.Mean,
                                                   Train.sd, Test.sd, Test.3M.sd, Test.10M.sd, Test.15M.sd)

data1 <- read.csv('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data2 <- read.csv('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data3 <- read.csv('Results/Figure3/Files/MouseNIS.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data <- left_join(data1, data2, by = 'method')
data <- left_join(data, data3, by = 'method')
rownames(data) <- data$method
data.NIS <- data %>% dplyr::select(method, Train, Train.x, Train.y, Test, Test.x, Test.y,
                                     Test.3M, Test.3M.x, Test.3M.y, Test.10M, Test.10M.x, Test.10M.y,
                                     Test.15M, Test.15M.x, Test.15M.y)
data.NIS$Train.Mean <- rowMeans(data.NIS[2:4])
data.NIS$Train.sd <- apply(data.NIS[2:4], 1, sd)
data.NIS$Test.Mean <- rowMeans(data.NIS[5:7])
data.NIS$Test.sd <- apply(data.NIS[5:7], 1, sd)
data.NIS$Test.3M.Mean <- rowMeans(data.NIS[8:10])
data.NIS$Test.3M.sd <- apply(data.NIS[8:10], 1, sd)
data.NIS$Test.10M.Mean <- rowMeans(data.NIS[11:13])
data.NIS$Test.10M.sd <- apply(data.NIS[11:13], 1, sd)
data.NIS$Test.15M.Mean <- rowMeans(data.NIS[14:16])
data.NIS$Test.15M.sd <- apply(data.NIS[14:16], 1, sd)
data.NIS.heatmap <- data.NIS %>% dplyr::select(Train.Mean, Test.Mean, Test.3M.Mean, Test.10M.Mean, Test.15M.Mean,
                                                   Train.sd, Test.sd, Test.3M.sd, Test.10M.sd, Test.15M.sd)

data1 <- read.csv('Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data2 <- read.csv('Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data3 <- read.csv('Results/Figure3/Files/MouseSerum.pseudo10.CPM.quantiles.train_2022_LASSO_Logistic.Pos_Top20_caret.ROC.csv', row.names = 1)
data <- left_join(data1, data2, by = 'method')
data <- left_join(data, data3, by = 'method')
rownames(data) <- data$method
data.Serum <- data %>% dplyr::select(method, Train, Train.x, Train.y, Test, Test.x, Test.y,
                                     Test.3M, Test.3M.x, Test.3M.y, Test.10M, Test.10M.x, Test.10M.y,
                                     Test.15M, Test.15M.x, Test.15M.y)
data.Serum$Train.Mean <- rowMeans(data.Serum[2:4])
data.Serum$Train.sd <- apply(data.Serum[2:4], 1, sd)
data.Serum$Test.Mean <- rowMeans(data.Serum[5:7])
data.Serum$Test.sd <- apply(data.Serum[5:7], 1, sd)
data.Serum$Test.3M.Mean <- rowMeans(data.Serum[8:10])
data.Serum$Test.3M.sd <- apply(data.Serum[8:10], 1, sd)
data.Serum$Test.10M.Mean <- rowMeans(data.Serum[11:13])
data.Serum$Test.10M.sd <- apply(data.Serum[11:13], 1, sd)
data.Serum$Test.15M.Mean <- rowMeans(data.Serum[14:16])
data.Serum$Test.15M.sd <- apply(data.Serum[14:16], 1, sd)
data.Serum.heatmap <- data.Serum %>% dplyr::select(Train.Mean, Test.Mean, Test.3M.Mean, Test.10M.Mean, Test.15M.Mean,
                                                   Train.sd, Test.sd, Test.3M.sd, Test.10M.sd, Test.15M.sd)
# 作图
data.heatmap <- rbind(data.Urine.heatmap, data.NIS.heatmap, data.Serum.heatmap)
data.heatmap <- data.heatmap %>% round(., 3)
data.heatmap$Train.Show <- str_c(data.heatmap$Train.Mean, '±', data.heatmap$Train.sd, sep = ' ')
data.heatmap$Test.Show <- str_c(data.heatmap$Test.Mean, '±', data.heatmap$Test.sd, sep = ' ')
data.heatmap$Test.3M.Show <- str_c(data.heatmap$Test.3M.Mean, '±', data.heatmap$Test.3M.sd, sep = ' ')
data.heatmap$Test.10M.Show <- str_c(data.heatmap$Test.10M.Mean, '±', data.heatmap$Test.10M.sd, sep = ' ')
data.heatmap$Test.15M.Show <- str_c(data.heatmap$Test.15M.Mean, '±', data.heatmap$Test.15M.sd, sep = ' ')
data.heatmap[1:5] %>% range()
breakList <- seq(0.3, 1, 0.1)
pheatmap(data.heatmap[1:5], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[11:15], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(11, 22, 33),
         filename = 'Results/Figure3/Mouse.pseudo10_caret.ROC_2.png') # 保存，自动调整纸张大小
pheatmap(data.heatmap[1:5], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[11:15], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(11, 22, 33),
         filename = 'Results/Figure3/Mouse.pseudo10_caret.ROC_2.pdf') # 保存，自动调整纸张大小
data.heatmap[1:2] %>% range()
breakList <- seq(0.55, 1, 0.05)
pheatmap(data.heatmap[1:2], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[11:12], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(11, 22, 33),
         filename = 'Results/Figure3/Mouse.pseudo10_caret.ROC_Train.Test_2.png') # 保存，自动调整纸张大小
pheatmap(data.heatmap[1:2], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[11:12], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(11, 22, 33),
         filename = 'Results/Figure3/Mouse.pseudo10_caret.ROC_Train.Test_2.pdf') # 保存，自动调整纸张大小
data.heatmap[1:11, 1:5] %>% range()
breakList <- seq(0.55, 1, 0.05)
pheatmap(data.heatmap[1:11, 1:5], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[1:11, 11:15], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         filename = 'Results/Figure3/MouseUrine.pseudo10_caret.ROC_2.png') # 保存，自动调整纸张大小
pheatmap(data.heatmap[1:11, 1:5], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[1:11, 11:15], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         filename = 'Results/Figure3/MouseUrine.pseudo10_caret.ROC_2.pdf') # 保存，自动调整纸张大小

# -- ------------------------------------------------------------------------------------------------------------------
# 5.展示每种体液建模效果最好的ROC图和Feature排序
# -- ------------------------------------------------------------------------------------------------------------------
# MouseUrine的3-10-15MROC图
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20_caret.RData')
Fit <- readRDS('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20.glmnet.Fit.rds')
var.importance <- varImp(Fit, scale = T)
png('Results/Figure3/MouseUrine.pseudo10.CPM.quantiles.glmnet.Fit.VarImportance.png', width = 300, height = 400)
plot(var.importance, xlab = 'Importance', top = 10, main = 'glmnet')
dev.off()
pdf('Results/Figure3/MouseUrine.pseudo10.CPM.quantiles.glmnet.Fit.VarImportance.pdf', width = 3, height = 4)
plot(var.importance, xlab = 'Importance', top = 10, main = 'glmnet')
dev.off()
# ROC图
predict.train <- predict(Fit, newdata = train.transformed, type = 'prob')
predict.test <- predict(Fit, newdata = test.transformed, type = 'prob')
test.transformed.3M <- test.transformed %>% filter(Class == 3)
predict.test.3M <- predict(Fit, newdata = test.transformed.3M, type = 'prob')
test.transformed.10M <- test.transformed %>% filter(Class == 10)
predict.test.10M <- predict(Fit, newdata = test.transformed.10M, type = 'prob')
test.transformed.15M <- test.transformed %>% filter(Class == 15)
predict.test.15M <- predict(Fit, newdata = test.transformed.15M, type = 'prob')
roc.train <- roc(train.transformed$Group, predict.train$AD, levels = c('Con','AD'), direction = '<')
roc.test <- roc(test.transformed$Group, predict.test$AD, levels = c('Con','AD'), direction = '<')
roc.test.3M <- roc(test.transformed.3M$Group, predict.test.3M$AD, levels = c('Con','AD'), direction = '<')
roc.test.10M <- roc(test.transformed.10M$Group, predict.test.10M$AD, levels = c('Con','AD'), direction = '<')
roc.test.15M <- roc(test.transformed.15M$Group, predict.test.15M$AD, levels = c('Con','AD'), direction = '<')
png('Results/Figure3/MouseUrine.pseudo10.glmnet.Fit.ROC.png', height = 600, width = 600)
plot.roc(roc.test, col = "#3C5488FF", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#3C5488FF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.5,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.test.3M, col = "#F39B7FFF", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F39B7FFF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.45,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.test.10M, col = "#4DBBD5FF", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#4DBBD5FF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.4,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.test.15M, col = "#DC0000FF", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#DC0000FF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.35,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
pdf('Results/Figure3/MouseUrine.pseudo10.glmnet.Fit.ROC.pdf', height = 6, width = 6)
plot.roc(roc.test, col = "#3C5488FF", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#3C5488FF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.5,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.test.3M, col = "#F39B7FFF", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F39B7FFF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.45,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.test.10M, col = "#4DBBD5FF", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#4DBBD5FF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.4,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.test.15M, col = "#DC0000FF", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#DC0000FF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.35,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
# -- ------------------------------------------------------------------------------------------------------------------
# Mouse三种体液的ROC图
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20_caret.RData')
Fit <- readRDS('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_old_LASSO_Logistic.Pos_Top20.glmnet.Fit.rds')
Urine.predict.test <- predict(Fit, newdata = test.transformed, type = 'prob')
Urine.roc.test <- roc(test.transformed$Group, Urine.predict.test$AD, levels = c('Con','AD'), direction = '<')

load('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20_caret.RData')
Fit <- readRDS('Results/Figure3/RDS/MouseNIS.pseudo10.CPM.quantiles.train_2021_LASSO_Logistic.Pos_Top20.kknn.Fit.rds')
NIS.predict.test <- predict(Fit, newdata = test.transformed, type = 'prob')
NIS.roc.test <- roc(test.transformed$Group, NIS.predict.test$AD, levels = c('Con','AD'), direction = '<')
var.importance <- varImp(Fit, scale = T)
png('Results/Figure3/MouseNIS.pseudo10.CPM.quantiles.kknn.Fit.2022.VarImportance.png', width = 300, height = 400)
plot(var.importance, xlab = 'Importance', top = 10, main = 'kknn')
dev.off()
pdf('Results/Figure3/MouseNIS.pseudo10.CPM.quantiles.kknn.Fit.2022.VarImportance.pdf', width = 3, height = 4)
plot(var.importance, xlab = 'Importance', top = 10, main = 'kknn')
dev.off()

load('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
Fit <- readRDS('Results/Figure3/RDS/MouseSerum.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.rf.Fit.rds')
Serum.predict.test <- predict(Fit, newdata = test.transformed, type = 'prob')
Serum.roc.test <- roc(test.transformed$Group, Serum.predict.test$AD, levels = c('Con','AD'), direction = '<')
var.importance <- varImp(Fit, scale = T)
png('Results/Figure3/MouseSerum.pseudo10.CPM.quantiles.rf.Fit.2022.VarImportance.png', width = 300, height = 400)
plot(var.importance, xlab = 'Importance', top = 10, main = 'rf')
dev.off()
pdf('Results/Figure3/MouseSerum.pseudo10.CPM.quantiles.rf.Fit.2022.VarImportance.pdf', width = 3, height = 4)
plot(var.importance, xlab = 'Importance', top = 10, main = 'rf')
dev.off()
# ROC图
png('Results/Figure3/Mouse.pseudo10.Best.Fit.ROC_2.png', height = 600, width = 600)
plot.roc(Serum.roc.test, col = "#D595A7FF", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#D595A7FF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.4,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(NIS.roc.test, col = "#6BD76BFF", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#6BD76BFF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.45,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(Urine.roc.test, col = "#F0E685FF", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F0E685FF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.5,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
pdf('Results/Figure3/Mouse.pseudo10.Best.Fit.ROC_2.pdf', height = 6, width = 6)
plot.roc(Serum.roc.test, col = "#D595A7FF", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#D595A7FF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.4,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(NIS.roc.test, col = "#6BD76BFF", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#6BD76BFF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.45,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(Urine.roc.test, col = "#F0E685FF", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F0E685FF", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.5,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()


# -- ------------------------------------------------------------------------------------------------------------------
# 20230302 用于建模的蛋白与小鼠Aβ浓度的correlation
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse) # 处理数据神包
library(ggsci) # 期刊配色，单个离散配色模块中最多只有51种颜色
library(scales) # show_col
library(ggthemes) # 主题设置
library(patchwork) # 拼接图
library(ggbreak) # 坐标轴截断
library(ggsignif) # 计算差异显著性并添加显著性线段
library(ggpubr) # stat_cor函数添加R和P
library(ggpmisc) # stat_poly_eq函数添加趋势线和方程
library(ggrepel)
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/MouseUrine_Raw.Counts_CPM.quantiles.csv',row.names = 1)
data <- data %>% t() %>% as.data.frame()
data$Group <- ifelse(str_detect(rownames(data), 'AD'), 'AD', 'Con')
data$Month <- ifelse(str_detect(rownames(data), '15'), 15,
                     ifelse(str_detect(rownames(data), '10'), 10,3))
data <- data %>% relocate(Group, Month)
data <- data %>% filter(Group == 'AD')
names(data) <- names(data) %>% str_replace_all('/', '.') # 文件名中有/的时候会报错
for (i in 1:193) {
  p <- ggplot(data, aes_(x = as.name(colnames(data)[i+2]))) + theme_bw() +
    geom_jitter(aes(y = Month, color = factor(Month)), width = 0.1, size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#F39B7FFF','#4DBBD5FF','#DC0000FF')) +
    stat_smooth(aes(y = Month), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
    stat_cor(aes(y = Month), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = 'darkgreen') +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none") +
    scale_y_continuous(breaks = c(3,10,15), labels = c(3,10,15))
  fn <- str_c('Results/Figure3/Month_Correlation/MouseUrine_Correlation_Abeta vs ', colnames(data)[i+2], '.png')
  ggsave(p, filename = fn, width = 4, height = 4)
}
# data$Abeta <- c(0.056,0.073,0.146,0.036,4.495,4.891,4.242,4.329,6.419,7.364,9.249,5.817)  # 4X
data$Abeta <- c(0.315,0.456,0.601,0.345,14.054,12.29,15.804,13.85,26.376,29.812,34.586,31.771)  # Cortex
data <- data %>% relocate(Group, Month, Abeta)
names(data) <- names(data) %>% str_replace_all('/', '.') # 文件名中有/的时候会报错
for (i in 1:193) {
  p <- ggplot(data, aes_(x = as.name(colnames(data)[i+3]))) + theme_bw() +
    geom_point(aes(y = Abeta, color = factor(Month)), size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#F39B7FFF','#4DBBD5FF','#DC0000FF')) +
    stat_smooth(aes(y = Abeta), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
    stat_cor(aes(y = Abeta), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = 'darkgreen') +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  # fn <- str_c('Results/Figure3/Abeta_Correlation/Cortex//MouseUrine_Correlation_Abeta vs ', colnames(data)[i+3], '.png')
  # ggsave(p, filename = fn, width = 4, height = 4)
  fn <- str_c('Results/Figure3/Abeta_Correlation/Cortex/MouseUrine_Correlation_Abeta vs ', colnames(data)[i+3], '.pdf')
  ggsave(p, filename = fn, width = 4, height = 4)
}
data$Abeta <- c(3.968,5.773,8.295,7.672,14.186,15.225,15.074,14.63,24.715,25.079,26.693,22.708)  # Hippo
data <- data %>% relocate(Group, Month, Abeta)
names(data) <- names(data) %>% str_replace_all('/', '.') # 文件名中有/的时候会报错
for (i in 1:193) {
  p <- ggplot(data, aes_(x = as.name(colnames(data)[i+3]))) + theme_bw() +
    geom_point(aes(y = Abeta, color = factor(Month)), size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#F39B7FFF','#4DBBD5FF','#DC0000FF')) +
    stat_smooth(aes(y = Abeta), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
    stat_cor(aes(y = Abeta), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = 'darkgreen') +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  # fn <- str_c('Results/Figure3/Abeta_Correlation/Hippo//MouseUrine_Correlation_Abeta vs ', colnames(data)[i+3], '.png')
  # ggsave(p, filename = fn, width = 4, height = 4)
  fn <- str_c('Results/Figure3/Abeta_Correlation/Hippo/MouseUrine_Correlation_Abeta vs ', colnames(data)[i+3], '.pdf')
  ggsave(p, filename = fn, width = 4, height = 4)
}
# Tang说把皮层海马的放在一个图中
rm(list=ls())
gc()
data <- read.csv('Results/MouseUrine_Raw.Counts_CPM.quantiles.csv',row.names = 1)
data <- data %>% t() %>% as.data.frame()
data$Group <- ifelse(str_detect(rownames(data), 'AD'), 'AD', 'Con')
data$Month <- ifelse(str_detect(rownames(data), '15'), 15,
                     ifelse(str_detect(rownames(data), '10'), 10,3))
data <- data %>% relocate(Group, Month)
data <- data %>% filter(Group == 'AD')
data$AbetaCortex <- c(0.315,0.456,0.601,0.345,14.054,12.29,15.804,13.85,26.376,29.812,34.586,31.771)  # Cortex
data$AbetaHippo <- c(3.968,5.773,8.295,7.672,14.186,15.225,15.074,14.63,24.715,25.079,26.693,22.708)  # Hippo
data <- data %>% relocate(Group, Month, AbetaCortex, AbetaHippo)
names(data) <- names(data) %>% str_replace_all('/', '.') # 文件名中有/的时候会报错
show_col(pal_aaas(alpha = 1)(10))
i = 19 # CD26
i = 116 # ITGB2
p1 <- ggplot(data, aes_(x = as.name(colnames(data)[i+5]))) + theme_bw() +
  geom_point(aes(y = AbetaCortex), colour = '#EE0000FF', size = 3, alpha = 0.8) + 
  stat_smooth(aes(y = AbetaCortex), method='lm', formula = y~x, colour='#EE0000FF', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(aes(y = AbetaCortex), method = "spearman", digits = 3, size = 4.5, label.y = 4, colour = '#EE0000FF') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.position = "none")
p2 <- p1 + geom_point(data= data, aes(y = AbetaHippo), colour = '#5F559BFF', size = 3, alpha = 0.8) + 
  stat_smooth(aes(y = AbetaHippo), method='lm', formula = y~x, colour='#5F559BFF', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(aes(y = AbetaHippo), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = '#5F559BFF') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.position = "none") + ylab('Aβ Mean IntDen')
fn <- str_c('Results/Figure3/Abeta_Correlation/MouseUrine_Correlation_Abeta vs ', colnames(data)[i+5], '.png')
ggsave(p2, filename = fn, width = 4, height = 4)
fn <- str_c('Results/Figure3/Abeta_Correlation/MouseUrine_Correlation_Abeta vs ', colnames(data)[i+5], '.pdf')
ggsave(p2, filename = fn, width = 4, height = 4)
i = 69 # DSC3
i = 142 # NECTIN1
p1 <- ggplot(data, aes_(x = as.name(colnames(data)[i+5]))) + theme_bw() +
  geom_point(aes(y = AbetaCortex), colour = '#EE0000FF', size = 3, alpha = 0.8) + 
  stat_smooth(aes(y = AbetaCortex), method='lm', formula = y~x, colour='#EE0000FF', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(aes(y = AbetaCortex), method = "spearman", digits = 3, size = 4.5, label.y = 40, colour = '#EE0000FF') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.position = "none")
p2 <- p1 + geom_point(data= data, aes(y = AbetaHippo), colour = '#5F559BFF', size = 3, alpha = 0.8) + 
  stat_smooth(aes(y = AbetaHippo), method='lm', formula = y~x, colour='#5F559BFF', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(aes(y = AbetaHippo), method = "spearman", digits = 3, size = 4.5, label.y = 36, colour = '#5F559BFF') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.position = "none") + ylab('Aβ Mean IntDen')
fn <- str_c('Results/Figure3/Abeta_Correlation/MouseUrine_Correlation_Abeta vs ', colnames(data)[i+5], '.png')
ggsave(p2, filename = fn, width = 4, height = 4)
fn <- str_c('Results/Figure3/Abeta_Correlation/MouseUrine_Correlation_Abeta vs ', colnames(data)[i+5], '.pdf')
ggsave(p2, filename = fn, width = 4, height = 4)
# -- ------------------------------------------------------------------------------------------------------------------
# NIS
# -- ------------------------------------------------------------------------------------------------------------------
# 把皮层海马的放在一个图中
rm(list=ls())
gc()
data <- read.csv('Results/MouseNIS_Raw.Counts_CPM.quantiles.csv',row.names = 1)
data <- data %>% t() %>% as.data.frame()
data$Group <- ifelse(str_detect(rownames(data), 'AD'), 'AD', 'Con')
data$Month <- ifelse(str_detect(rownames(data), '15'), 15,
                     ifelse(str_detect(rownames(data), '10'), 10,3))
data <- data %>% relocate(Group, Month)
data <- data %>% filter(Group == 'AD')
data$AbetaCortex <- c(0.315,0.456,0.601,0.345,14.054,12.29,15.804,13.85,26.376,29.812,34.586,31.771)  # Cortex
data$AbetaHippo <- c(3.968,5.773,8.295,7.672,14.186,15.225,15.074,14.63,24.715,25.079,26.693,22.708)  # Hippo
data <- data %>% relocate(Group, Month, AbetaCortex, AbetaHippo)
names(data) <- names(data) %>% str_replace_all('/', '.') # 文件名中有/的时候会报错
show_col(pal_aaas(alpha = 1)(10))
for (i in 1:193) {
  p1 <- ggplot(data, aes_(x = as.name(colnames(data)[i+5]))) + theme_bw() +
    geom_point(aes(y = AbetaCortex), colour = '#EE0000FF', size = 3, alpha = 0.8) + 
    stat_smooth(aes(y = AbetaCortex), method='lm', formula = y~x, colour='#EE0000FF', linetype = 2, size = 0.5, alpha = 0.2) +
    stat_cor(aes(y = AbetaCortex), method = "spearman", digits = 3, size = 4.5, label.y = 4, colour = '#EE0000FF') +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  p2 <- p1 + geom_point(data= data, aes(y = AbetaHippo), colour = '#5F559BFF', size = 3, alpha = 0.8) + 
    stat_smooth(aes(y = AbetaHippo), method='lm', formula = y~x, colour='#5F559BFF', linetype = 2, size = 0.5, alpha = 0.2) +
    stat_cor(aes(y = AbetaHippo), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = '#5F559BFF') +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none") + ylab('Aβ Mean IntDen')
  fn <- str_c('Results/Figure3/Abeta_Correlation/NIS/MouseNIS_Correlation_Abeta vs ', colnames(data)[i+3], '.png')
  ggsave(p2, filename = fn, width = 4, height = 4)
  # fn <- str_c('Results/Figure3/Abeta_Correlation/NIS/MouseNIS_Correlation_Abeta vs ', colnames(data)[i+3], '.pdf')
  # ggsave(p2, filename = fn, width = 4, height = 4)
}
i = 53 # CLDN19
i = 86 # GPC1
p1 <- ggplot(data, aes_(x = as.name(colnames(data)[i+5]))) + theme_bw() +
  geom_point(aes(y = AbetaCortex), colour = '#EE0000FF', size = 3, alpha = 0.8) + 
  stat_smooth(aes(y = AbetaCortex), method='lm', formula = y~x, colour='#EE0000FF', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(aes(y = AbetaCortex), method = "spearman", digits = 3, size = 4.5, label.y = 40, colour = '#EE0000FF') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.position = "none")
p2 <- p1 + geom_point(data= data, aes(y = AbetaHippo), colour = '#5F559BFF', size = 3, alpha = 0.8) + 
  stat_smooth(aes(y = AbetaHippo), method='lm', formula = y~x, colour='#5F559BFF', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(aes(y = AbetaHippo), method = "spearman", digits = 3, size = 4.5, label.y = 36, colour = '#5F559BFF') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.position = "none") + ylab('Aβ Mean IntDen')
fn <- str_c('Results/Figure3/Abeta_Correlation/MouseNIS_Correlation_Abeta vs ', colnames(data)[i+5], '.png')
ggsave(p2, filename = fn, width = 4, height = 4)
fn <- str_c('Results/Figure3/Abeta_Correlation/MouseNIS_Correlation_Abeta vs ', colnames(data)[i+5], '.pdf')
ggsave(p2, filename = fn, width = 4, height = 4)
i = 168 # SIGLEC1
i = 163 # PROM1
p1 <- ggplot(data, aes_(x = as.name(colnames(data)[i+5]))) + theme_bw() +
  geom_point(aes(y = AbetaCortex), colour = '#EE0000FF', size = 3, alpha = 0.8) + 
  stat_smooth(aes(y = AbetaCortex), method='lm', formula = y~x, colour='#EE0000FF', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(aes(y = AbetaCortex), method = "spearman", digits = 3, size = 4.5, label.y = 4, colour = '#EE0000FF') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.position = "none")
p2 <- p1 + geom_point(data= data, aes(y = AbetaHippo), colour = '#5F559BFF', size = 3, alpha = 0.8) + 
  stat_smooth(aes(y = AbetaHippo), method='lm', formula = y~x, colour='#5F559BFF', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(aes(y = AbetaHippo), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = '#5F559BFF') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.position = "none") + ylab('Aβ Mean IntDen')
fn <- str_c('Results/Figure3/Abeta_Correlation/MouseNIS_Correlation_Abeta vs ', colnames(data)[i+5], '.png')
ggsave(p2, filename = fn, width = 4, height = 4)
fn <- str_c('Results/Figure3/Abeta_Correlation/MouseNIS_Correlation_Abeta vs ', colnames(data)[i+5], '.pdf')
ggsave(p2, filename = fn, width = 4, height = 4)
# -- ------------------------------------------------------------------------------------------------------------------
# Serum
# -- ------------------------------------------------------------------------------------------------------------------
# 把皮层海马的放在一个图中
rm(list=ls())
gc()
data <- read.csv('Results/MouseSerum_Raw.Counts_CPM.quantiles.csv',row.names = 1)
data <- data %>% t() %>% as.data.frame()
data$Group <- ifelse(str_detect(rownames(data), 'AD'), 'AD', 'Con')
data$Month <- ifelse(str_detect(rownames(data), '15'), 15,
                     ifelse(str_detect(rownames(data), '10'), 10,3))
data <- data %>% relocate(Group, Month)
data <- data %>% filter(Group == 'AD')
data$AbetaCortex <- c(0.315,0.456,0.601,0.345,14.054,12.29,15.804,13.85,26.376,29.812,34.586,31.771)  # Cortex
data$AbetaHippo <- c(3.968,5.773,8.295,7.672,14.186,15.225,15.074,14.63,24.715,25.079,26.693,22.708)  # Hippo
data <- data %>% relocate(Group, Month, AbetaCortex, AbetaHippo)
names(data) <- names(data) %>% str_replace_all('/', '.') # 文件名中有/的时候会报错
show_col(pal_aaas(alpha = 1)(10))
for (i in 1:193) {
  p1 <- ggplot(data, aes_(x = as.name(colnames(data)[i+5]))) + theme_bw() +
    geom_point(aes(y = AbetaCortex), colour = '#EE0000FF', size = 3, alpha = 0.8) + 
    stat_smooth(aes(y = AbetaCortex), method='lm', formula = y~x, colour='#EE0000FF', linetype = 2, size = 0.5, alpha = 0.2) +
    stat_cor(aes(y = AbetaCortex), method = "spearman", digits = 3, size = 4.5, label.y = 4, colour = '#EE0000FF') +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  p2 <- p1 + geom_point(data= data, aes(y = AbetaHippo), colour = '#5F559BFF', size = 3, alpha = 0.8) + 
    stat_smooth(aes(y = AbetaHippo), method='lm', formula = y~x, colour='#5F559BFF', linetype = 2, size = 0.5, alpha = 0.2) +
    stat_cor(aes(y = AbetaHippo), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = '#5F559BFF') +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none") + ylab('Aβ Mean IntDen')
  fn <- str_c('Results/Figure3/Abeta_Correlation/Serum/MouseSerum_Correlation_Abeta vs ', colnames(data)[i+3], '.png')
  ggsave(p2, filename = fn, width = 4, height = 4)
  # fn <- str_c('Results/Figure3/Abeta_Correlation/Serum/MouseSerum_Correlation_Abeta vs ', colnames(data)[i+3], '.pdf')
  # ggsave(p2, filename = fn, width = 4, height = 4)
}
i = 47 # CEACAM8
p1 <- ggplot(data, aes_(x = as.name(colnames(data)[i+5]))) + theme_bw() +
  geom_point(aes(y = AbetaCortex), colour = '#EE0000FF', size = 3, alpha = 0.8) + 
  stat_smooth(aes(y = AbetaCortex), method='lm', formula = y~x, colour='#EE0000FF', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(aes(y = AbetaCortex), method = "spearman", digits = 3, size = 4.5, label.y = 40, colour = '#EE0000FF') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.position = "none")
p2 <- p1 + geom_point(data= data, aes(y = AbetaHippo), colour = '#5F559BFF', size = 3, alpha = 0.8) + 
  stat_smooth(aes(y = AbetaHippo), method='lm', formula = y~x, colour='#5F559BFF', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(aes(y = AbetaHippo), method = "spearman", digits = 3, size = 4.5, label.y = 36, colour = '#5F559BFF') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.position = "none") + ylab('Aβ Mean IntDen')
fn <- str_c('Results/Figure3/Abeta_Correlation/MouseSerum_Correlation_Abeta vs ', colnames(data)[i+5], '.png')
ggsave(p2, filename = fn, width = 4, height = 4)
fn <- str_c('Results/Figure3/Abeta_Correlation/MouseSerum_Correlation_Abeta vs ', colnames(data)[i+5], '.pdf')
ggsave(p2, filename = fn, width = 4, height = 4)
i = 87 # GSN
i = 115 # ITGB1
p1 <- ggplot(data, aes_(x = as.name(colnames(data)[i+5]))) + theme_bw() +
  geom_point(aes(y = AbetaCortex), colour = '#EE0000FF', size = 3, alpha = 0.8) + 
  stat_smooth(aes(y = AbetaCortex), method='lm', formula = y~x, colour='#EE0000FF', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(aes(y = AbetaCortex), method = "spearman", digits = 3, size = 4.5, label.y = 4, colour = '#EE0000FF') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.position = "none")
p2 <- p1 + geom_point(data= data, aes(y = AbetaHippo), colour = '#5F559BFF', size = 3, alpha = 0.8) + 
  stat_smooth(aes(y = AbetaHippo), method='lm', formula = y~x, colour='#5F559BFF', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(aes(y = AbetaHippo), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = '#5F559BFF') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(colour = 'black', size = 16),
        legend.position = "none") + ylab('Aβ Mean IntDen')
fn <- str_c('Results/Figure3/Abeta_Correlation/MouseSerum_Correlation_Abeta vs ', colnames(data)[i+5], '.png')
ggsave(p2, filename = fn, width = 4, height = 4)
fn <- str_c('Results/Figure3/Abeta_Correlation/MouseSerum_Correlation_Abeta vs ', colnames(data)[i+5], '.pdf')
ggsave(p2, filename = fn, width = 4, height = 4)

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
  fn2 <- str_c('Results/Figure3/HumanUrine_pseudo10.CPM.quantiles_2021.nb_Train_', myTitle, '.png')
  ggsave(p2, filename = fn2, width = 4, height = 4)
  fn2 <- str_c('Results/Figure3/HumanUrine_pseudo10.CPM.quantiles_2021.nb_Train_', myTitle, '.pdf')
  ggsave(p2, filename = fn2, width = 4, height = 4)
}



