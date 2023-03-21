# 设置最大的内存
if(.Platform$OS.type == "windows") withAutoprint({
  memory.size()
  memory.size(TRUE)
  memory.limit()
})

# -- ------------------------------------------------------------------------------------------------------------------
# 1.Pseudo processing of samples
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
sampleID <- read.csv('Mouse_SampleID_SampleName.csv')
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

# -- ------------------------------------------------------------------------------------------------------------------
# 2.Correlation of pseudo samples with raw
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read_csv('Results/Mouse_Raw.Counts.csv', show_col_types = F)
names(data)[1] <- 'Pro'
data.anno <- read_csv('Mouse_SampleID_SampleName.csv', show_col_types = F)
for (i in 1:nrow(data.anno)) {
  dt <- data %>% dplyr::select(Pro, data.anno$SampleName[i])
  fn1 <- str_c('pseudoData/pseudo10/', data.anno$SampleID[i], '.pseudo10.csv')
  temp <- read.csv(fn1, row.names = 1)
  temp <- t(temp) %>% as.data.frame()
  temp$Pro <- rownames(temp)
  dt <- inner_join(dt, temp, by = 'Pro')
  dt.cor <- cor(dt[-1])
  fn2 <- str_c('Results/FigureS4/', data.anno$SampleName[i], '_Raw.pseudo10.png')
  pheatmap(dt.cor, color = pal_material(palette = c("orange"), n = 20)(20), cellwidth = 10, cellheight = 6,
           display_numbers = T, number_color = 'black', number_format = '%.3f', 
           fontsize_number = 3, fontsize_row = 6, fontsize_col = 8,
           cluster_rows = F, cluster_cols = F,
           filename = fn2)
  fn2 <- str_c('Results/FigureS4/', data.anno$SampleName[i], '_Raw.pseudo10.pdf')
  pheatmap(dt.cor, color = pal_material(palette = c("orange"), n = 20)(20), cellwidth = 10, cellheight = 6,
           display_numbers = T, number_color = 'black', number_format = '%.3f', 
           fontsize_number = 3, fontsize_row = 6, fontsize_col = 8,
           cluster_rows = F, cluster_cols = F,
           filename = fn2)
}

# -- ------------------------------------------------------------------------------------------------------------------
# 3.Machine learning
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse) 
library(glmnet)
library(edgeR) 
library(preprocessCore) 
library(caret)
library(pROC)
library(gbm)
# -- ------------------------------------------------------------------------------------------------------------------
# MouseUrine
# -- ------------------------------------------------------------------------------------------------------------------
# Synthesize the total expression matrix
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
# Standardization of data and annotation of sample groups
rm(list=ls())
gc()
data <- read.csv('pseudoData/MouseUrine.AllSample.pseudo10.csv', row.names = 1)
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
data$Class <- ifelse(str_detect(data$ID, '15'), 15,
                        ifelse(str_detect(data$ID, '10'), 10, 3))
data <- relocate(data, ID, Group, Class)
write.csv(data, 'pseudoData/MouseUrine.AllSample.pseudo10.CPM.quantiles.Group.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# Machine learning with the caret package
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('pseudoData/MouseUrine.AllSample.pseudo10.CPM.quantiles.Group.csv', row.names = 1)
data$ID <- data$ID %>% factor()
data$Group <- data$Group %>% factor()
data$Class <- data$Class %>% factor(levels = c(3,10,15))
# -- ------------------------------------------------------------------------------------------------------------------
# a.Dividing the training set and the test set:
# -- ------------------------------------------------------------------------------------------------------------------
data.ID <- levels(data$ID)
data.ID <- data.frame(ID = data.ID)
data.ID$Group <- ifelse(str_detect(data.ID$ID, 'D'), 'AD', 'Con')
data.ID$Class <- ifelse(str_detect(data.ID$ID, '15'), 15,
                     ifelse(str_detect(data.ID$ID, '10'), 10, 3))
data.ID$Class <- factor(data.ID$Class, levels = c(3,10,15), labels = c('3M', '10M', '15M'))
data.ID$Split <- str_c(data.ID$Class, data.ID$Group, sep = '.')
data.ID$Split <- factor(data.ID$Split, levels = c('3M.AD', '3M.Con', '10M.AD', '10M.Con', '15M.AD', '15M.Con'))
set.seed(1234)
train.ID <- createDataPartition(data.ID$Split, p = 0.5, list = FALSE)
train.ID <- data.ID$ID[train.ID]
data.train <- data %>% filter(ID %in% train.ID) %>% dplyr::select(!ID)
data.test <- data %>% filter(!ID %in% train.ID) %>% dplyr::select(!ID)
# -- ------------------------------------------------------------------------------------------------------------------
# b.Filtering variables with LASSO
# -- ------------------------------------------------------------------------------------------------------------------
X <- as.matrix(data.train[3:ncol(data.train)])
Y <- as.matrix(data.train[1])
# glmnet() syntax with alpha=0 for ridge regression and 1 for LASSO regression
myLasso <- glmnet(X, Y, alpha = 1, family = 'binomial', nlambda = 200) 
NROW(myLasso$lambda)
min(myLasso$lambda) # Choose the lambda value of the optimal model
pdf('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO.pdf', width = 6, height = 8)
plot(myLasso, xvar = 'lambda', label = T)
dev.off()
lasso.coef <- coef(myLasso, s = min(myLasso$lambda)) # Regression coefficients under the optimal solution
my.lasso.coef <- as.matrix(lasso.coef) %>% as.data.frame()
my.lasso.coef$abs <- abs(my.lasso.coef$s1)
my.lasso.coef <- arrange(my.lasso.coef, desc(abs)) %>% filter(abs > 0)
write.csv(my.lasso.coef, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO.csv')
# Single-factor logistic regression was performed with the variables screened by LASSO to obtain the p-value of each facto
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
res2 <- res1 %>% filter(`Pr(>|z|)` < 0.05)
write.csv(res2, file = 'Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# c.Construction of Train and Test datasets with positive variables
# -- ------------------------------------------------------------------------------------------------------------------
SigPro <- read.csv('Results/Figure3/Files/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos.csv', row.names = 1)
myPro <- str_replace_all(rownames(SigPro)[1:20], '/', '.') %>% str_replace_all('-', '.')
myPro <- intersect(myPro, names(data.train))
data.train <- data.train %>% dplyr::select(Class, Group, all_of(myPro))
data.test <- data.test %>% dplyr::select(Class, Group, all_of(myPro))
# Normalize the data and fill in the missing values, this can be done with the preProcess command
preProcValues <- preProcess(data.train, method = c('center','scale')) # Normalization for validation datasets
train.transformed <- predict(preProcValues, data.train)
test.transformed <- predict(preProcValues, data.test)
str(train.transformed)
# SAve data
save(data.ID, train.ID, data.train, data.test, preProcValues, train.transformed, test.transformed,
     file = 'Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
# -- ------------------------------------------------------------------------------------------------------------------
# d.Model Training and Tuning
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20_caret.RData')
myMethod <- c('glmnet','kknn','C5.0','rf','AdaBag','gbm','nnet','svmLinear','svmPoly','svmRadial','nb')
fitControl = trainControl(method = 'cv', number = 5,
                          classProbs = T, summaryFunction = twoClassSummary,
                          search = 'random') # Random transfer of parameters
for (i in 1:length(myMethod)) {
  set.seed(2022)
  Fit <- train(Group ~ ., data = train.transformed[-1], method = myMethod[i],
               trControl = fitControl, 
               tuneLength = 30, # Random number of random transfers of parameters
               metric = 'ROC', verbose = F) 
  fn <- str_c('Results/Figure3/RDS/MouseUrine.pseudo10.CPM.quantiles.train_1234_LASSO_Logistic.Pos_Top20.', myMethod[i], '.Fit.rds')
  saveRDS(Fit, fn)
}
# -- ------------------------------------------------------------------------------------------------------------------
# e.Model prediction and evaluation
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

# The analysis of other groups is similar to Urine

# -- ------------------------------------------------------------------------------------------------------------------
# Correlation of protein used for modeling with mouse Aβ concentration
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse) 
library(ggsci) 
library(scales) # show_col
library(ggthemes) 
library(ggsignif) # Calculate the significance of differences and add significant line segments
library(ggpubr) # The stat_cor function for R and P
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
data$AbetaCortex <- c(0.315,0.456,0.601,0.345,14.054,12.29,15.804,13.85,26.376,29.812,34.586,31.771)  # Cortex
data$AbetaHippo <- c(3.968,5.773,8.295,7.672,14.186,15.225,15.074,14.63,24.715,25.079,26.693,22.708)  # Hippo
data <- data %>% relocate(Group, Month, AbetaCortex, AbetaHippo)
names(data) <- names(data) %>% str_replace_all('/', '.') # An error will be reported if the file name has / in it
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



