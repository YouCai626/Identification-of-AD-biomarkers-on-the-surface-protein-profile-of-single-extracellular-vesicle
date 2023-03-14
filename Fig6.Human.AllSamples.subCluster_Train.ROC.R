# -- ------------------------------------------------------------------------------------------------------------------
getwd()
# 设置最大的内存
if(.Platform$OS.type == "windows") withAutoprint({
  memory.size()
  memory.size(TRUE)
  memory.limit()
})
# options(future.globals.maxSize = 30*1024^3) #设置为30G
setwd("D:/2.AD.PBA.CAM.Article.V3.Proj")

# -- ------------------------------------------------------------------------------------------------------------------
# 对每一个亚群进行ROC分析。
library(tidyverse) # 处理数据神包
library(glmnet) # 进行正则技术所需
library(pROC) # ROC分析与绘图
library(ggsci) # 期刊配色，单个离散配色模块中最多只有51种颜色
library(scales) # show_col
library(patchwork) # 拼接图
library(ggbreak) # 坐标轴截断
library(ggsignif) # 计算差异显著性并添加显著性线段
library(ggpubr) # stat_cor函数添加R和P
library(ggpmisc) # stat_poly_eq函数添加趋势线和方程
library(ggrepel)
library(ggthemes) # 主题设置
library(pheatmap) # 热图

# -- ------------------------------------------------------------------------------------------------------------------
# 根据Tang的要求对Train的数据进行ROC分析，并做1-5个蛋白的联合ROC的AUC趋势图
dir.create('Results/Figure6')
dir.create('Results/Figure6/Files')
# -- ------------------------------------------------------------------------------------------------------------------
# HumanUrine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
for( j in 0:20) {
  fn1 <- str_c('Results/HumanUrine/subCluster/Expression/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster',j,'.Expression.SCT.csv')
  data <- read.csv(fn1, row.names = 1)
  data <- data %>% t() %>% as.data.frame()
  data$ID <- rownames(data) %>% factor
  data$Group <- ifelse(str_detect(rownames(data), 'D'),'AD','Con') %>% factor(levels = c('AD','Con'))
  data <- relocate(data, ID, Group)
  # 划分训练集与测试集:
  data.ID <- levels(data$ID) # 只有factor才能用levels这个函数
  data.ID <- data.frame(ID = data.ID)
  data.ID$Group <- ifelse(str_detect(data.ID$ID, 'D'), 'AD', 'Con')
  set.seed(2022)
  train.ID <- createDataPartition(data.ID$Group, p = 0.6, list = FALSE)
  train.ID <- data.ID$ID[train.ID]
  data.train <- data %>% filter(ID %in% train.ID) %>% dplyr::select(!ID)
  data.test <- data %>% filter(!ID %in% train.ID) %>% dplyr::select(!ID)
  # 去掉全为0的列：
  data.train <- data.train %>% select_if(is.numeric) %>% select_if(~sum(.) > 0)
  data.train$Group <- ifelse(str_detect(rownames(data.train), 'D'),'AD','Con') %>% factor(levels = c('AD','Con'))
  data.train <- relocate(data.train, Group)
  fn2 <- str_c('Results/HumanUrine/subCluster/caret/RDS/pseudoPer10.Res0.1.Cluster',j,'_Train.2022_Logistic.csv')
  dt <- read.csv(fn2, row.names = 1)
  myPro <- rownames(dt)[1:5]
  roc.data <- data.train %>% select(Group, all_of(myPro))
  roc.data$Group <- roc.data$Group %>% factor(levels = c('Con','AD'))
  roc.data <- relocate(roc.data, Group)
  f2 <- glm(Group ~ roc.data[,2] + roc.data[,3], data = roc.data, family = binomial())
  roc.data$Pro2 <- predict(f2, newdata = roc.data, type = "response")
  f3 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4], data = roc.data, family = binomial())
  roc.data$Pro3 <- predict(f3, newdata = roc.data, type = "response")
  f4 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4] + roc.data[,5], data = roc.data, family = binomial())
  roc.data$Pro4 <- predict(f4, newdata = roc.data, type = "response")
  f5 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4] + roc.data[,5] + roc.data[,6], data = roc.data, family = binomial())
  roc.data$Pro5 <- predict(f5, newdata = roc.data, type = "response")
  auc <- c()
  z <- c()
  pvalue <- c()
  rocPlot1 <- roc(roc.data[,1] ~ roc.data[,2], levels = c('Con','AD'), direction = '<', data = roc.data)
  auc[1] <- rocPlot1$auc
  se <- sqrt(pROC::var(rocPlot1))
  b <- auc[1] - 0.5
  z[1] <- (b / se)
  pvalue[1] <- 2 * pt(-abs(z[1]), df=Inf)
  for(a in 2:5){
    rocPlot <- roc(Group ~ roc.data[, a + 5], levels = c('Con','AD'), direction = '<', data = roc.data)
    auc[a] <- rocPlot$auc
    se <- sqrt(pROC::var(rocPlot))
    b <- auc[a] - 0.5
    z[a] <- (b / se)
    pvalue[a] <- 2 * pt(-abs(z[a]), df=Inf)
  }
  auc_pvalue <- data.frame(auc, z, pvalue)
  rownames(auc_pvalue) <- str_c('Pro', 1:5)
  fn3 <- str_c('Results/Figure6/Files/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster',j,'.Train_2022.ROC.csv')
  write.csv(auc_pvalue, file = fn3)
}
# -- ------------------------------------------------------------------------------------------------------------------
# 将各个亚群的Train ROC合并成一个总的csv
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
filenames <- Sys.glob('Results/Figure6/Files/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster*.Train_2022.ROC.csv')
clusters <- filenames %>% str_remove('Results/Figure6/Files/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster') %>%
  str_remove('.Train_2022.ROC.csv')
data <- read.csv(filenames[1])
data$Cluster <- clusters[1]
for (i in 2:length(clusters)) {
  temp <- read.csv(filenames[i])
  if (nrow(temp) > 0) {
    temp$Cluster <- clusters[i]
    data <- rbind(data, temp)
  }
}
write_csv(data, file = 'Results/Figure6/HumanUrine.pseudoPer10_SCT.Harmony.Dim1-20.0.1.subClusters_Train_2022.ROC.csv')

# -- ------------------------------------------------------------------------------------------------------------------
# 将各个亚群的Train ROC合并成一个总的csv
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data1 <- read.csv('Results/Figure6/HumanUrine.pseudoPer10_SCT.Harmony.Dim1-20.0.1.subClusters_Train_5678.ROC.csv')
data1$ID <- str_c(data1$Cluster, data1$X, sep = '_')
data2 <- read.csv('Results/Figure6/HumanUrine.pseudoPer10_SCT.Harmony.Dim1-20.0.1.subClusters_Train_1234.ROC.csv')
data2$ID <- str_c(data2$Cluster, data2$X, sep = '_')
data3 <- read.csv('Results/Figure6/HumanUrine.pseudoPer10_SCT.Harmony.Dim1-20.0.1.subClusters_Train_2020.ROC.csv')
data3$ID <- str_c(data3$Cluster, data3$X, sep = '_')
data <- full_join(data1, data2, by = 'ID')
data <- full_join(data, data3, by = 'ID')
data <- data %>% dplyr::select(ID, Cluster, X, auc, auc.x, auc.y)
data$Mean <- rowMeans(data[4:6])
data$Mean2 <- ifelse(data$Mean > 0.5, data$Mean, 1-data$Mean)
range(data$Mean2)
data <- data %>% arrange(Cluster, X)
dt <- data %>% group_by(X) %>% mutate(Mean.Mean = mean(Mean))
# 作图
show_col(pal_igv(alpha = 0.8)(56))
# ggplot2绘图点的形状不够用怎么办？
shape_level <- nlevels(factor(data$Cluster))
if (shape_level < 15){
  shapes = (0:shape_level) %% 15
} else{
  shapes = c(0:14, c((15:shape_level) %% 110 + 18))
} 
p <- ggplot(data[1:50,], mapping = aes(x = X, y = Mean2)) + theme_bw() + 
  geom_point(aes(colour = factor(Cluster), shape = factor(Cluster)), size = 2) + 
  scale_shape_manual(values = shapes) +
  geom_line(aes(colour = factor(Cluster), group = factor(Cluster)), size = 0.5) + 
  scale_color_igv() + 
  xlab('Protein Number') + ylab('Mean AUC of 3 Random') + 
  theme(legend.position="none")
ggsave(p, filename = 'Results/Figure6/HumanUrine.pseudoPer10.20-0.1.subClusters.10_Train.ROC.AUC.png', 
       width = 3.6, height = 3.6)
ggsave(p, filename = 'Results/Figure6/HumanUrine.pseudoPer10.20-0.1.subClusters.10_Train.ROC.AUC.pdf', 
       width = 3.6, height = 3.6)
p <- ggplot(data, mapping = aes(x = X, y = Mean2)) + theme_bw() + 
  geom_point(aes(colour = factor(Cluster), shape = factor(Cluster)), size = 2) + 
  scale_shape_manual(values = shapes) +
  geom_line(aes(colour = factor(Cluster), group = factor(Cluster)), size = 0.5) + 
  scale_color_igv() + 
  xlab('Protein Number') + ylab('Mean AUC of 3 Random') + 
  theme(legend.position="none")
ggsave(p, filename = 'Results/Figure6/HumanUrine.pseudoPer10.20-0.1.subClusters_Train.ROC.AUC.png', 
       width = 3.6, height = 3.6)
ggsave(p, filename = 'Results/Figure6/HumanUrine.pseudoPer10.20-0.1.subClusters_Train.ROC.AUC.pdf', 
       width = 3.6, height = 3.6)
p <- ggplot(dt, mapping = aes(x = X)) + theme_bw() + 
  geom_point(aes(y = Mean2, colour = factor(Cluster), shape = factor(Cluster)), size = 2) + 
  scale_shape_manual(values = shapes) +
  geom_line(aes(y = Mean2, colour = factor(Cluster), group = factor(Cluster)), size = 0.5) + 
  scale_color_igv() + 
  xlab('Protein Number') + ylab('Mean AUC of 3 Random') + 
  theme(legend.position="none")
dt2 <- dt %>% select(X, Mean.Mean) %>% unique()
p <- p + geom_col(data = dt2, aes(x = X, y = Mean.Mean), width = 0.5, alpha = 0.5) + 
  coord_cartesian(ylim = c(0.6, 1))
ggsave(p, filename = 'Results/Figure6/HumanUrine.pseudoPer10.20-0.1.subClusters_Train.ROC.AUC.png', 
       width = 3.6, height = 3.6)
ggsave(p, filename = 'Results/Figure6/HumanUrine.pseudoPer10.20-0.1.subClusters_Train.ROC.AUC.pdf', 
       width = 3.6, height = 3.6)











# -- ------------------------------------------------------------------------------------------------------------------
# 将各个亚群的caret ROC合并成一个总的csv，求Mean
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data1 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.5678.caret.ROC.csv')
data1$ID <- str_c(data1$Cluster, data1$method, sep = '_')
data2 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.2020.caret.ROC.csv')
data2$ID <- str_c(data2$Cluster, data2$method, sep = '_')
data3 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.1234.caret.ROC.csv')
data3$ID <- str_c(data3$Cluster, data3$method, sep = '_')
# data4 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.2021.caret.ROC.csv')
# data4$ID <- str_c(data4$Cluster, data4$method, sep = '_')
# data5 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.2022.caret.ROC.csv')
# data5$ID <- str_c(data5$Cluster, data5$method, sep = '_')
data <- full_join(data1, data2, by = 'ID')
data <- full_join(data, data3, by = 'ID')
# data <- full_join(data, data4, by = 'ID')
# data <- full_join(data, data5, by = 'ID')
# data <- data %>% dplyr::select(ID, Cluster, method, Train, Train.x, Train.y, Train.x.x, Train.y.y, Test, Test.x, Test.y, Test.x.x, Test.y.y)
# -- ------------------------------------------------------------------------------------------------------------------
# 取几个数中最大的3个
# df <- apply(data[9:13], 1, function(x){sort(x, decreasing=T)[1:3]}) %>% t() %>% as.data.frame()
# names(df) <- str_c('Test', 1:3, sep = '.')
# df$Test.Mean <- rowMeans(df)
# dt <- cbind(data, df)
# dt <- dt %>% arrange(desc(Test.Mean))
# -- ------------------------------------------------------------------------------------------------------------------
data <- data %>% dplyr::select(ID, Cluster, method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data$Train.Mean <- rowMeans(data[4:6])
data$Train.sd <- apply(data[4:6], 1, sd)
data$Test.Mean <- rowMeans(data[7:9])
data$Test.sd <- apply(data[7:9], 1, sd)
data <- data %>% arrange(desc(Test.Mean), Cluster, method)
data.Urine <- filter(data, Test.Mean > 0.75)
rownames(data.Urine) <- data.Urine$ID
data.Urine.heatmap <- data.Urine %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)


data1 <- read.csv('Results/HumanSaliva/subCluster/caret/AllSamples.pseudoPer20.Res0.1.subClusters_Train.1234.caret.ROC.csv')
data1$ID <- str_c(data1$Cluster, data1$method, sep = '_')
data2 <- read.csv('Results/HumanSaliva/subCluster/caret/AllSamples.pseudoPer20.Res0.1.subClusters_Train.2022.caret.ROC.csv')
data2$ID <- str_c(data2$Cluster, data2$method, sep = '_')
data3 <- read.csv('Results/HumanSaliva/subCluster/caret/AllSamples.pseudoPer20.Res0.1.subClusters_Train.5678.caret.ROC.csv')
data3$ID <- str_c(data3$Cluster, data3$method, sep = '_')
data <- full_join(data1, data2, by = 'ID')
data <- full_join(data, data3, by = 'ID')
data <- data %>% dplyr::select(ID, Cluster, method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data$Train.Mean <- rowMeans(data[4:6])
data$Train.sd <- apply(data[4:6], 1, sd)
data$Test.Mean <- rowMeans(data[7:9])
data$Test.sd <- apply(data[7:9], 1, sd)
data <- data %>% arrange(desc(Test.Mean), Cluster, method)
data.Saliva <- filter(data, Test.Mean > 0.6)
rownames(data.Saliva) <- data.Saliva$ID
data.Saliva.heatmap <- data.Saliva %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)


data1 <- read.csv('Results/HumanPlasma/subCluster/caret/AllSamples.pseudoPer30.Res0.3.subClusters_Train.2020.caret.ROC.csv')
data1$ID <- str_c(data1$Cluster, data1$method, sep = '_')
data2 <- read.csv('Results/HumanPlasma/subCluster/caret/AllSamples.pseudoPer30.Res0.3.subClusters_Train.2022.caret.ROC.csv')
data2$ID <- str_c(data2$Cluster, data2$method, sep = '_')
data3 <- read.csv('Results/HumanPlasma/subCluster/caret/AllSamples.pseudoPer30.Res0.3.subClusters_Train.5678.caret.ROC.csv')
data3$ID <- str_c(data3$Cluster, data3$method, sep = '_')
data <- full_join(data1, data2, by = 'ID')
data <- full_join(data, data3, by = 'ID')
data <- data %>% dplyr::select(ID, Cluster, method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data$Train.Mean <- rowMeans(data[4:6])
data$Train.sd <- apply(data[4:6], 1, sd)
data$Test.Mean <- rowMeans(data[7:9])
data$Test.sd <- apply(data[7:9], 1, sd)
data <- data %>% arrange(desc(Test.Mean), Cluster, method)
data.Plasma <- filter(data, Test.Mean > 0.65)
rownames(data.Plasma) <- data.Plasma$ID
data.Plasma.heatmap <- data.Plasma %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)


data1 <- read.csv('Results/HumanTear/subCluster/caret/AllSamples.pseudoPer10.Res0.3.subClusters_Train.12345.caret.ROC.csv')
data1$ID <- str_c(data1$Cluster, data1$method, sep = '_')
data2 <- read.csv('Results/HumanTear/subCluster/caret/AllSamples.pseudoPer10.Res0.3.subClusters_Train.2022.caret.ROC.csv')
data2$ID <- str_c(data2$Cluster, data2$method, sep = '_')
data3 <- read.csv('Results/HumanTear/subCluster/caret/AllSamples.pseudoPer10.Res0.3.subClusters_Train.5678.caret.ROC.csv')
data3$ID <- str_c(data3$Cluster, data3$method, sep = '_')
data <- full_join(data1, data2, by = 'ID')
data <- full_join(data, data3, by = 'ID')
data <- data %>% dplyr::select(ID, Cluster, method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data$Train.Mean <- rowMeans(data[4:6])
data$Train.sd <- apply(data[4:6], 1, sd)
data$Test.Mean <- rowMeans(data[7:9])
data$Test.sd <- apply(data[7:9], 1, sd)
data <- data %>% arrange(desc(Test.Mean), Cluster, method)
data.Tear <- filter(data, Test.Mean > 0.65)
rownames(data.Tear) <- data.Tear$ID
data.Tear.heatmap <- data.Tear %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)

# 作图
data.heatmap <- rbind(data.Urine.heatmap, data.Saliva.heatmap, data.Plasma.heatmap, data.Tear.heatmap)
data.heatmap <- data.heatmap %>% round(., 3)
data.heatmap$Train.Show <- str_c(data.heatmap$Train.Mean, '±', data.heatmap$Train.sd, sep = ' ')
data.heatmap$Test.Show <- str_c(data.heatmap$Test.Mean, '±', data.heatmap$Test.sd, sep = ' ')
data.heatmap[1:2] %>% range()
breakList <- seq(0.55, 1, 0.05)
pheatmap(data.heatmap[1:2], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         # display_numbers = T, 
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(10, 24, 44),
         filename = 'Results/Figure6/Human.subClusters_caret.ROC.png') # 保存，自动调整纸张大小
pheatmap(data.heatmap[1:2], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         # display_numbers = T, 
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(10, 24, 44),
         filename = 'Results/Figure6/Human.subClusters_caret.ROC.pdf') # 保存，自动调整纸张大小


# -- ------------------------------------------------------------------------------------------------------------------
# 展示最大AUC亚群的AD vs Con的差异
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data1 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.1234_Logistic.Pos.csv')
data2 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.5678_Logistic.Pos.csv')
data3 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.2020_Logistic.Pos.csv')
# Cluster == 3
data1 <- data1 %>% filter(Cluster == 3)
data2 <- data2 %>% filter(Cluster == 3)
data3 <- data3 %>% filter(Cluster == 3)
my_features <- c(data1$X[1:3], data2$X[1:3], data3$X[1:3]) %>% unique()
dt <- read.csv('Results/HumanUrine/subCluster/Expression/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster3.Expression.SCT.csv', row.names = 1)
dt <- dt %>% t() %>% as.data.frame()
dt.sub <- dt %>% select(all_of(my_features))
dt.sub$Group <- ifelse(str_detect(rownames(dt.sub), 'D'), 'AD', 'Con')
dt.sub <- relocate(dt.sub, Group)
for (i in 2:ncol(dt.sub)) {
  myTitle <- colnames(dt.sub)[i]
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(dt.sub, aes(x = Group, fill = Group)) + theme_bw() +
    geom_jitter(aes_(y = as.name(colnames(dt.sub)[i])), alpha = 0.6, size = 2) + 
    geom_violin(aes_(y = as.name(colnames(dt.sub)[i])), scale = 'width', width = 0.7, alpha = 0.7) +
    geom_boxplot(aes_(y = as.name(colnames(dt.sub)[i])), width = 0.5, size = 0.6) +  
    geom_signif(aes_(y = as.name(colnames(dt.sub)[i])), comparisons = list(c("AD","Con")), 
                map_signif_level = F, step_increase = 0.05, tip_length = 0.01, 
                test = "t.test", textsize = 2.5) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure6/HumanUrine.pseudoPer10.Dim1-20.0.1.Cluster3.',myTitle, '.Pro.png')
  ggsave(p2, filename = fn2, width = 3.5, height = 3.5)
  fn2 <- str_c('Results/Figure6/HumanUrine.pseudoPer10.Dim1-20.0.1.Cluster3.',myTitle, '.Pro..pdf')
  ggsave(p2, filename = fn2, width = 3.5, height = 3.5)
}
# Protein与MMSE和MoCA-B进行Correlation----Tang20230223
dt.sub$ID <- str_remove(rownames(dt.sub), 'u')
clinical <- read.csv('Human_ClinicalData.csv')
data <- left_join(dt.sub, clinical, by = 'ID')
data <- data %>% relocate(Group, ID, MMSE, MoCA.B)
data <- data %>% filter(Group == 'AD')
for (i in 1:5) {
  p <- ggplot(data, aes_(x = as.name(colnames(data)[i+4]))) + theme_bw() +
    geom_point(aes(y = MMSE, color = Group), size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
    stat_smooth(aes(y = MMSE), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
    # stat_poly_eq(aes(y = MMSE, label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x, parse=T,size=5)
    stat_cor(aes(y = MMSE), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = 'darkgreen') +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  fn <- str_c('Results/Figure6/HumanUrine_Cluster3_AD_Correlation_MMSE vs ', colnames(data)[i+4], '.png')
  ggsave(p, filename = fn, width = 4, height = 4)
  fn <- str_c('Results/Figure6/HumanUrine_Cluster3_AD_Correlation_MMSE vs ', colnames(data)[i+4], '.pdf')
  ggsave(p, filename = fn, width = 4, height = 4)
  p <- ggplot(data, aes_(x = as.name(colnames(data)[i+4]))) + theme_bw() +
    geom_point(aes(y = MoCA.B, color = Group), size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
    stat_smooth(aes(y = MoCA.B), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
    # stat_poly_eq(aes(y = MoCA.B, label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x, parse=T,size=5)
    stat_cor(aes(y = MoCA.B), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = 'darkblue') +
    ylab('MoCA-B') + 
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  fn <- str_c('Results/Figure6/HumanUrine_Cluster3_AD_Correlation_MoCA-B vs ', colnames(data)[i+4], '.png')
  ggsave(p, filename = fn, width = 4, height = 4)
  fn <- str_c('Results/Figure6/HumanUrine_Cluster3_AD_Correlation_MoCA-B vs ', colnames(data)[i+4], '.pdf')
  ggsave(p, filename = fn, width = 4, height = 4)
}




# -- ------------------------------------------------------------------------------------------------------------------
# Saliva
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
# Cluster == 14
data1 <- read.csv('Results/HumanSaliva/subCluster/caret/RDS/pseudoPer20.Res0.1.Cluster14_Train.1234_Logistic.csv')
data2 <- read.csv('Results/HumanSaliva/subCluster/caret/RDS/pseudoPer20.Res0.1.Cluster14_Train.5678_Logistic.csv')
data3 <- read.csv('Results/HumanSaliva/subCluster/caret/RDS/pseudoPer20.Res0.1.Cluster14_Train.2022_Logistic.csv')
my_features <- c(data1$X[1:3], data2$X[1:3], data3$X[1:3]) %>% unique()
dt <- read.csv('Results/HumanSaliva/subCluster/Expression/HumanSaliva.pseudoPer20.Merge_SCT_Harmony.Dim1-20.0.1.Cluster14.Expression.SCT.csv', row.names = 1)
dt <- dt %>% t() %>% as.data.frame()
dt.sub <- dt %>% select(all_of(my_features))
dt.sub$Group <- ifelse(str_detect(rownames(dt.sub), 'D'), 'AD', 'Con')
dt.sub <- relocate(dt.sub, Group)
for (i in 2:ncol(dt.sub)) {
  myTitle <- colnames(dt.sub)[i]
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(dt.sub, aes(x = Group, fill = Group)) + theme_bw() +
    geom_boxplot(aes_(y = as.name(colnames(dt.sub)[i])), width = 0.5, size = 0.6) +  
    geom_jitter(aes_(y = as.name(colnames(dt.sub)[i])), width = 0.2, alpha = 0.6, size = 2) + 
    # geom_violin(aes_(y = as.name(colnames(dt.sub)[i])), scale = 'width', width = 0.7, alpha = 0.7) +
    geom_signif(aes_(y = as.name(colnames(dt.sub)[i])), comparisons = list(c("AD","Con")), 
                map_signif_level = F, step_increase = 0.05, tip_length = 0.01, 
                test = "t.test", textsize = 2.5) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure6/HumanSaliva.pseudoPer20.Dim1-20.0.1.Cluster14.',myTitle, '.Pro.png')
  ggsave(p2, filename = fn2, width = 3, height = 4)
  fn2 <- str_c('Results/Figure6/HumanSaliva.pseudoPer20.Dim1-20.0.1.Cluster14.',myTitle, '.Pro..pdf')
  ggsave(p2, filename = fn2, width = 3, height = 4)
}


# -- ------------------------------------------------------------------------------------------------------------------
# Plasma
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
# Cluster == 10
data1 <- read.csv('Results/HumanPlasma/subCluster/caret/RDS/pseudoPer30.Res0.3.Cluster10_Train.1234_Logistic.csv')
data2 <- read.csv('Results/HumanPlasma/subCluster/caret/RDS/pseudoPer30.Res0.3.Cluster10_Train.5678_Logistic.csv')
data3 <- read.csv('Results/HumanPlasma/subCluster/caret/RDS/pseudoPer30.Res0.3.Cluster10_Train.2022_Logistic.csv')
my_features <- c(data1$X[1:3], data2$X[1:3], data3$X[1:3]) %>% unique()
dt <- read.csv('Results/HumanPlasma/subCluster/Expression/HumanPlasma.AllSamples.pseudoPer30.Merge_SCT_Harmony.Dim1-20.0.3.Cluster10.Expression.SCT.csv', row.names = 1)
dt <- dt %>% t() %>% as.data.frame()
dt.sub <- dt %>% select(all_of(my_features))
dt.sub$Group <- ifelse(str_detect(rownames(dt.sub), 'D'), 'AD', 'Con')
dt.sub <- relocate(dt.sub, Group)
for (i in 2:ncol(dt.sub)) {
  myTitle <- colnames(dt.sub)[i]
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(dt.sub, aes(x = Group, fill = Group)) + theme_bw() +
    geom_jitter(aes_(y = as.name(colnames(dt.sub)[i])), alpha = 0.6, size = 2) + 
    geom_violin(aes_(y = as.name(colnames(dt.sub)[i])), scale = 'width', width = 0.7, alpha = 0.7) +
    geom_boxplot(aes_(y = as.name(colnames(dt.sub)[i])), width = 0.5, size = 0.6) +  
    geom_signif(aes_(y = as.name(colnames(dt.sub)[i])), comparisons = list(c("AD","Con")), 
                map_signif_level = F, step_increase = 0.05, tip_length = 0.01, 
                test = "t.test", textsize = 2.5) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure6/HumanPlasma.pseudoPer30.Dim1-20.0.3.Cluster10.',myTitle, '.Pro.png')
  ggsave(p2, filename = fn2, width = 3.5, height = 3.5)
  fn2 <- str_c('Results/Figure6/HumanPlasma.pseudoPer30.Dim1-20.0.3.Cluster10.',myTitle, '.Pro..pdf')
  ggsave(p2, filename = fn2, width = 3.5, height = 3.5)
}


# -- ------------------------------------------------------------------------------------------------------------------
# Tear
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
# Cluster == 7
data1 <- read.csv('Results/HumanTear/subCluster/caret/RDS/pseudoPer10.Res0.3.Cluster7_Train.1234_Logistic.csv')
data2 <- read.csv('Results/HumanTear/subCluster/caret/RDS/pseudoPer10.Res0.3.Cluster7_Train.5678_Logistic.csv')
data3 <- read.csv('Results/HumanTear/subCluster/caret/RDS/pseudoPer10.Res0.3.Cluster10_Train.2022_Logistic.csv')
my_features <- c(data1$X[1:3], data2$X[1:3], data3$X[1:3]) %>% unique()
dt <- read.csv('Results/HumanTear/subCluster/Expression/HumanTear.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.3.Cluster7.Expression.SCT.csv', row.names = 1)
dt <- dt %>% t() %>% as.data.frame()
dt.sub <- dt %>% select(all_of(my_features))
dt.sub$Group <- ifelse(str_detect(rownames(dt.sub), 'D'), 'AD', 'Con')
dt.sub <- relocate(dt.sub, Group)
for (i in 2:ncol(dt.sub)) {
  myTitle <- colnames(dt.sub)[i]
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(dt.sub, aes(x = Group, fill = Group)) + theme_bw() +
    geom_jitter(aes_(y = as.name(colnames(dt.sub)[i])), alpha = 0.6, size = 2) + 
    geom_violin(aes_(y = as.name(colnames(dt.sub)[i])), scale = 'width', width = 0.7, alpha = 0.7) +
    geom_boxplot(aes_(y = as.name(colnames(dt.sub)[i])), width = 0.5, size = 0.6) +  
    geom_signif(aes_(y = as.name(colnames(dt.sub)[i])), comparisons = list(c("AD","Con")), 
                map_signif_level = F, step_increase = 0.05, tip_length = 0.01, 
                test = "t.test", textsize = 2.5) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure6/HumanTear.pseudoPer10.Dim1-20.0.3.Cluster7.',myTitle, '.Pro.png')
  ggsave(p2, filename = fn2, width = 3.5, height = 3.5)
  fn2 <- str_c('Results/Figure6/HumanTear.pseudoPer10.Dim1-20.0.3.Cluster7.',myTitle, '.Pro..pdf')
  ggsave(p2, filename = fn2, width = 3.5, height = 3.5)
}


