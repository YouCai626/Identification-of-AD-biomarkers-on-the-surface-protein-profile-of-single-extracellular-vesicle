getwd()
# 设置最大的内存
if(.Platform$OS.type == "windows") withAutoprint({
  memory.size()
  memory.size(TRUE)
  memory.limit()
})
# options(future.globals.maxSize = 30*1024^3) #设置为30G
# 可视化常用包
library(tidyverse) # 包含了ggplot2等包
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
# 1.对人的样本原始数据进行处理, 从原始数据获取Raw表达矩阵
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
sampleID <- read.csv('SampleID_SampleName.csv')
head(sampleID)
sampleID <- sampleID %>% filter(Species == 'Human') %>% filter(Batch == 1) %>% 
  filter(!SampleName %in% c('DMp_3', 'DMu_3', 'DMs_3', 'DMt_3'))
sampleID <- sampleID %>% arrange(Source, Diseases, Gender, SampleName)
write_csv(sampleID, 'Human_SampleID_SampleName.csv')
table(sampleID$Source,sampleID$Diseases)

rm(list=ls())
gc()
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID$EVcounts <- c()
setwd("D:/2.AD.PBA.CAM.Article.V1.Proj")
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('RawData/',sampleID[i,1],'.total_ev_protein.csv')
  data <- read.csv(fn1, header = T)
  wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
  # 获得每个样本的EV总数
  sampleID$EVcounts[i] <- nrow(wide_data)
  rownames(wide_data) <- wide_data$ev
  wide_data <- wide_data[, -c(1:2)]
  pro.exp <- as.data.frame(colSums(wide_data))
  colnames(pro.exp) <- sampleID[i,2]
  pro.exp$Pro <- rownames(pro.exp)
  fn2 <- str_c('pseudoData/Raw/', sampleID[i,1], '.Raw.csv')
  write_csv(pro.exp, fn2)
  rm(data, wide_data)
}
head(sampleID)
setwd("D:/2.AD.PBA.CAM.Article.V3.Proj")
write_csv(sampleID, 'Human_SampleID_SampleName EV counts.csv')

# 由于之前已经跑过这部分的内容了，所以直接把之前的结果COPY过来
# R语言中文件夹中文件COPY
rm(list=ls())
gc()
setwd("D:/2.AD.PBA.CAM.Article.V3.Proj")
sampleID <- read.csv('Human_SampleID_SampleName.csv')
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('D:/2.AD.PBA.CAM.Article.V1.Proj/RawData/',sampleID[i,1],'.total_ev_protein.csv')
  fn2 <- str_c('D:/2.AD.PBA.CAM.Article.V3.Proj/RawData/',sampleID[i,1],'.total_ev_protein.csv')
  file.copy(from = fn1, to = fn2)
}
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('D:/2.AD.PBA.CAM.Article.V1.Proj/pseudoData/Raw/',sampleID[i,1],'.Raw.csv')
  fn2 <- str_c('D:/2.AD.PBA.CAM.Article.V3.Proj/pseudoData/Raw/',sampleID[i,1],'.Raw.csv')
  file.copy(from = fn1, to = fn2)
}
EV.Counts <- read.csv('D:/2.AD.PBA.CAM.Article.V1.Proj/SampleID_SampleName EV counts.csv')
Human.EV.Counts <- left_join(sampleID, EV.Counts, by = 'SampleID')
Human.EV.Counts <- Human.EV.Counts %>% select(!SampleName.y:Batch.y)
Human.EV.Counts <- Human.EV.Counts %>% select(!starts_with('pseudoSize')) %>% select(!starts_with('sEV'))
names(Human.EV.Counts) <- names(Human.EV.Counts) %>% str_remove('.x')
write_csv(Human.EV.Counts, 'Human_SampleID_SampleName EV counts.csv')
# 将人的临床数据copy过来
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/ClinicalData.csv',
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/Human_ClinicalData.csv')

# -- ------------------------------------------------------------------------------------------------------------------
# 2.EV数进行统计展示
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("D:/2.AD.PBA.CAM.Article.V3.Proj")
data <- read.csv('Human_SampleID_SampleName EV counts.csv')
data$Diseases <- data$Diseases %>% factor(levels = c('AD', 'Con'))
data$Source <- data$Source %>% factor(levels = c('Plasma', 'Saliva', 'Tear', 'Urine'))
p <- ggplot(data, aes(x = Diseases)) + theme_bw() +
  # geom_violin(aes(y = EVcounts, fill = Diseases), scale = 'width', width = 0.8, alpha = 0.8) +
  geom_boxplot(aes(y = EVcounts, fill = Diseases), width = 0.8, alpha = 1, size = 0.8) +
  geom_jitter(aes(y = EVcounts, colour = Source), width = 0.3, size = 2, alpha = 0.8) +  
  geom_signif(aes(y = EVcounts), 
              comparisons = list(c('AD', 'Con')), 
              map_signif_level = F, 
              tip_length = 0.02, test = "t.test", textsize = 3) +
  xlab('') + ylab('') + ggtitle('EVcounts') +
  scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
  scale_color_manual(values = c('#D62728E5','#9467BDE5','#1F77B4E5','#BCBD22E5')) +
  theme(legend.position = "none") + facet_wrap(~Source, scales = 'free_y', nrow = 1) +
  theme(strip.text.x = element_text(size = 14)) # 设置分面的字字体大小、颜色、背景、边框
ggsave(p, filename = 'Results/Figure4/Human_EVcounts.png', width = 12, height = 4)
ggsave(p, filename = 'Results/Figure4/Human_EVcounts.pdf', width = 12, height = 4)

# -- ------------------------------------------------------------------------------------------------------------------
# 3.Human样本的临床信息可视化：Age；Gender；MMSE；MoCA-B
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Human_ClinicalData.csv')
data$Group <- ifelse(str_detect(data$ID, 'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
p <- ggplot(data, aes(x = Group)) + theme_bw() +
  # geom_violin(aes(y = Age, fill = Group), scale = 'width', width = 0.7, alpha = 0.7) +
  geom_boxplot(aes(y = Age, fill = Group), width = 0.8, alpha = 0.8, size = 0.8) +  
  geom_jitter(aes(y = Age), width = 0.4, size = 2, alpha = 0.9) + 
  geom_signif(aes(y = Age), 
              comparisons = list(c('AD', 'Con')), 
              map_signif_level = F, step_increase = 0.06, 
              tip_length = 0.01, test = "t.test", textsize = 2) +
  xlab('') + ylab('') + ggtitle('Age') +
  scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
  theme(legend.position = "none")
ggsave(p, filename = 'Results/Figure4/Human_Clinical_Age.png', width = 4, height = 4)
ggsave(p, filename = 'Results/Figure4/Human_Clinical_Age.pdf', width = 4, height = 4)

p <- ggplot(data, aes(x = Group)) + theme_bw() +
  geom_bar(aes(fill = Group), width = 0.8, alpha = 0.8, size = 0.8) + 
  xlab('') + ylab('') + ggtitle('Gender') +
  scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) + facet_grid(~Gender) +
  theme(strip.text.x = element_text(size = 14)) + # 设置分面的字字体大小、颜色、背景、边框，
  scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
  theme(legend.position = "none")
ggsave(p, filename = 'Results/Figure4/Human_Clinical_Gender.png', width = 4, height = 4)
ggsave(p, filename = 'Results/Figure4/Human_Clinical_Gender.pdf', width = 4, height = 4)

p <- ggplot(data, aes(x = Group)) + theme_bw() +
  geom_violin(aes(y = MMSE, fill = Group), scale = 'width', width = 0.8, alpha = 0.8) +
  # geom_boxplot(aes(y = MMSE, fill = Group), width = 0.8, alpha = 0.8, size = 0.8) +
  geom_jitter(aes(y = MMSE), width = 0.4, size = 2, alpha = 0.9) + 
  geom_signif(aes(y = MMSE), 
              comparisons = list(c('AD', 'Con')), 
              map_signif_level = F, step_increase = 0.06, 
              tip_length = 0.01, test = "wilcox.test", textsize = 2) +
  geom_hline(yintercept = 24, lty = 2, lwd = 0.6, col = "darkgreen") +
  xlab('') + ylab('') + ggtitle('MMSE') +
  scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
  theme(legend.position = "none")
ggsave(p, filename = 'Results/Figure4/Human_Clinical_MMSE.png', width = 4, height = 4)
ggsave(p, filename = 'Results/Figure4/Human_Clinical_MMSE.pdf', width = 4, height = 4)

p <- ggplot(data, aes(x = Group)) + theme_bw() +
  geom_violin(aes(y = MoCA.B, fill = Group), scale = 'width', width = 0.8, alpha = 0.8) +
  # geom_boxplot(aes(y = MMSE, fill = Group), width = 0.8, alpha = 0.8, size = 0.8) +
  geom_jitter(aes(y = MoCA.B), width = 0.4, size = 2, alpha = 0.9) + 
  geom_signif(aes(y = MoCA.B), 
              comparisons = list(c('AD', 'Con')), 
              map_signif_level = F, step_increase = 0.06, 
              tip_length = 0.01, test = "wilcox.test", textsize = 2) +
  geom_hline(yintercept = 22, lty = 2, lwd = 0.6, col = "darkgreen") +
  xlab('') + ylab('') + ggtitle('MoCA-B') +
  scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
  theme(legend.position = "none")
ggsave(p, filename = 'Results/Figure4/Human_Clinical_MoCA-B.png', width = 4, height = 4)
ggsave(p, filename = 'Results/Figure4/Human_Clinical_MoCA-B.pdf', width = 4, height = 4)

# -- ------------------------------------------------------------------------------------------------------------------
# 4.将每个样本的蛋白表达总量构建成表达矩阵
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
sampleID <- read.csv('Human_SampleID_SampleName.csv')
fn <- str_c('pseudoData/Raw/', sampleID[1,1], '.Raw.csv')
Raw.Counts <- read_csv(fn, show_col_types = F)
for (i in 2:nrow(sampleID)) {
  fn <- str_c('pseudoData/Raw/', sampleID[i,1], '.Raw.csv')
  data <- read_csv(fn,show_col_types = F)
  Raw.Counts <- full_join(Raw.Counts, data, by = 'Pro', fill = 0) # 取并集
}
Raw.Counts <- relocate(Raw.Counts, Pro)
Raw.Counts[is.na(Raw.Counts)] <- 0  # 将NA替换成0
write_csv(Raw.Counts, 'Results/Human_Raw.Counts.csv')
# 将每一种样本的单独拧出来，因为每一种样本要单独分析
rm(list=ls())
gc()
data <- read.csv('Results/Human_Raw.Counts.csv', row.names = 1)
sampleID <- read_csv('Human_SampleID_SampleName.csv', show_col_types = F)
# Plasma
sub.ID <- sampleID %>% filter(Source == 'Plasma')
sub.data <- data %>% dplyr::select(sub.ID$SampleName)
sub.data$mean <- rowMeans(sub.data)
sub.data <- sub.data %>% filter(mean > 0) %>% dplyr::select(!mean) # 将表达量为0的Pro去除
write_csv(sub.ID, 'Results/HumanPlasma_Samples.Annotation.csv')
write.csv(sub.data, 'Results/HumanPlasma_Raw.Counts.csv')
# Urine
sub.ID <- sampleID %>% filter(Source == 'Urine')
sub.data <- data %>% dplyr::select(sub.ID$SampleName)
sub.data$mean <- rowMeans(sub.data)
sub.data <- sub.data %>% filter(mean > 0) %>% dplyr::select(!mean) # 将表达量为0的Pro去除
write_csv(sub.ID, 'Results/HumanUrine_Samples.Annotation.csv')
write.csv(sub.data, 'Results/HumanUrine_Raw.Counts.csv')
# Saliva
sub.ID <- sampleID %>% filter(Source == 'Saliva')
sub.data <- data %>% dplyr::select(sub.ID$SampleName)
sub.data$mean <- rowMeans(sub.data)
sub.data <- sub.data %>% filter(mean > 0) %>% dplyr::select(!mean) # 将表达量为0的Pro去除
write_csv(sub.ID, 'Results/HumanSaliva_Samples.Annotation.csv')
write.csv(sub.data, 'Results/HumanSaliva_Raw.Counts.csv')
# Tear
sub.ID <- sampleID %>% filter(Source == 'Tear')
sub.data <- data %>% dplyr::select(sub.ID$SampleName)
sub.data$mean <- rowMeans(sub.data)
sub.data <- sub.data %>% filter(mean > 0) %>% dplyr::select(!mean) # 将表达量为0的Pro去除
write_csv(sub.ID, 'Results/HumanTear_Samples.Annotation.csv')
write.csv(sub.data, 'Results/HumanTear_Raw.Counts.csv')

# -- ------------------------------------------------------------------------------------------------------------------
# 5.对所有种类的样本进行热图展示，查看是否有组织特异性
# -- ------------------------------------------------------------------------------------------------------------------
# 对表达值进行标准化
library(edgeR) # 用于cpm标准化
library(preprocessCore) # CPM.quantiles函数进行标准化
library(pheatmap)
data <- read.csv('Results/Human_Raw.Counts.csv', row.names = 1)
data$mean <- rowMeans(data)
data <- data %>% filter(mean > 0) %>% dplyr::select(!mean) # 将表达量为0的Pro去除
data.nor <- data %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(data.nor) <- rownames(data)
names(data.nor) <- names(data)
write.csv(data.nor, 'Results/Human_Raw.Counts_cpm.quantiles.csv')
coldata <- names(data.nor) %>% as.data.frame()
names(coldata) <- 'Samples'
coldata$Group <- ifelse(str_detect(coldata$Samples,'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
coldata$Type <- ifelse(str_detect(coldata$Samples,'p'), 'Plasma',
                       ifelse(str_detect(coldata$Samples,'u'), 'Urine',
                              ifelse(str_detect(coldata$Samples,'s'), 'Saliva', 'Tear'))) %>% 
  factor(levels = c('Plasma', 'Saliva', 'Tear', 'Urine'))
rownames(coldata) <- coldata$Samples
coldata <- coldata[-1]
ann_colors = list(Group = c('AD' ='#CC0000FF', 'Con' ='#6699FFFF'),
                  Type = c('Plasma' = '#D62728E5', 'Saliva' = '#9467BDE5',
                           'Tear' = '#1F77B4E5', 'Urine' = '#BCBD22E5')) #注意ann_colors是列表
pheatmap(data.nor, scale = 'row', color = pal_gsea(palette = c("default"), n = 20)(20), cellwidth = 3, cellheight = 1.5,
         display_numbers = F, number_color = 'black', number_format = '%.3f', 
         fontsize_number = 4, show_rownames = F, fontsize_row = 4, fontsize_col = 3, border_color = NA,
         cluster_rows = T, cluster_cols = F, treeheight_row = 0,
         annotation_col = coldata, annotation_colors = ann_colors,
         cutree_rows = 4, gaps_col = c(38,75,112),
         filename = 'Results/Figure4/Human_Raw.Counts_CPM.quantiles_pheatmap.png') # 保存，自动调整纸张大小

# -- ------------------------------------------------------------------------------------------------------------------
# 6.ROC分析+多个Pro的联合ROC
# -- ------------------------------------------------------------------------------------------------------------------
library(edgeR) # 用于cpm标准化
library(preprocessCore) # CPM.quantiles函数进行标准化
library(pROC) # ROC分析与作图
# -- ------------------------------------------------------------------------------------------------------------------
# Plasma
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/HumanPlasma_Raw.Counts.csv', row.names = 1)
data.nor <- data %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(data.nor) <- rownames(data)
names(data.nor) <- names(data)
write.csv(data.nor, 'Results/HumanPlasma_Raw.Counts_cpm.quantiles.csv')
roc.data <- t(data.nor) %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
roc.data <- relocate(roc.data, Group)
auc <- c()
z <- c()
pvalue <- c()
for (j in 1:(ncol(roc.data)-1)){
  rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('AD', 'Con'), direction = '>', data = roc.data)
  auc[j] <- rocPlot$auc
  se <- sqrt(pROC::var(rocPlot))
  b <- auc[j] - 0.5
  z[j] <- (b / se)
  pvalue[j] <- 2 * pt(-abs(z[j]), df=Inf)
  # 对AUC > 0.7 或 < 0.3 的蛋白进行绘图
  if (pvalue[j] < 0.05 & auc[j] >= 0.7) {
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    fn6 <- paste('Results/Figure4/HumanPlasma_Raw.Counts_cpm.quantiles_', proteinname, '.png', sep = '')
    png(fn6, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure4/HumanPlasma_Raw.Counts_cpm.quantiles_', proteinname, '.pdf', sep = '')
    pdf(fn6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
  }
  if (pvalue[j] < 0.05 & auc[j] <= 0.3) {
    rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('AD', 'Con'), direction = '<', data = roc.data)
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    fn6 <- paste('Results/Figure4/HumanPlasma_Raw.Counts_cpm.quantiles_', proteinname, '.png', sep = '')
    png(fn6, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "yellow",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure4/HumanPlasma_Raw.Counts_cpm.quantiles_', proteinname, '.pdf', sep = '')
    pdf(fn6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "yellow",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
  }
}
auc_pvalue <- data.frame(auc, z, pvalue)
rownames(auc_pvalue) <- colnames(roc.data)[-1]
auc_pvalue <- arrange(auc_pvalue, desc(auc))
write.csv(auc_pvalue, file = 'Results/Figure4/Files/HumanPlasma_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.csv')
auc_pvalue <- filter(auc_pvalue, pvalue < 0.05)
write.csv(auc_pvalue, file = 'Results/Figure4/Files/HumanPlasma_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv')
# 将ROC阳性的Protein与MMSE和MoCA-B进行Correlation
rm(list=ls())
gc()
data.nor <- read.csv('Results/HumanPlasma_Raw.Counts_cpm.quantiles.csv', row.names = 1)
roc.data <- t(data.nor) %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
roc.data <- relocate(roc.data, Group)
roc.data$ID <- str_remove(rownames(roc.data), 'p')
roc.data$Samples <- rownames(roc.data)
clinical <- read.csv('Human_ClinicalData.csv')
data <- left_join(roc.data, clinical, by = 'ID')
data.pos <- read.csv('Results/Figure4/Files/HumanPlasma_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv', row.names = 1)
myPro <- c(rownames(data.pos))
dt <- select(data, Group, Samples, MMSE, MoCA.B, all_of(myPro))
names(dt) <- names(dt) %>% str_replace_all('/', '.') # 文件名中有/的时候会报错
i = length(myPro) # 为NECTIN3添加最佳阈值线
for (i in 1:length(myPro)) {
  p <- ggplot(dt, aes_(x = as.name(colnames(dt)[i+4]))) + theme_bw() +
    geom_point(aes(y = MMSE, color = Group), size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
    stat_smooth(aes(y = MMSE), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
    # stat_poly_eq(aes(y = MMSE, label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x, parse=T,size=5)
    stat_cor(aes(y = MMSE), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = 'darkgreen') +
    # geom_vline(xintercept = 6339.312, color = "darkgreen", linetype = 2, size = 0.5) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  fn <- str_c('Results/Figure4/HumanPlasma_Correlation_MMSE vs ', colnames(dt)[i+4], '.png')
  ggsave(p, filename = fn, width = 4, height = 4)
  fn <- str_c('Results/Figure4/HumanPlasma_Correlation_MMSE vs ', colnames(dt)[i+4], '.pdf')
  ggsave(p, filename = fn, width = 4, height = 4)
  p <- ggplot(dt, aes_(x = as.name(colnames(dt)[i+4]))) + theme_bw() +
    geom_point(aes(y = MoCA.B, color = Group), size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
    stat_smooth(aes(y = MoCA.B), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
    # stat_poly_eq(aes(y = MoCA.B, label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x, parse=T,size=5)
    stat_cor(aes(y = MoCA.B), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = 'darkblue') +
    ylab('MoCA-B') + 
    # geom_vline(xintercept = 6339.312, color = "darkgreen", linetype = 2, size = 0.5) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  fn <- str_c('Results/Figure4/HumanPlasma_Correlation_MoCA-B vs ', colnames(dt)[i+4], '.png')
  ggsave(p, filename = fn, width = 4, height = 4)
  fn <- str_c('Results/Figure4/HumanPlasma_Correlation_MoCA-B vs ', colnames(dt)[i+4], '.pdf')
  ggsave(p, filename = fn, width = 4, height = 4)
}
# 将ROC阳性的Protein的表达值图示出来
i = ncol(dt) # 给NECTIN3添加阈值线
for (i in 5:ncol(dt)) {
  myTitle <- gsub("/", "_", colnames(dt)[i])
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(dt, aes(x = Group, fill = Group)) + theme_bw() +
    geom_boxplot(aes_(y = as.name(colnames(dt)[i])), width = 0.6, size = 0.8) +
    geom_violin(aes_(y = as.name(colnames(dt)[i])), scale = 'width', width = 0.8, alpha = 0.6) +
    geom_jitter(aes_(y = as.name(colnames(dt)[i])), width = 0.3, alpha = 0.8, size = 2) + 
    geom_signif(aes_(y = as.name(colnames(dt)[i])), 
                comparisons = list(c('AD', 'Con')), 
                map_signif_level = F, step_increase = 0.01, 
                tip_length = 0.01, test = "wilcox.test", textsize = 2) +
    # geom_hline(yintercept = 6339.312, color = "darkgreen", linetype = 2, size = 0.5) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) + 
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure4/HumanPlasma_Raw.Counts_CPM.quantiles_wilcox.test_', myTitle, '.png')
  ggsave(p2, filename = fn2, width = 4, height = 4)
  fn2 <- str_c('Results/Figure4/HumanPlasma_Raw.Counts_CPM.quantiles_wilcox.test_', myTitle, '.pdf')
  ggsave(p2, filename = fn2, width = 4, height = 4)
}
# 2-4个Pro的联合ROC
data.pos$auc <- ifelse(data.pos$auc > 0.5, data.pos$auc, 1-data.pos$auc)
data.pos <- data.pos %>% arrange(desc(auc))
roc.data <- roc.data %>% select(Group, rownames(data.pos)[1:5])
f2 <- glm(Group ~ roc.data[,2] + roc.data[,3], data = roc.data, family = binomial())
roc.data$Pro2 <- predict(f2, newdata = roc.data, type = "response")
f3 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4], data = roc.data, family = binomial())
roc.data$Pro3 <- predict(f3, newdata = roc.data, type = "response")
f4 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4] + roc.data[,5], data = roc.data, family = binomial())
roc.data$Pro4 <- predict(f4, newdata = roc.data, type = "response")
f5 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4] + roc.data[,5] + roc.data[,6], data = roc.data, family = binomial())
roc.data$Pro5 <- predict(f5, newdata = roc.data, type = "response")
roc.pro2 <- roc(roc.data$Group, roc.data$Pro2)
roc.pro2$auc
roc.pro3 <- roc(roc.data$Group, roc.data$Pro3)
roc.pro3$auc
roc.pro4 <- roc(roc.data$Group, roc.data$Pro4)
roc.pro4$auc
roc.pro5 <- roc(roc.data$Group, roc.data$Pro5)
roc.pro5$auc
show_col(pal_npg(alpha = 0.8)(10))
png('Results/Figure4/HumanPlasma_Raw.Counts_CPM.quantiles_ROC_2-5.png', height = 400, width = 400)
plot.roc(roc.pro2, col = "#3C5488CC", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#3C5488CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.35,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro3, col = "#F39B7FCC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F39B7FCC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.3,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro4, col = "#4DBBD5CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#4DBBD5CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.25,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro5, col = "#DC0000CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#DC0000CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.2,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
pdf('Results/Figure4/HumanPlasma_Raw.Counts_CPM.quantiles_ROC_2-5.pdf', height = 8, width = 8)
plot.roc(roc.pro2, col = "#3C5488CC", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#3C5488CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.35,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro3, col = "#F39B7FCC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F39B7FCC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.3,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro4, col = "#4DBBD5CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#4DBBD5CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.25,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro5, col = "#DC0000CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#DC0000CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.2,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
# -- ------------------------------------------------------------------------------------------------------------------
# Saliva
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/HumanSaliva_Raw.Counts.csv', row.names = 1)
data.nor <- data %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(data.nor) <- rownames(data)
names(data.nor) <- names(data)
write.csv(data.nor, 'Results/HumanSaliva_Raw.Counts_cpm.quantiles.csv')
roc.data <- t(data.nor) %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
roc.data <- relocate(roc.data, Group)
auc <- c()
z <- c()
pvalue <- c()
for (j in 1:(ncol(roc.data)-1)){
  rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('AD', 'Con'), direction = '>', data = roc.data)
  auc[j] <- rocPlot$auc
  se <- sqrt(pROC::var(rocPlot))
  b <- auc[j] - 0.5
  z[j] <- (b / se)
  pvalue[j] <- 2 * pt(-abs(z[j]), df=Inf)
  # 对AUC > 0.7 或 < 0.3 的蛋白进行绘图
  if (pvalue[j] < 0.05 & auc[j] >= 0.7) {
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    fn6 <- paste('Results/Figure4/HumanSaliva_Raw.Counts_cpm.quantiles_', proteinname, '.png', sep = '')
    png(fn6, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure4/HumanSaliva_Raw.Counts_cpm.quantiles_', proteinname, '.pdf', sep = '')
    pdf(fn6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
  }
  if (pvalue[j] < 0.05 & auc[j] <= 0.3) {
    rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('AD', 'Con'), direction = '<', data = roc.data)
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    fn6 <- paste('Results/Figure4/HumanSaliva_Raw.Counts_cpm.quantiles_', proteinname, '.png', sep = '')
    png(fn6, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "yellow",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure4/HumanSaliva_Raw.Counts_cpm.quantiles_', proteinname, '.pdf', sep = '')
    pdf(fn6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "yellow",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
  }
}
auc_pvalue <- data.frame(auc, z, pvalue)
rownames(auc_pvalue) <- colnames(roc.data)[-1]
auc_pvalue <- arrange(auc_pvalue, desc(auc))
write.csv(auc_pvalue, file = 'Results/Figure4/Files/HumanSaliva_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.csv')
auc_pvalue <- filter(auc_pvalue, pvalue < 0.05)
write.csv(auc_pvalue, file = 'Results/Figure4/Files/HumanSaliva_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv')
# 将ROC阳性的Protein与MMSE和MoCA-B进行Correlation
rm(list=ls())
gc()
data.nor <- read.csv('Results/HumanSaliva_Raw.Counts_cpm.quantiles.csv', row.names = 1)
roc.data <- t(data.nor) %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
roc.data <- relocate(roc.data, Group)
roc.data$ID <- str_remove(rownames(roc.data), 's')
roc.data$Samples <- rownames(roc.data)
clinical <- read.csv('Human_ClinicalData.csv')
data <- left_join(roc.data, clinical, by = 'ID')
data.pos <- read.csv('Results/Figure4/Files/HumanSaliva_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv', row.names = 1)
myPro <- rownames(data.pos)
dt <- select(data, Group, Samples, MMSE, MoCA.B, all_of(myPro))
names(dt) <- names(dt) %>% str_replace_all('/', '.') # 文件名中有/的时候会报错
i = 1 # THY1
for (i in 1:length(myPro)) {
  p <- ggplot(dt, aes_(x = as.name(colnames(dt)[i+4]))) + theme_bw() +
    geom_point(aes(y = MMSE, color = Group), size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
    stat_smooth(aes(y = MMSE), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
    # stat_poly_eq(aes(y = MMSE, label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x, parse=T,size=5)
    stat_cor(aes(y = MMSE), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = 'darkgreen') +
    # geom_vline(xintercept = 1662.642, color = "darkgreen", linetype = 2, size = 0.5) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  fn <- str_c('Results/Figure4/HumanSaliva_Correlation_MMSE vs ', colnames(dt)[i+4], '.png')
  ggsave(p, filename = fn, width = 4, height = 4)
  fn <- str_c('Results/Figure4/HumanSaliva_Correlation_MMSE vs ', colnames(dt)[i+4], '.pdf')
  ggsave(p, filename = fn, width = 4, height = 4)
  p <- ggplot(dt, aes_(x = as.name(colnames(dt)[i+4]))) + theme_bw() +
    geom_point(aes(y = MoCA.B, color = Group), size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
    stat_smooth(aes(y = MoCA.B), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
    # stat_poly_eq(aes(y = MoCA.B, label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x, parse=T,size=5)
    stat_cor(aes(y = MoCA.B), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = 'darkblue') +
    ylab('MoCA-B') + 
    # geom_vline(xintercept = 1662.642, color = "darkgreen", linetype = 2, size = 0.5) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  fn <- str_c('Results/Figure4/HumanSaliva_Correlation_MoCA-B vs ', colnames(dt)[i+4], '.png')
  ggsave(p, filename = fn, width = 4, height = 4)
  fn <- str_c('Results/Figure4/HumanSaliva_Correlation_MoCA-B vs ', colnames(dt)[i+4], '.pdf')
  ggsave(p, filename = fn, width = 4, height = 4)
}
# 将ROC阳性的Protein的表达值图示出来
i = 5
for (i in 5:ncol(dt)) {
  myTitle <- gsub("/", "_", colnames(dt)[i])
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(dt, aes(x = Group, fill = Group)) + theme_bw() +
    geom_boxplot(aes_(y = as.name(colnames(dt)[i])), width = 0.6, size = 0.8) +  
    geom_violin(aes_(y = as.name(colnames(dt)[i])), scale = 'width', width = 0.8, alpha = 0.6) +
    geom_jitter(aes_(y = as.name(colnames(dt)[i])), width = 0.3, alpha = 0.8, size = 2) + 
    geom_signif(aes_(y = as.name(colnames(dt)[i])), 
                comparisons = list(c('AD', 'Con')), 
                map_signif_level = F, step_increase = 0.01, 
                tip_length = 0.01, test = "wilcox.test", textsize = 2) +
    # geom_hline(yintercept = 1662.642, color = "darkgreen", linetype = 2, size = 0.5) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure4/HumanSaliva_Raw.Counts_CPM.quantiles_wilcox.test_', myTitle, '.png')
  ggsave(p2, filename = fn2, width = 4, height = 4)
  fn2 <- str_c('Results/Figure4/HumanSaliva_Raw.Counts_CPM.quantiles_wilcox.test_', myTitle, '.pdf')
  ggsave(p2, filename = fn2, width = 4, height = 4)
}
# 2-4个Pro的联合ROC
data.pos$auc <- ifelse(data.pos$auc > 0.5, data.pos$auc, 1-data.pos$auc)
data.pos <- data.pos %>% arrange(desc(auc))
roc.data <- roc.data %>% select(Group, rownames(data.pos)[1:5])
f2 <- glm(Group ~ roc.data[,2] + roc.data[,3], data = roc.data, family = binomial())
roc.data$Pro2 <- predict(f2, newdata = roc.data, type = "response")
f3 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4], data = roc.data, family = binomial())
roc.data$Pro3 <- predict(f3, newdata = roc.data, type = "response")
f4 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4] + roc.data[,5], data = roc.data, family = binomial())
roc.data$Pro4 <- predict(f4, newdata = roc.data, type = "response")
f5 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4] + roc.data[,5] + roc.data[,6], data = roc.data, family = binomial())
roc.data$Pro5 <- predict(f5, newdata = roc.data, type = "response")
roc.pro2 <- roc(roc.data$Group, roc.data$Pro2)
roc.pro2$auc
roc.pro3 <- roc(roc.data$Group, roc.data$Pro3)
roc.pro3$auc
roc.pro4 <- roc(roc.data$Group, roc.data$Pro4)
roc.pro4$auc
roc.pro5 <- roc(roc.data$Group, roc.data$Pro5)
roc.pro5$auc
show_col(pal_npg(alpha = 0.8)(10))
png('Results/Figure4/HumanSaliva_Raw.Counts_CPM.quantiles_ROC_2-5.png', height = 400, width = 400)
plot.roc(roc.pro2, col = "#3C5488CC", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#3C5488CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.35,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro3, col = "#F39B7FCC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F39B7FCC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.3,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro4, col = "#4DBBD5CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#4DBBD5CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.25,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro5, col = "#DC0000CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#DC0000CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.2,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
pdf('Results/Figure4/HumanSaliva_Raw.Counts_CPM.quantiles_ROC_2-5.pdf', height = 8, width = 8)
plot.roc(roc.pro2, col = "#3C5488CC", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#3C5488CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.35,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro3, col = "#F39B7FCC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F39B7FCC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.3,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro4, col = "#4DBBD5CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#4DBBD5CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.25,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro5, col = "#DC0000CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#DC0000CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.2,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
# -- ------------------------------------------------------------------------------------------------------------------
# Tear
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/HumanTear_Raw.Counts.csv', row.names = 1)
data.nor <- data %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(data.nor) <- rownames(data)
names(data.nor) <- names(data)
write.csv(data.nor, 'Results/HumanTear_Raw.Counts_cpm.quantiles.csv')
roc.data <- t(data.nor) %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
roc.data <- relocate(roc.data, Group)
auc <- c()
z <- c()
pvalue <- c()
for (j in 1:(ncol(roc.data)-1)){
  rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('AD', 'Con'), direction = '>', data = roc.data)
  auc[j] <- rocPlot$auc
  se <- sqrt(pROC::var(rocPlot))
  b <- auc[j] - 0.5
  z[j] <- (b / se)
  pvalue[j] <- 2 * pt(-abs(z[j]), df=Inf)
  # 对AUC > 0.7 或 < 0.3 的蛋白进行绘图
  if (pvalue[j] < 0.05 & auc[j] >= 0.7) {
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    fn6 <- paste('Results/Figure4/HumanTear_Raw.Counts_cpm.quantiles_', proteinname, '.png', sep = '')
    png(fn6, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure4/HumanTear_Raw.Counts_cpm.quantiles_', proteinname, '.pdf', sep = '')
    pdf(fn6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
  }
  if (pvalue[j] < 0.05 & auc[j] <= 0.3) {
    rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('AD', 'Con'), direction = '<', data = roc.data)
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    fn6 <- paste('Results/Figure4/HumanTear_Raw.Counts_cpm.quantiles_', proteinname, '.png', sep = '')
    png(fn6, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "yellow",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure4/HumanTear_Raw.Counts_cpm.quantiles_', proteinname, '.pdf', sep = '')
    pdf(fn6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "yellow",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
  }
}
auc_pvalue <- data.frame(auc, z, pvalue)
rownames(auc_pvalue) <- colnames(roc.data)[-1]
auc_pvalue <- arrange(auc_pvalue, desc(auc))
write.csv(auc_pvalue, file = 'Results/Figure4/Files/HumanTear_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.csv')
auc_pvalue <- filter(auc_pvalue, pvalue < 0.05)
write.csv(auc_pvalue, file = 'Results/Figure4/Files/HumanTear_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv')
# 将ROC阳性的Protein与MMSE和MoCA-B进行Correlation
rm(list=ls())
gc()
data.nor <- read.csv('Results/HumanTear_Raw.Counts_cpm.quantiles.csv', row.names = 1)
roc.data <- t(data.nor) %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
roc.data <- relocate(roc.data, Group)
roc.data$ID <- str_remove(rownames(roc.data), 't')
roc.data$Samples <- rownames(roc.data)
clinical <- read.csv('Human_ClinicalData.csv')
data <- left_join(roc.data, clinical, by = 'ID')
data.pos <- read.csv('Results/Figure4/Files/HumanTear_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv', row.names = 1)
myPro <- rownames(data.pos)
dt <- select(data, Group, Samples, MMSE, MoCA.B, all_of(myPro))
names(dt) <- names(dt) %>% str_replace_all('/', '.') # 文件名中有/的时候会报错
i = 1
for (i in 1:length(myPro)) {
  p <- ggplot(dt, aes_(x = as.name(colnames(dt)[i+4]))) + theme_bw() +
    geom_point(aes(y = MMSE, color = Group), size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
    stat_smooth(aes(y = MMSE), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
    # stat_poly_eq(aes(y = MMSE, label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x, parse=T,size=5)
    stat_cor(aes(y = MMSE), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = 'darkgreen') +
    geom_vline(xintercept = 4742.664, color = "darkgreen", linetype = 2, size = 0.5) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  fn <- str_c('Results/Figure4/HumanTear_Correlation_MMSE vs ', colnames(dt)[i+4], '.png')
  ggsave(p, filename = fn, width = 4, height = 4)
  fn <- str_c('Results/Figure4/HumanTear_Correlation_MMSE vs ', colnames(dt)[i+4], '.pdf')
  ggsave(p, filename = fn, width = 4, height = 4)
  p <- ggplot(dt, aes_(x = as.name(colnames(dt)[i+4]))) + theme_bw() +
    geom_point(aes(y = MoCA.B, color = Group), size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
    stat_smooth(aes(y = MoCA.B), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
    # stat_poly_eq(aes(y = MoCA.B, label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x, parse=T,size=5)
    stat_cor(aes(y = MoCA.B), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = 'darkblue') +
    ylab('MoCA-B') + 
    geom_vline(xintercept = 4742.664, color = "darkgreen", linetype = 2, size = 0.5) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  fn <- str_c('Results/Figure4/HumanTear_Correlation_MoCA-B vs ', colnames(dt)[i+4], '.png')
  ggsave(p, filename = fn, width = 4, height = 4)
  fn <- str_c('Results/Figure4/HumanTear_Correlation_MoCA-B vs ', colnames(dt)[i+4], '.pdf')
  ggsave(p, filename = fn, width = 4, height = 4)
}
# 将ROC阳性的Protein的表达值图示出来
i = 5
for (i in 5:ncol(dt)) {
  myTitle <- gsub("/", "_", colnames(dt)[i])
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(dt, aes(x = Group, fill = Group)) + theme_bw() +
    geom_boxplot(aes_(y = as.name(colnames(dt)[i])), width = 0.6, size = 0.8) +  
    geom_violin(aes_(y = as.name(colnames(dt)[i])), scale = 'width', width = 0.8, alpha = 0.6) +
    geom_jitter(aes_(y = as.name(colnames(dt)[i])), width = 0.3, alpha = 0.8, size = 2) + 
    geom_signif(aes_(y = as.name(colnames(dt)[i])), 
                comparisons = list(c('AD', 'Con')), 
                map_signif_level = F, step_increase = 0.01, 
                tip_length = 0.01, test = "wilcox.test", textsize = 2) +
    geom_hline(yintercept = 4742.664, color = "darkgreen", linetype = 2, size = 0.5) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure4/HumanTear_Raw.Counts_CPM.quantiles_wilcox.test_', myTitle, '.png')
  ggsave(p2, filename = fn2, width = 4, height = 4)
  fn2 <- str_c('Results/Figure4/HumanTear_Raw.Counts_CPM.quantiles_wilcox.test_', myTitle, '.pdf')
  ggsave(p2, filename = fn2, width = 4, height = 4)
}
# 2-4个Pro的联合ROC
data.pos$auc <- ifelse(data.pos$auc > 0.5, data.pos$auc, 1-data.pos$auc)
data.pos <- data.pos %>% arrange(desc(auc))
roc.data <- roc.data %>% select(Group, rownames(data.pos)[1:5])
f2 <- glm(Group ~ roc.data[,2] + roc.data[,3], data = roc.data, family = binomial())
roc.data$Pro2 <- predict(f2, newdata = roc.data, type = "response")
f3 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4], data = roc.data, family = binomial())
roc.data$Pro3 <- predict(f3, newdata = roc.data, type = "response")
f4 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4] + roc.data[,5], data = roc.data, family = binomial())
roc.data$Pro4 <- predict(f4, newdata = roc.data, type = "response")
f5 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4] + roc.data[,5] + roc.data[,6], data = roc.data, family = binomial())
roc.data$Pro5 <- predict(f5, newdata = roc.data, type = "response")
roc.pro2 <- roc(roc.data$Group, roc.data$Pro2)
roc.pro2$auc
roc.pro3 <- roc(roc.data$Group, roc.data$Pro3)
roc.pro3$auc
roc.pro4 <- roc(roc.data$Group, roc.data$Pro4)
roc.pro4$auc
roc.pro5 <- roc(roc.data$Group, roc.data$Pro5)
roc.pro5$auc
png('Results/Figure4/HumanTear_Raw.Counts_CPM.quantiles_ROC_2-5.png', height = 400, width = 400)
plot.roc(roc.pro2, col = "#3C5488CC", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#3C5488CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.35,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro3, col = "#F39B7FCC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F39B7FCC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.3,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro4, col = "#4DBBD5CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#4DBBD5CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.25,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro5, col = "#DC0000CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#DC0000CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.2,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
pdf('Results/Figure4/HumanTear_Raw.Counts_CPM.quantiles_ROC_2-5.pdf', height = 8, width = 8)
plot.roc(roc.pro2, col = "#3C5488CC", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#3C5488CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.35,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro3, col = "#F39B7FCC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F39B7FCC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.3,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro4, col = "#4DBBD5CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#4DBBD5CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.25,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro5, col = "#DC0000CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#DC0000CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.2,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/HumanUrine_Raw.Counts.csv', row.names = 1)
data.nor <- data %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(data.nor) <- rownames(data)
names(data.nor) <- names(data)
write.csv(data.nor, 'Results/HumanUrine_Raw.Counts_cpm.quantiles.csv')
roc.data <- t(data.nor) %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
roc.data <- relocate(roc.data, Group)
auc <- c()
z <- c()
pvalue <- c()
for (j in 1:(ncol(roc.data)-1)){
  rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('AD', 'Con'), direction = '>', data = roc.data)
  auc[j] <- rocPlot$auc
  se <- sqrt(pROC::var(rocPlot))
  b <- auc[j] - 0.5
  z[j] <- (b / se)
  pvalue[j] <- 2 * pt(-abs(z[j]), df=Inf)
  # 对AUC > 0.7 或 < 0.3 的蛋白进行绘图
  if (pvalue[j] < 0.05 & auc[j] >= 0.7) {
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    fn6 <- paste('Results/Figure4/HumanUrine_Raw.Counts_cpm.quantiles_', proteinname, '.png', sep = '')
    png(fn6, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure4/HumanUrine_Raw.Counts_cpm.quantiles_', proteinname, '.pdf', sep = '')
    pdf(fn6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
  }
  if (pvalue[j] < 0.05 & auc[j] <= 0.3) {
    rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('AD', 'Con'), direction = '<', data = roc.data)
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    fn6 <- paste('Results/Figure4/HumanUrine_Raw.Counts_cpm.quantiles_', proteinname, '.png', sep = '')
    png(fn6, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "yellow",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure4/HumanUrine_Raw.Counts_cpm.quantiles_', proteinname, '.pdf', sep = '')
    pdf(fn6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "yellow",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
  }
}
auc_pvalue <- data.frame(auc, z, pvalue)
rownames(auc_pvalue) <- colnames(roc.data)[-1]
auc_pvalue <- arrange(auc_pvalue, desc(auc))
write.csv(auc_pvalue, file = 'Results/Figure4/Files/HumanUrine_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.csv')
auc_pvalue <- filter(auc_pvalue, pvalue < 0.05)
write.csv(auc_pvalue, file = 'Results/Figure4/Files/HumanUrine_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv')
# 将ROC阳性的Protein与MMSE和MoCA-B进行Correlation
rm(list=ls())
gc()
data.nor <- read.csv('Results/HumanUrine_Raw.Counts_cpm.quantiles.csv', row.names = 1)
roc.data <- t(data.nor) %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
roc.data <- relocate(roc.data, Group)
roc.data$ID <- str_remove(rownames(roc.data), 'u')
roc.data$Samples <- rownames(roc.data)
clinical <- read.csv('Human_ClinicalData.csv')
data <- left_join(roc.data, clinical, by = 'ID')
data.pos <- read.csv('Results/Figure4/Files/HumanUrine_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv', row.names = 1)
myPro <- rownames(data.pos)
dt <- select(data, Group, Samples, MMSE, MoCA.B, all_of(myPro))
names(dt) <- names(dt) %>% str_replace_all('/', '.') # 文件名中有/的时候会报错
i = 1
for (i in 1:length(myPro)) {
  p <- ggplot(dt, aes_(x = as.name(colnames(dt)[i+4]))) + theme_bw() +
    geom_point(aes(y = MMSE, color = Group), size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
    stat_smooth(aes(y = MMSE), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
    # stat_poly_eq(aes(y = MMSE, label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x, parse=T,size=5)
    stat_cor(aes(y = MMSE), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = 'darkgreen') +
    geom_vline(xintercept = 1447.946, color = "darkgreen", linetype = 2, size = 0.5) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  fn <- str_c('Results/Figure4/HumanUrine_Correlation_MMSE vs ', colnames(dt)[i+4], '.png')
  ggsave(p, filename = fn, width = 4, height = 4)
  fn <- str_c('Results/Figure4/HumanUrine_Correlation_MMSE vs ', colnames(dt)[i+4], '.pdf')
  ggsave(p, filename = fn, width = 4, height = 4)
  p <- ggplot(dt, aes_(x = as.name(colnames(dt)[i+4]))) + theme_bw() +
    geom_point(aes(y = MoCA.B, color = Group), size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
    stat_smooth(aes(y = MoCA.B), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
    # stat_poly_eq(aes(y = MoCA.B, label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x, parse=T,size=5)
    stat_cor(aes(y = MoCA.B), method = "spearman", digits = 3, size = 4.5, label.y = 0, colour = 'darkblue') +
    ylab('MoCA-B') + 
    geom_vline(xintercept = 1447.946, color = "darkgreen", linetype = 2, size = 0.5) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(colour = 'black', size = 16),
          legend.position = "none")
  fn <- str_c('Results/Figure4/HumanUrine_Correlation_MoCA-B vs ', colnames(dt)[i+4], '.png')
  ggsave(p, filename = fn, width = 4, height = 4)
  fn <- str_c('Results/Figure4/HumanUrine_Correlation_MoCA-B vs ', colnames(dt)[i+4], '.pdf')
  ggsave(p, filename = fn, width = 4, height = 4)
}
# 将ROC阳性的Protein的表达值图示出来
i = 5
for (i in 5:ncol(dt)) {
  myTitle <- gsub("/", "_", colnames(dt)[i])
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(dt, aes(x = Group, fill = Group)) + theme_bw() +
    geom_boxplot(aes_(y = as.name(colnames(dt)[i])), width = 0.6, size = 0.8) +  
    geom_violin(aes_(y = as.name(colnames(dt)[i])), scale = 'width', width = 0.8, alpha = 0.6) +
    geom_jitter(aes_(y = as.name(colnames(dt)[i])), width = 0.3, alpha = 0.8, size = 2) + 
    geom_signif(aes_(y = as.name(colnames(dt)[i])), 
                comparisons = list(c('AD', 'Con')), 
                map_signif_level = F, step_increase = 0.01, 
                tip_length = 0.01, test = "wilcox.test", textsize = 2) +
    geom_hline(yintercept = 1447.946, color = "darkgreen", linetype = 2, size = 0.5) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure4/HumanUrine_Raw.Counts_CPM.quantiles_wilcox.test_', myTitle, '.png')
  ggsave(p2, filename = fn2, width = 4, height = 4)
  fn2 <- str_c('Results/Figure4/HumanUrine_Raw.Counts_CPM.quantiles_wilcox.test_', myTitle, '.pdf')
  ggsave(p2, filename = fn2, width = 4, height = 4)
}
# 2-4个Pro的联合ROC
data.pos$auc <- ifelse(data.pos$auc > 0.5, data.pos$auc, 1-data.pos$auc)
data.pos <- data.pos %>% arrange(desc(auc))
roc.data <- roc.data %>% select(Group, rownames(data.pos)[1:5])
f2 <- glm(Group ~ roc.data[,2] + roc.data[,3], data = roc.data, family = binomial())
roc.data$Pro2 <- predict(f2, newdata = roc.data, type = "response")
f3 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4], data = roc.data, family = binomial())
roc.data$Pro3 <- predict(f3, newdata = roc.data, type = "response")
f4 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4] + roc.data[,5], data = roc.data, family = binomial())
roc.data$Pro4 <- predict(f4, newdata = roc.data, type = "response")
f5 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4] + roc.data[,5] + roc.data[,6], data = roc.data, family = binomial())
roc.data$Pro5 <- predict(f5, newdata = roc.data, type = "response")
roc.pro2 <- roc(roc.data$Group, roc.data$Pro2)
roc.pro2$auc
roc.pro3 <- roc(roc.data$Group, roc.data$Pro3)
roc.pro3$auc
roc.pro4 <- roc(roc.data$Group, roc.data$Pro4)
roc.pro4$auc
roc.pro5 <- roc(roc.data$Group, roc.data$Pro5)
roc.pro5$auc
png('Results/Figure4/HumanUrine_Raw.Counts_CPM.quantiles_ROC_2-5.png', height = 400, width = 400)
plot.roc(roc.pro2, col = "#3C5488CC", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#3C5488CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.35,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro3, col = "#F39B7FCC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F39B7FCC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.3,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro4, col = "#4DBBD5CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#4DBBD5CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.25,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro5, col = "#DC0000CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#DC0000CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.2,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
pdf('Results/Figure4/HumanUrine_Raw.Counts_CPM.quantiles_ROC_2-5.pdf', height = 8, width = 8)
plot.roc(roc.pro2, col = "#3C5488CC", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#3C5488CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.35,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro3, col = "#F39B7FCC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F39B7FCC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.3,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro4, col = "#4DBBD5CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#4DBBD5CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.25,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro5, col = "#DC0000CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#DC0000CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.2,# AUC的坐标为（x，y）
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()

# -- ------------------------------------------------------------------------------------------------------------------
# 7.PCA分析
# -- ------------------------------------------------------------------------------------------------------------------
library(FactoMineR) # PCA分析
library(factoextra) # PCA分析可视化
# -- ------------------------------------------------------------------------------------------------------------------
# 单独做每一种的
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data.Plasma <- read.csv('Results/HumanPlasma_Raw.Counts_cpm.quantiles.csv', row.names = 1)
# 读取阳性的ROC结果
SigPro <- read.csv('Results/Figure4/Files/HumanPlasma_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv', row.names = 1)
myPro <- intersect(rownames(SigPro), rownames(data.Plasma))
data <- data.Plasma %>% filter(rownames(data.Plasma) %in% myPro)
data.pca <- data %>% t() %>% as.data.frame()
data.pca$Group <- ifelse(str_detect(rownames(data.pca),'D'), 'AD', 'Con') %>% factor(levels = c('AD','Con'))
data.pca <- data.pca %>% relocate(Group)
pca <- PCA(data.pca[-1], graph = F)
show_col(c('#AE1F63CC','#D595A7CC','#CC9900CC','#F0E685CC','#6BD76BCC','#749B58CC'))
show_col(pal_lancet(alpha = 0.8)(n = 9))
p <- fviz_pca_ind(pca,
                  mean.point = F, # 去除分组的中心点，否则每个群中间会有一个比较大的点
                  # label = "ind", # 展示每一个样本的标签
                  geom.ind = "point", # 样本采用点来展示
                  repel = T, # 避免重叠
                  pointsize = 3,      # 点的大小
                  pointshape = 21,    # 点的形状
                  fill.ind = data.pca$Group, # 根据样本类型来着色
                  alpha.ind = 0.8,
                  palette = c('#ED0000CC','#00468BCC'), # n个组设置n种颜色
                  title = 'HumanPlasma PCA') + theme_bw() + theme(legend.position = 'none')
ggsave(p, filename = 'Results/Figure4/HumanPlasma_ROC.Pos_PCA.png', width = 4, height = 4.2)
ggsave(p, filename = 'Results/Figure4/HumanPlasma_ROC.Pos_PCA.pdf', width = 4, height = 4.2)

rm(list=ls())
gc()
data.Urine <- read.csv('Results/HumanUrine_Raw.Counts_cpm.quantiles.csv', row.names = 1)
# 读取阳性的ROC结果
SigPro <- read.csv('Results/Figure4/Files/HumanUrine_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv', row.names = 1)
myPro <- intersect(rownames(SigPro), rownames(data.Urine))
data <- data.Urine %>% filter(rownames(data.Urine) %in% myPro)
data.pca <- data %>% t() %>% as.data.frame()
data.pca$Group <- ifelse(str_detect(rownames(data.pca),'D'), 'AD', 'Con') %>% factor(levels = c('AD','Con'))
data.pca <- data.pca %>% relocate(Group)
pca <- PCA(data.pca[-1], graph = F)
show_col(c('#AE1F63CC','#D595A7CC','#CC9900CC','#F0E685CC','#6BD76BCC','#749B58CC'))
show_col(pal_lancet(alpha = 0.8)(n = 9))
p <- fviz_pca_ind(pca,
                  mean.point = F, # 去除分组的中心点，否则每个群中间会有一个比较大的点
                  # label = "ind", # 展示每一个样本的标签
                  geom.ind = "point", # 样本采用点来展示
                  repel = T, # 避免重叠
                  pointsize = 3,      # 点的大小
                  pointshape = 21,    # 点的形状
                  fill.ind = data.pca$Group, # 根据样本类型来着色
                  alpha.ind = 0.8,
                  palette = c('#ED0000CC','#00468BCC'), # n个组设置n种颜色
                  title = 'HumanUrine PCA') + theme_bw() + theme(legend.position = 'none')
ggsave(p, filename = 'Results/Figure4/HumanUrine_ROC.Pos_PCA.png', width = 4, height = 4.2)
ggsave(p, filename = 'Results/Figure4/HumanUrine_ROC.Pos_PCA.pdf', width = 4, height = 4.2)

rm(list=ls())
gc()
data.Saliva <- read.csv('Results/HumanSaliva_Raw.Counts_cpm.quantiles.csv', row.names = 1)
# 读取阳性的ROC结果
SigPro <- read.csv('Results/Figure4/Files/HumanSaliva_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv', row.names = 1)
myPro <- intersect(rownames(SigPro), rownames(data.Saliva))
data <- data.Saliva %>% filter(rownames(data.Saliva) %in% myPro)
data.pca <- data %>% t() %>% as.data.frame()
data.pca$Group <- ifelse(str_detect(rownames(data.pca),'D'), 'AD', 'Con') %>% factor(levels = c('AD','Con'))
data.pca <- data.pca %>% relocate(Group)
pca <- PCA(data.pca[-1], graph = F)
show_col(c('#AE1F63CC','#D595A7CC','#CC9900CC','#F0E685CC','#6BD76BCC','#749B58CC'))
show_col(pal_lancet(alpha = 0.8)(n = 9))
p <- fviz_pca_ind(pca,
                  mean.point = F, # 去除分组的中心点，否则每个群中间会有一个比较大的点
                  # label = "ind", # 展示每一个样本的标签
                  geom.ind = "point", # 样本采用点来展示
                  repel = T, # 避免重叠
                  pointsize = 3,      # 点的大小
                  pointshape = 21,    # 点的形状
                  fill.ind = data.pca$Group, # 根据样本类型来着色
                  alpha.ind = 0.8,
                  palette = c('#ED0000CC','#00468BCC'), # n个组设置n种颜色
                  title = 'HumanSaliva PCA') + theme_bw() + theme(legend.position = 'none')
ggsave(p, filename = 'Results/Figure4/HumanSaliva_ROC.Pos_PCA.png', width = 4, height = 4.2)
ggsave(p, filename = 'Results/Figure4/HumanSaliva_ROC.Pos_PCA.pdf', width = 4, height = 4.2)

rm(list=ls())
gc()
data.Tear <- read.csv('Results/HumanTear_Raw.Counts_cpm.quantiles.csv', row.names = 1)
# 读取阳性的ROC结果
SigPro <- read.csv('Results/Figure4/Files/HumanTear_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv', row.names = 1)
myPro <- intersect(rownames(SigPro), rownames(data.Tear))
data <- data.Tear %>% filter(rownames(data.Tear) %in% myPro)
data.pca <- data %>% t() %>% as.data.frame()
data.pca$Group <- ifelse(str_detect(rownames(data.pca),'D'), 'AD', 'Con') %>% factor(levels = c('AD','Con'))
data.pca <- data.pca %>% relocate(Group)
pca <- PCA(data.pca[-1], graph = F)
show_col(c('#AE1F63CC','#D595A7CC','#CC9900CC','#F0E685CC','#6BD76BCC','#749B58CC'))
show_col(pal_lancet(alpha = 0.8)(n = 9))
p <- fviz_pca_ind(pca,
                  mean.point = F, # 去除分组的中心点，否则每个群中间会有一个比较大的点
                  # label = "ind", # 展示每一个样本的标签
                  geom.ind = "point", # 样本采用点来展示
                  repel = T, # 避免重叠
                  pointsize = 3,      # 点的大小
                  pointshape = 21,    # 点的形状
                  fill.ind = data.pca$Group, # 根据样本类型来着色
                  alpha.ind = 0.8,
                  palette = c('#ED0000CC','#00468BCC'), # n个组设置n种颜色
                  title = 'HumanTear PCA') + theme_bw() + theme(legend.position = 'none')
ggsave(p, filename = 'Results/Figure4/HumanTear_ROC.Pos_PCA.png', width = 4, height = 4.2)
ggsave(p, filename = 'Results/Figure4/HumanTear_ROC.Pos_PCA.pdf', width = 4, height = 4.2)

# -- ------------------------------------------------------------------------------------------------------------------
# 4种体液合在一起做
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/Human_Raw.Counts_cpm.quantiles.csv', row.names = 1)
# 读取阳性的ROC结果
SigPro1 <- read.csv('Results/Figure4/Files/HumanPlasma_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv', row.names = 1)
SigPro2 <- read.csv('Results/Figure4/Files/HumanSaliva_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv', row.names = 1)
SigPro3 <- read.csv('Results/Figure4/Files/HumanUrine_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv', row.names = 1)
SigPro4 <- read.csv('Results/Figure4/Files/HumanTear_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv', row.names = 1)
myPro <- c(rownames(SigPro1),rownames(SigPro2),rownames(SigPro3),rownames(SigPro4)) %>% unique()
myPro <- intersect(myPro, rownames(data))
data <- data %>% filter(rownames(data) %in% myPro)
data.pca <- data %>% t() %>% as.data.frame()
data.pca$Group <- ifelse(str_detect(rownames(data.pca),'D'), 'AD', 'Con') %>% factor(levels = c('AD','Con'))
data.pca$Type <- ifelse(str_detect(rownames(data.pca),'p'), 'Plasma',
                        ifelse(str_detect(rownames(data.pca),'u'), 'Urine',
                               ifelse(str_detect(rownames(data.pca),'s'), 'Saliva','Tear'))) %>% 
  factor(levels = c('Plasma','Saliva', 'Urine', 'Tear'))
data.pca$color <- str_c(data.pca$Type, data.pca$Group, sep = '.') %>% 
  factor(levels = c('Plasma.AD','Plasma.Con','Saliva.AD','Saliva.Con','Urine.AD','Urine.Con','Tear.AD','Tear.Con'))
data.pca <- data.pca %>% relocate(Type, Group, color)
pca <- PCA(data.pca[, -c(1,2,3)], graph = F)
show_col(c('#AE1F63CC','#D595A7CC','#925E9FFF','#925E9FAA','#CC9900CC','#F0E685CC','#6BD76BCC','#749B58CC'))
show_col(pal_lancet(alpha = 0.8)(n = 9))
p <- fviz_pca_ind(pca,
                  mean.point = F, # 去除分组的中心点，否则每个群中间会有一个比较大的点
                  # label = "ind", # 展示每一个样本的标签
                  geom.ind = "point", # 样本采用点来展示
                  repel = T, # 避免重叠
                  pointsize = 3,      # 点的大小
                  pointshape = 21,    # 点的形状
                  fill.ind = data.pca$color, # 根据样本类型来着色
                  alpha.ind = 0.8,
                  palette = c('#AE1F63CC','#D595A7CC','#925E9FFF','#ADB6B6CC','#CC9900CC','#F0E685CC','#6BD76BCC','#749B58CC'), # n个组设置n种颜色
                  title = 'Human PCA') + theme_bw() + theme(legend.title=element_blank())
ggsave(p, filename = 'Results/Figure4/Human_ROC.Pos_PCA.png', width = 5, height = 4)
ggsave(p, filename = 'Results/Figure4/Human_ROC.Pos_PCA.pdf', width = 5, height = 4)

rm(list=ls())
gc()
data <- read.csv('Results/Human_Raw.Counts_cpm.quantiles.csv', row.names = 1)
data.pca <- data %>% t() %>% as.data.frame()
data.pca$Group <- ifelse(str_detect(rownames(data.pca),'D'), 'AD', 'Con') %>% factor(levels = c('AD','Con'))
data.pca$Type <- ifelse(str_detect(rownames(data.pca),'p'), 'Plasma',
                        ifelse(str_detect(rownames(data.pca),'u'), 'Urine',
                               ifelse(str_detect(rownames(data.pca),'s'), 'Saliva','Tear'))) %>% 
  factor(levels = c('Plasma','Saliva', 'Urine', 'Tear'))
data.pca$color <- str_c(data.pca$Type, data.pca$Group, sep = '.') %>% 
  factor(levels = c('Plasma.AD','Plasma.Con','Saliva.AD','Saliva.Con','Urine.AD','Urine.Con','Tear.AD','Tear.Con'))
data.pca <- data.pca %>% relocate(Type, Group, color)
pca <- PCA(data.pca[, -c(1,2,3)], graph = F)
show_col(c('#AE1F63CC','#D595A7CC','#925E9FFF','#925E9FAA','#CC9900CC','#F0E685CC','#6BD76BCC','#749B58CC'))
show_col(pal_lancet(alpha = 0.8)(n = 9))
p <- fviz_pca_ind(pca,
                  mean.point = F, # 去除分组的中心点，否则每个群中间会有一个比较大的点
                  # label = "ind", # 展示每一个样本的标签
                  geom.ind = "point", # 样本采用点来展示
                  repel = T, # 避免重叠
                  pointsize = 3,      # 点的大小
                  pointshape = 21,    # 点的形状
                  fill.ind = data.pca$color, # 根据样本类型来着色
                  alpha.ind = 0.8,
                  palette = c('#AE1F63CC','#D595A7CC','#925E9FFF','#ADB6B6CC','#CC9900CC','#F0E685CC','#6BD76BCC','#749B58CC'), # n个组设置n种颜色
                  title = 'Human PCA') + theme_bw() + theme(legend.title=element_blank())
ggsave(p, filename = 'Results/Figure4/Human_PCA.png', width = 5, height = 4)
ggsave(p, filename = 'Results/Figure4/Human_PCA.pdf', width = 5, height = 4)


# -- ------------------------------------------------------------------------------------------------------------------
# 8.20230218 Tang想看两组表达值的FoldChange
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/HumanPlasma_Raw.Counts_cpm.quantiles.csv', row.names = 1)
conditions <- ifelse(str_detect(colnames(data), 'D'), 'AD', 'Con') %>% factor()
pvalues <- sapply(1:nrow(data),function(i){
  temp <- cbind.data.frame(gene = as.numeric(t(data[i,])), conditions)
  p = wilcox.test(gene~conditions, temp)$p.value
  return(p)
})
fdr <- p.adjust(pvalues, method = "fdr")
conditionsLevel <- levels(conditions)
data.AD <- data[,c(which(conditions==conditionsLevel[1]))]
data.Con <- data[,c(which(conditions==conditionsLevel[2]))]
foldChanges <- rowMeans(data.AD)/rowMeans(data.Con)
Log2FC <- log2(rowMeans(data.AD)/rowMeans(data.Con))
outRst <- data.frame(FoldChange = foldChanges, log2FC = Log2FC, pValues = pvalues, FDR = fdr)
outRst <- na.omit(outRst) %>% arrange(pValues)
write.csv(outRst, 'Results/Figure4/Files/HumanPlasma_ADvsCon_DEpro_wilcox.test.csv')







