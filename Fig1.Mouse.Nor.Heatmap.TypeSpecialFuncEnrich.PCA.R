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
# 1.对Mouse的表达值进行标准化
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse) # 处理数据神包
library(edgeR) # 用于cpm标准化
library(preprocessCore) # normalize.quantiles函数进行分位数标准化
rm(list=ls())
gc()

data <- read.csv('Results/Mouse_Raw.Counts.csv',row.names = 1)
data$mean <- rowMeans(data)
data <- data %>% filter(mean > 0) %>% dplyr::select(!mean) # 将表达量为0的Pro去除
boxplot(data)
data.nor <- data %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(data.nor) <- rownames(data)
names(data.nor) <- names(data)
boxplot(data.nor)
write.csv(data.nor, 'Results/Mouse_Raw.Counts_CPM.quantiles.csv')

# -- ------------------------------------------------------------------------------------------------------------------
# 2.标准化的表达值绘制热图heatmap
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse) # 处理数据神包
library(ggsci) # 期刊配色，单个离散配色模块中最多只有51种颜色
library(scales) # show_col
library(patchwork) # 拼接图
library(ggthemes) # 主题设置
library(pheatmap) # 热图
rm(list=ls())
gc()

data.nor <- read.csv('Results/Mouse_Raw.Counts_CPM.quantiles.csv', row.names = 1)
coldata <- names(data.nor) %>% as.data.frame()
names(coldata) <- 'Samples'
coldata$Group <- ifelse(str_detect(coldata$Samples,'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
coldata$Type <- ifelse(str_detect(coldata$Samples,'s'), 'Serum',
                       ifelse(str_detect(coldata$Samples,'u'), 'Urine', 'NIS')) %>% factor(levels = c('Serum', 'Urine', 'NIS'))
ID <- str_split(coldata$Samples, '_') %>% as.data.frame() %>% t() %>% as.data.frame()
coldata$Month <- ifelse(str_detect(ID$V1,'3'), '3M',
                        ifelse(str_detect(ID$V1,'10'), '10M','15M')) %>% factor(levels = c('3M', '10M', '15M'))
rownames(coldata) <- coldata$Samples
coldata <- coldata[-1]
# show_col(pal_igv(alpha = 0.8)(25))
ann_colors = list(Group = c('AD' ='#CC0000FF', 'Con' ='#6699FFFF'),
                  Type = c('Serum' = '#D595A7CC', 'Urine' = '#F0E685CC', 'NIS' = '#6BD76BCC'), 
                  Month = c('3M' = '#8491B4FF',
                            '10M' = '#4DBBD5FF',
                            '15M' = '#E64B35FF')) #注意ann_colors是列表
pheatmap(data.nor, scale = 'row', color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 6, cellheight = 4,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = T,
         fontsize_number = 4, fontsize_row = 4, fontsize_col = 6, border_color = NA,
         cluster_rows = T, cluster_cols = F, treeheight_row = 40,
         annotation_col = coldata, annotation_colors = ann_colors,
         cutree_rows = 3, gaps_col = c(24,48),
         filename = 'Results/Figure1/Mouse_Raw.Counts_CPM.quantiles_pheatmap.TreeRow.RowName.png') # 保存，自动调整纸张大小
# 获取scale之后的表达矩阵
# 由于scale函数默认对列进行归一化，因此这里做了两次转置；
data.nor_scaled <- t(scale(t(data.nor))) %>% as.data.frame()
#查看归一化后的数据前6行；
head(data.nor_scaled)
write.csv(data.nor_scaled, 'Results/Figure1/Files/Mouse_Raw.Counts_CPM.quantiles_pheatmap.scaleRow.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# Tang说按照组内scale之后表达值大于1的多于n个，组外表达值大于1的小于n分为标准去定义Type Special
# -- ------------------------------------------------------------------------------------------------------------------
# library(tidyverse) # 处理数据神包
# library(edgeR) # 用于cpm标准化
# library(preprocessCore) # normalize.quantiles函数进行分位数标准化
# rm(list=ls())
# gc()
# data.nor <- read.csv('Results/Mouse_Raw.Counts_CPM.quantiles_scaleRow.csv', row.names = 1)
# range(data.nor)
# data.Serum <- data.nor[25:48]
# data.Con <- data.nor[,c(1:24,49:72)]
# data.Serum.T.F <- data.Serum >= 0.5
# data.Serum$Count <- rowSums(data.Serum.T.F)
# data.Serum <- data.Serum %>% filter(Count >= 16)
# data.Con.T.F <- data.Con >= 0.5
# data.Con$Count <- rowSums(data.Con.T.F)
# data.Con <- data.Con %>% filter(Count <= 16)
# data.Serum$Pro <- rownames(data.Serum)
# data.Con$Pro <- rownames(data.Con)
# data.Serum.T <- inner_join(data.Serum, data.Con, by = 'Pro')
# write_csv(data.Serum.T, 'Results/Mouse.PT_Samples_Raw.Counts_CPM.quantiles_scaleRow_SerumSpecial.Tang.csv')
# 
# data.Urine <- data.nor[49:72]
# data.Con <- data.nor[1:48]
# data.Urine.T.F <- data.Urine >= 0.5
# data.Urine$Count <- rowSums(data.Urine.T.F)
# data.Urine <- data.Urine %>% filter(Count >= 16)
# data.Con.T.F <- data.Con >= 0.5
# data.Con$Count <- rowSums(data.Con.T.F)
# data.Con <- data.Con %>% filter(Count <= 16)
# data.Urine$Pro <- rownames(data.Urine)
# data.Con$Pro <- rownames(data.Con)
# data.Urine.T <- inner_join(data.Urine, data.Con, by = 'Pro')
# write_csv(data.Urine.T, 'Results/Mouse.PT_Samples_Raw.Counts_CPM.quantiles_scaleRow_UrineSpecial.Tang.csv')
# 
# data.NIS <- data.nor[1:24]
# data.Con <- data.nor[25:72]
# data.NIS.T.F <- data.NIS >= 0.5
# data.NIS$Count <- rowSums(data.NIS.T.F)
# data.NIS <- data.NIS %>% filter(Count > 16)
# data.Con.T.F <- data.Con >= 0.5
# data.Con$Count <- rowSums(data.Con.T.F)
# data.Con <- data.Con %>% filter(Count < 16)
# data.NIS$Pro <- rownames(data.NIS)
# data.Con$Pro <- rownames(data.Con)
# data.NIS.T <- inner_join(data.NIS, data.Con, by = 'Pro')
# write_csv(data.NIS.T, 'Results/Mouse.PT_Samples_Raw.Counts_CPM.quantiles_scaleRow_NISSpecial.Tang.csv')
# 
# myPro <- c(data.Serum.T$Pro, data.Urine.T$Pro, data.NIS.T$Pro)
# myPro %>% unique()
# data.nor <- data.nor %>% filter(rownames(data.nor) %in% myPro)
# coldata <- names(data.nor) %>% as.data.frame()
# names(coldata) <- 'Samples'
# coldata$Group <- ifelse(str_detect(coldata$Samples,'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
# coldata$Type <- ifelse(str_detect(coldata$Samples,'s'), 'Serum',
#                        ifelse(str_detect(coldata$Samples,'u'), 'Urine', 'NIS')) %>% factor(levels = c('Serum', 'Urine', 'NIS'))
# rownames(coldata) <- coldata$Samples
# coldata <- coldata[-1]
# # show_col(pal_igv(alpha = 0.8)(25))
# ann_colors = list(Group = c('AD' ='#CC0000FF', 'Con' ='#6699FFFF'),
#                   Type = c('Serum' = '#D595A7CC', 'Urine' = '#F0E685CC', 'NIS' = '#6BD76BCC')) #注意ann_colors是列表
# pheatmap(data.nor, color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 6,
#          display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = T,
#          fontsize_number = 4, fontsize_row = 6, fontsize_col = 8, border_color = NA,
#          cluster_rows = T, cluster_cols = F, treeheight_row = 40,
#          annotation_col = coldata, annotation_colors = ann_colors,
#          cutree_rows = 3,
#          gaps_col = c(24,48),
#          # gaps_row = c(9,11),
#          filename = 'Results/Mouse.PT_Raw.Counts_CPM.quantiles.TreeRow.RowName_TypeSpecial_Tang.png') # 保存，自动调整纸张大小
# pheatmap(data.nor, color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 6,
#          display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = T,
#          fontsize_number = 4, fontsize_row = 6, fontsize_col = 8, border_color = NA,
#          cluster_rows = T, cluster_cols = F, treeheight_row = 40,
#          annotation_col = coldata, annotation_colors = ann_colors,
#          cutree_rows = 3,
#          gaps_col = c(24,48),
#          # gaps_row = c(9,11),
#          filename = 'Results/Mouse.PT_Raw.Counts_CPM.quantiles.TreeRow.RowName_TypeSpecial_Tang.pdf') # 保存，自动调整纸张大小
# -- ------------------------------------------------------------------------------------------------------------------
# 整体的画在一起太丑了，分开每个时间点画
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/Mouse_Raw.Counts.csv', row.names = 1)
coldata <- names(data) %>% as.data.frame()
names(coldata) <- 'Samples'
coldata$Group <- ifelse(str_detect(coldata$Samples,'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
coldata$Type <- ifelse(str_detect(coldata$Samples,'s'), 'Serum',
                       ifelse(str_detect(coldata$Samples,'u'), 'Urine', 'NIS')) %>% factor(levels = c('Serum', 'Urine', 'NIS'))
ID <- str_split(coldata$Samples, '_') %>% as.data.frame() %>% t() %>% as.data.frame()
coldata$Month <- ifelse(str_detect(ID$V1,'3'), '3M',
                        ifelse(str_detect(ID$V1,'10'), '10M','15M')) %>% factor(levels = c('3M', '10M', '15M'))
rownames(coldata) <- coldata$Samples
coldata <- coldata[-1]
# 15M
coldata.sub <- coldata %>% filter(Month == '15M') %>% dplyr::select(!Month)
dt <- data %>% dplyr::select(rownames(coldata.sub))
dt$mean <- rowMeans(dt)
range(dt$mean)
dt <- dt %>% filter(mean > 0) %>% dplyr::select(!mean) # 将表达量为0的Pro去除
dt.nor <- dt %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(dt.nor) <- rownames(dt)
names(dt.nor) <- names(dt)
write.csv(dt.nor, 'Results/Mouse.15M_Raw.Counts_CPM.quantiles.csv')
dt.nor$var <- apply(dt.nor, 1, var) # 对每一行求方差
range(dt.nor$var)
dt.nor <- dt.nor %>% filter(var > 0) %>% dplyr::select(!var) # 将表达量全相等的Pro去除
ann_colors.sub = list(Group = c('AD' ='#CC0000FF', 'Con' ='#6699FFFF'),
                      Type = c('Serum' = '#D595A7CC', 'Urine' = '#F0E685CC', 'NIS' = '#6BD76BCC')) #注意ann_colors是列表
pheatmap(dt.nor, scale = 'row', color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 2.5,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = F,
         fontsize_number = 4, fontsize_row = 2.5, fontsize_col = 10,
         cluster_rows = T, cluster_cols = T, treeheight_row = 30, treeheight_col = 20,
         annotation_col = coldata.sub, annotation_colors = ann_colors.sub,
         # cutree_rows = 3, gaps_col = c(8,16),
         filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_pheatmap.TreeRow.RowName 2.png') # 保存，自动调整纸张大小
pheatmap(dt.nor, scale = 'row', color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 2.5,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = F,
         fontsize_number = 4, fontsize_row = 2.5, fontsize_col = 10,
         cluster_rows = T, cluster_cols = T, treeheight_row = 30, treeheight_col = 20,
         annotation_col = coldata.sub, annotation_colors = ann_colors.sub,
         # cutree_rows = 3, gaps_col = c(8,16),
         filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_pheatmap.TreeRow.RowName 2.pdf') # 保存，自动调整纸张大小
pheatmap(dt.nor, scale = 'row', color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 2.5,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = T,
         fontsize_number = 4, fontsize_row = 2.5, fontsize_col = 10,
         cluster_rows = T, cluster_cols = F, treeheight_row = 30,
         annotation_col = coldata.sub, annotation_colors = ann_colors.sub,
         cutree_rows = 3, gaps_col = c(8,16),
         filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_pheatmap.TreeRow.RowName.png') # 保存，自动调整纸张大小
pheatmap(dt.nor, scale = 'row', color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 2.5,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = T,
         fontsize_number = 4, fontsize_row = 2.5, fontsize_col = 10,
         cluster_rows = T, cluster_cols = F, treeheight_row = 30,
         annotation_col = coldata.sub, annotation_colors = ann_colors.sub,
         cutree_rows = 3, gaps_col = c(8,16),
         filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_pheatmap.TreeRow.RowName.pdf') # 保存，自动调整纸张大小
# 获取scale之后的表达矩阵
# 由于scale函数默认对列进行归一化，因此这里做了两次转置；
dt.nor_scaled <- t(scale(t(dt.nor))) %>% as.data.frame()
#查看归一化后的数据前6行；
head(dt.nor_scaled)
write.csv(dt.nor_scaled, 'Results/Figure1/Files/Mouse.15M_Raw.Counts_CPM.quantiles_scaleRow.csv')
# 10M
coldata.sub <- coldata %>% filter(Month == '10M') %>% dplyr::select(!Month)
dt <- data %>% dplyr::select(rownames(coldata.sub))
dt$mean <- rowMeans(dt)
range(dt$mean)
dt <- dt %>% filter(mean > 0) %>% dplyr::select(!mean) # 将表达量为0的Pro去除
dt.nor <- dt %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(dt.nor) <- rownames(dt)
names(dt.nor) <- names(dt)
write.csv(dt.nor, 'Results/Mouse.10M_Raw.Counts_CPM.quantiles.csv')
dt.nor$var <- apply(dt.nor, 1, var) # 对每一行求方差
range(dt.nor$var)
dt.nor <- dt.nor %>% filter(var > 0) %>% dplyr::select(!var) # 将表达量全相等的Pro去除
ann_colors.sub = list(Group = c('AD' ='#CC0000FF', 'Con' ='#6699FFFF'),
                      Type = c('Serum' = '#D595A7CC', 'Urine' = '#F0E685CC', 'NIS' = '#6BD76BCC')) #注意ann_colors是列表
pheatmap(dt.nor, scale = 'row', color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 2.5,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = T,
         fontsize_number = 4, fontsize_row = 2.5, fontsize_col = 10,
         cluster_rows = T, cluster_cols = F, treeheight_row = 30,
         annotation_col = coldata.sub, annotation_colors = ann_colors.sub,
         cutree_rows = 3, gaps_col = c(8,16),
         filename = 'Results/Figure1/Mouse.10M_Raw.Counts_CPM.quantiles_pheatmap.TreeRow.RowName.png') # 保存，自动调整纸张大小
pheatmap(dt.nor, scale = 'row', color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 2.5,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = T,
         fontsize_number = 4, fontsize_row = 2.5, fontsize_col = 10,
         cluster_rows = T, cluster_cols = F, treeheight_row = 30,
         annotation_col = coldata.sub, annotation_colors = ann_colors.sub,
         cutree_rows = 3, gaps_col = c(8,16),
         filename = 'Results/Figure1/Mouse.10M_Raw.Counts_CPM.quantiles_pheatmap.TreeRow.RowName.pdf') # 保存，自动调整纸张大小
# 获取scale之后的表达矩阵
# 由于scale函数默认对列进行归一化，因此这里做了两次转置；
dt.nor_scaled <- t(scale(t(dt.nor))) %>% as.data.frame()
#查看归一化后的数据前6行；
head(dt.nor_scaled)
write.csv(dt.nor_scaled, 'Results/Figure1/Files/Mouse.10M_Raw.Counts_CPM.quantiles_scaleRow.csv')

# -- ------------------------------------------------------------------------------------------------------------------
# Tang说按照组内scale之后表达值大于1的多于5个，组外表达值大于1的小于6分为标准去定义Type Special
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse) # 处理数据神包
library(edgeR) # 用于cpm标准化
library(preprocessCore) # normalize.quantiles函数进行分位数标准化
# -- ------------------------------------------------------------------------------------------------------------------
# 15M
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data.nor <- read.csv('Results/Figure1/Files/Mouse.15M_Raw.Counts_CPM.quantiles_scaleRow.csv', row.names = 1)
range(data.nor)
data.Serum <- data.nor[9:16]
data.Con <- data.nor[,c(1:8,17:24)]
data.Serum.T.F <- data.Serum >= 1
data.Serum$Count <- rowSums(data.Serum.T.F)
data.Serum <- data.Serum %>% filter(Count >= 5)
data.Con.T.F <- data.Con >= 1
data.Con$Count <- rowSums(data.Con.T.F)
data.Con <- data.Con %>% filter(Count <= 6)
data.Serum$Pro <- rownames(data.Serum)
data.Con$Pro <- rownames(data.Con)
data.Serum.T <- inner_join(data.Serum, data.Con, by = 'Pro')
write_csv(data.Serum.T, 'Results/Figure1/Files/Mouse.15M_Raw.Counts_CPM.quantiles_scaleRow_SerumSpecial.csv')

data.Urine <- data.nor[17:24]
data.Con <- data.nor[1:16]
data.Urine.T.F <- data.Urine >= 1
data.Urine$Count <- rowSums(data.Urine.T.F)
data.Urine <- data.Urine %>% filter(Count >= 5)
data.Con.T.F <- data.Con >= 1
data.Con$Count <- rowSums(data.Con.T.F)
data.Con <- data.Con %>% filter(Count <= 6)
data.Urine$Pro <- rownames(data.Urine)
data.Con$Pro <- rownames(data.Con)
data.Urine.T <- inner_join(data.Urine, data.Con, by = 'Pro')
write_csv(data.Urine.T, 'Results/Figure1/Files/Mouse.15M_Raw.Counts_CPM.quantiles_scaleRow_UrineSpecial.csv')

data.NIS <- data.nor[1:8]
data.Con <- data.nor[9:24]
data.NIS.T.F <- data.NIS >= 1
data.NIS$Count <- rowSums(data.NIS.T.F)
data.NIS <- data.NIS %>% filter(Count >= 5)
data.Con.T.F <- data.Con >= 1
data.Con$Count <- rowSums(data.Con.T.F)
data.Con <- data.Con %>% filter(Count <= 6)
data.NIS$Pro <- rownames(data.NIS)
data.Con$Pro <- rownames(data.Con)
data.NIS.T <- inner_join(data.NIS, data.Con, by = 'Pro')
write_csv(data.NIS.T, 'Results/Figure1/Files/Mouse.15M_Raw.Counts_CPM.quantiles_scaleRow_NISSpecial.csv')

myPro <- c(data.Serum.T$Pro, data.Urine.T$Pro, data.NIS.T$Pro)
myPro %>% unique()
data.nor <- data.nor %>% filter(rownames(data.nor) %in% myPro)
coldata <- names(data.nor) %>% as.data.frame()
names(coldata) <- 'Samples'
coldata$Group <- ifelse(str_detect(coldata$Samples,'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
coldata$Type <- ifelse(str_detect(coldata$Samples,'s'), 'Serum',
                       ifelse(str_detect(coldata$Samples,'u'), 'Urine', 'NIS')) %>% factor(levels = c('Serum', 'Urine', 'NIS'))
rownames(coldata) <- coldata$Samples
coldata <- coldata[-1]
# show_col(pal_igv(alpha = 0.8)(25))
ann_colors = list(Group = c('AD' ='#CC0000FF', 'Con' ='#6699FFFF'),
                  Type = c('Serum' = '#D595A7CC', 'Urine' = '#F0E685CC', 'NIS' = '#6BD76BCC')) #注意ann_colors是列表
pheatmap(data.nor, color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 6,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = T,
         fontsize_number = 4, fontsize_row = 6, fontsize_col = 8, border_color = NA,
         cluster_rows = T, cluster_cols = F, treeheight_row = 40,
         annotation_col = coldata, annotation_colors = ann_colors,
         cutree_rows = 3,
         gaps_col = c(8,16),
         filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_pheatmap.TreeRow.RowName_TypeSpecial.png') # 保存，自动调整纸张大小
pheatmap(data.nor, color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 6,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = T,
         fontsize_number = 4, fontsize_row = 6, fontsize_col = 8, border_color = NA,
         cluster_rows = T, cluster_cols = F, treeheight_row = 40,
         annotation_col = coldata, annotation_colors = ann_colors,
         cutree_rows = 3,
         gaps_col = c(8,16),
         filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_pheatmap.TreeRow.RowName_TypeSpecial.pdf') # 保存，自动调整纸张大小
# -- ------------------------------------------------------------------------------------------------------------------
# 组织特异性的蛋白进行功能分析
# -- ------------------------------------------------------------------------------------------------------------------
library(clusterProfiler)
library(enrichplot) #用于可视化的包
library(org.Hs.eg.db) #小鼠用org.Mm.eg.db
# -- ------------------------------------------------------------------------------------------------------------------
# Serum
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/Figure1/Files/Mouse.15M_Raw.Counts_CPM.quantiles_scaleRow_SerumSpecial.csv')
gene <- data$Pro
gene <-  str_replace(gene, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2') %>%
  str_replace('PAR2', 'F2RL1')
genelist <- bitr(gene, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
genelist$SYMBOL %>% unique()
myGene <- genelist$ENTREZID %>% unique()
# -- ------------------------------------------------------------------------
# GO富集分析：
# -- ------------------------------------------------------------------------
# GO over-representation analysis
# -- ------------------------------------------------------------------------
# Biological Process
egoBP_Serum <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoBP_Serum <- clusterProfiler::filter(egoBP_Serum, p.adjust < 0.05)
dim(egoBP_Serum@result)
head(egoBP_Serum@result$Description)
write_csv(egoBP_Serum@result, 'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-BP.Sig.csv')
egoBP_Serum.sim0.8 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.8)
egoBP_Serum.sim0.7 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.7)
egoBP_Serum.sim0.6 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.6)
egoBP_Serum.sim0.5 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.5)
egoBP_Serum.sim0.4 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.4)
egoBP_Serum.sim0.3 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.3)
egoBP_Serum.sim0.8@result$Description %>% head(n=10)
egoBP_Serum.sim0.7@result$Description %>% head(n=10)
egoBP_Serum.sim0.6@result$Description %>% head(n=10)
egoBP_Serum.sim0.5@result$Description %>% head(n=10)
egoBP_Serum.sim0.4@result$Description %>% head(n=10)
egoBP_Serum.sim0.3@result$Description %>% head(n=10)
write.csv(egoBP_Serum.sim0.6@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-BP.Sig.sim0.6.csv', row.names = F)
write.csv(egoBP_Serum.sim0.5@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-BP.Sig.sim0.5.csv', row.names = F)
# Cellular Component
egoCC_Serum <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoCC_Serum <- clusterProfiler::filter(egoCC_Serum, p.adjust < 0.05)
dim(egoCC_Serum@result)
head(egoCC_Serum@result$Description)
write_csv(egoCC_Serum@result, 'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-CC.Sig.csv')
egoCC_Serum.sim0.8 <- clusterProfiler::simplify(egoCC_Serum, cutoff = 0.8)
egoCC_Serum.sim0.7 <- clusterProfiler::simplify(egoCC_Serum, cutoff = 0.7)
egoCC_Serum.sim0.6 <- clusterProfiler::simplify(egoCC_Serum, cutoff = 0.6)
egoCC_Serum.sim0.5 <- clusterProfiler::simplify(egoCC_Serum, cutoff = 0.5)
egoCC_Serum.sim0.8@result$Description %>% head(n=10)
egoCC_Serum.sim0.7@result$Description %>% head(n=10)
egoCC_Serum.sim0.6@result$Description %>% head(n=10)
egoCC_Serum.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC_Serum.sim0.6@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-CC.Sig.sim0.6.csv', row.names = F)
write.csv(egoCC_Serum.sim0.5@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-CC.Sig.sim0.5.csv', row.names = F)
# Molecular Function
egoMF_Serum <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoMF_Serum <- clusterProfiler::filter(egoMF_Serum, p.adjust < 0.05)
dim(egoMF_Serum@result)
head(egoMF_Serum@result$Description)
write_csv(egoMF_Serum@result, 'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-MF.Sig.csv')
egoMF_Serum_5 <- gofilter(egoMF_Serum, level = 5)
egoMF_Serum_5@result$Description %>% head(n=10)
write.csv(egoMF_Serum_5@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-MF.Sig.Level5.csv', row.names = F)
egoMF_Serum.sim0.8 <- clusterProfiler::simplify(egoMF_Serum, cutoff = 0.8)
egoMF_Serum.sim0.7 <- clusterProfiler::simplify(egoMF_Serum, cutoff = 0.7)
egoMF_Serum.sim0.6 <- clusterProfiler::simplify(egoMF_Serum, cutoff = 0.6)
egoMF_Serum.sim0.5 <- clusterProfiler::simplify(egoMF_Serum, cutoff = 0.5)
egoMF_Serum.sim0.8@result$Description %>% head(n=10)
egoMF_Serum.sim0.7@result$Description %>% head(n=10)
egoMF_Serum.sim0.6@result$Description %>% head(n=10)
egoMF_Serum.sim0.5@result$Description %>% head(n=10)
write.csv(egoMF_Serum.sim0.6@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-MF.Sig.sim0.6.csv', row.names = F)
write.csv(egoMF_Serum.sim0.5@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-MF.Sig.sim0.5.csv', row.names = F)
# 保存数据，后续作图用
save(data, myGene, 
     egoBP_Serum, egoCC_Serum, egoMF_Serum, 
     egoBP_Serum.sim0.6, egoCC_Serum.sim0.6, egoMF_Serum.sim0.6,
     egoBP_Serum.sim0.5, egoCC_Serum.sim0.5, egoMF_Serum.sim0.5,
     file = "Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO.RData")
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/Figure1/Files/Mouse.15M_Raw.Counts_CPM.quantiles_scaleRow_UrineSpecial.csv')
gene <- data$Pro
gene <-  str_replace(gene, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2') %>%
  str_replace('PAR2', 'F2RL1')
genelist <- bitr(gene, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
genelist$SYMBOL %>% unique()
myGene <- genelist$ENTREZID %>% unique()
# -- ------------------------------------------------------------------------
# GO富集分析：
# -- ------------------------------------------------------------------------
# GO over-representation analysis
# -- ------------------------------------------------------------------------
# Biological Process
egoBP_Urine <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoBP_Urine <- clusterProfiler::filter(egoBP_Urine, p.adjust < 0.05)
dim(egoBP_Urine@result)
head(egoBP_Urine@result$Description)
write_csv(egoBP_Urine@result, 'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-BP.Sig.csv')
egoBP_Urine.sim0.8 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.8)
egoBP_Urine.sim0.7 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.7)
egoBP_Urine.sim0.6 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.6)
egoBP_Urine.sim0.5 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.5)
egoBP_Urine.sim0.4 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.4)
egoBP_Urine.sim0.3 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.3)
egoBP_Urine.sim0.8@result$Description %>% head(n=10)
egoBP_Urine.sim0.7@result$Description %>% head(n=10)
egoBP_Urine.sim0.6@result$Description %>% head(n=10)
egoBP_Urine.sim0.5@result$Description %>% head(n=10)
egoBP_Urine.sim0.4@result$Description %>% head(n=10)
egoBP_Urine.sim0.3@result$Description %>% head(n=10)
write.csv(egoBP_Urine.sim0.6@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-BP.Sig.sim0.6.csv', row.names = F)
write.csv(egoBP_Urine.sim0.5@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-BP.Sig.sim0.5.csv', row.names = F)
# Cellular Component
egoCC_Urine <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoCC_Urine <- clusterProfiler::filter(egoCC_Urine, p.adjust < 0.05)
dim(egoCC_Urine@result)
head(egoCC_Urine@result$Description)
write_csv(egoCC_Urine@result, 'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-CC.Sig.csv')
egoCC_Urine.sim0.8 <- clusterProfiler::simplify(egoCC_Urine, cutoff = 0.8)
egoCC_Urine.sim0.7 <- clusterProfiler::simplify(egoCC_Urine, cutoff = 0.7)
egoCC_Urine.sim0.6 <- clusterProfiler::simplify(egoCC_Urine, cutoff = 0.6)
egoCC_Urine.sim0.5 <- clusterProfiler::simplify(egoCC_Urine, cutoff = 0.5)
egoCC_Urine.sim0.8@result$Description %>% head(n=10)
egoCC_Urine.sim0.7@result$Description %>% head(n=10)
egoCC_Urine.sim0.6@result$Description %>% head(n=10)
egoCC_Urine.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC_Urine.sim0.6@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-CC.Sig.sim0.6.csv', row.names = F)
write.csv(egoCC_Urine.sim0.5@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-CC.Sig.sim0.5.csv', row.names = F)
# Molecular Function
egoMF_Urine <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoMF_Urine <- clusterProfiler::filter(egoMF_Urine, p.adjust < 0.05)
dim(egoMF_Urine@result)
head(egoMF_Urine@result$Description)
write_csv(egoMF_Urine@result, 'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-MF.Sig.csv')
egoMF_Urine_5 <- gofilter(egoMF_Urine, level = 5)
egoMF_Urine_5@result$Description %>% head(n=10)
write.csv(egoMF_Urine_5@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-MF.Sig.Level5.csv', row.names = F)
egoMF_Urine.sim0.8 <- clusterProfiler::simplify(egoMF_Urine, cutoff = 0.8)
egoMF_Urine.sim0.7 <- clusterProfiler::simplify(egoMF_Urine, cutoff = 0.7)
egoMF_Urine.sim0.6 <- clusterProfiler::simplify(egoMF_Urine, cutoff = 0.6)
egoMF_Urine.sim0.5 <- clusterProfiler::simplify(egoMF_Urine, cutoff = 0.5)
egoMF_Urine.sim0.8@result$Description %>% head(n=10)
egoMF_Urine.sim0.7@result$Description %>% head(n=10)
egoMF_Urine.sim0.6@result$Description %>% head(n=10)
egoMF_Urine.sim0.5@result$Description %>% head(n=10)
write.csv(egoMF_Urine.sim0.6@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-MF.Sig.sim0.6.csv', row.names = F)
write.csv(egoMF_Urine.sim0.5@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-MF.Sig.sim0.5.csv', row.names = F)
# 保存数据，后续作图用
save(data, myGene, 
     egoBP_Urine, egoCC_Urine, egoMF_Urine, 
     egoBP_Urine.sim0.6, egoCC_Urine.sim0.6, egoMF_Urine.sim0.6,
     egoBP_Urine.sim0.5, egoCC_Urine.sim0.5, egoMF_Urine.sim0.5,
     file = "Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO.RData")
# -- ------------------------------------------------------------------------------------------------------------------
# NIS
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/Figure1/Files/Mouse.15M_Raw.Counts_CPM.quantiles_scaleRow_NISSpecial.csv')
gene <- data$Pro
gene <-  str_replace(gene, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2') %>%
  str_replace('PAR2', 'F2RL1')
genelist <- bitr(gene, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
genelist$SYMBOL %>% unique()
myGene <- genelist$ENTREZID %>% unique()
# -- ------------------------------------------------------------------------
# GO富集分析：
# -- ------------------------------------------------------------------------
# GO over-representation analysis
# -- ------------------------------------------------------------------------
# Biological Process
egoBP_NIS <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoBP_NIS <- clusterProfiler::filter(egoBP_NIS, p.adjust < 0.05)
dim(egoBP_NIS@result)
head(egoBP_NIS@result$Description)
write_csv(egoBP_NIS@result, 'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-BP.Sig.csv')
egoBP_NIS.sim0.8 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.8)
egoBP_NIS.sim0.7 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.7)
egoBP_NIS.sim0.6 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.6)
egoBP_NIS.sim0.5 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.5)
egoBP_NIS.sim0.4 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.4)
egoBP_NIS.sim0.3 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.3)
egoBP_NIS.sim0.8@result$Description %>% head(n=10)
egoBP_NIS.sim0.7@result$Description %>% head(n=10)
egoBP_NIS.sim0.6@result$Description %>% head(n=10)
egoBP_NIS.sim0.5@result$Description %>% head(n=10)
egoBP_NIS.sim0.4@result$Description %>% head(n=10)
egoBP_NIS.sim0.3@result$Description %>% head(n=10)
write.csv(egoBP_NIS.sim0.6@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-BP.Sig.sim0.6.csv', row.names = F)
write.csv(egoBP_NIS.sim0.5@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-BP.Sig.sim0.5.csv', row.names = F)
# Cellular Component
egoCC_NIS <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoCC_NIS <- clusterProfiler::filter(egoCC_NIS, p.adjust < 0.05)
dim(egoCC_NIS@result)
head(egoCC_NIS@result$Description)
write_csv(egoCC_NIS@result, 'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-CC.Sig.csv')
egoCC_NIS.sim0.8 <- clusterProfiler::simplify(egoCC_NIS, cutoff = 0.8)
egoCC_NIS.sim0.7 <- clusterProfiler::simplify(egoCC_NIS, cutoff = 0.7)
egoCC_NIS.sim0.6 <- clusterProfiler::simplify(egoCC_NIS, cutoff = 0.6)
egoCC_NIS.sim0.5 <- clusterProfiler::simplify(egoCC_NIS, cutoff = 0.5)
egoCC_NIS.sim0.8@result$Description %>% head(n=10)
egoCC_NIS.sim0.7@result$Description %>% head(n=10)
egoCC_NIS.sim0.6@result$Description %>% head(n=10)
egoCC_NIS.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC_NIS.sim0.6@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-CC.Sig.sim0.6.csv', row.names = F)
write.csv(egoCC_NIS.sim0.5@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-CC.Sig.sim0.5.csv', row.names = F)
# Molecular Function
egoMF_NIS <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoMF_NIS <- clusterProfiler::filter(egoMF_NIS, p.adjust < 0.05)
dim(egoMF_NIS@result)
head(egoMF_NIS@result$Description)
write_csv(egoMF_NIS@result, 'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-MF.Sig.csv')
egoMF_NIS_5 <- gofilter(egoMF_NIS, level = 5)
egoMF_NIS_5@result$Description %>% head(n=10)
write.csv(egoMF_NIS_5@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-MF.Sig.Level5.csv', row.names = F)
egoMF_NIS.sim0.8 <- clusterProfiler::simplify(egoMF_NIS, cutoff = 0.8)
egoMF_NIS.sim0.7 <- clusterProfiler::simplify(egoMF_NIS, cutoff = 0.7)
egoMF_NIS.sim0.6 <- clusterProfiler::simplify(egoMF_NIS, cutoff = 0.6)
egoMF_NIS.sim0.5 <- clusterProfiler::simplify(egoMF_NIS, cutoff = 0.5)
egoMF_NIS.sim0.8@result$Description %>% head(n=10)
egoMF_NIS.sim0.7@result$Description %>% head(n=10)
egoMF_NIS.sim0.6@result$Description %>% head(n=10)
egoMF_NIS.sim0.5@result$Description %>% head(n=10)
write.csv(egoMF_NIS.sim0.6@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-MF.Sig.sim0.6.csv', row.names = F)
write.csv(egoMF_NIS.sim0.5@result, 
          'Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-MF.Sig.sim0.5.csv', row.names = F)
# 保存数据，后续作图用
save(data, myGene, 
     egoBP_NIS, egoCC_NIS, egoMF_NIS, 
     egoBP_NIS.sim0.6, egoCC_NIS.sim0.6, egoMF_NIS.sim0.6,
     egoBP_NIS.sim0.5, egoCC_NIS.sim0.5, egoMF_NIS.sim0.5,
     file = "Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO.RData")

# -- ------------------------------------------------------------------------
# GO over-representation analysis可视化
# -- ------------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db) #小鼠用org.Mm.eg.db
library(enrichplot) #用于可视化的包
library(ggsci) # 期刊配色，单个离散配色模块中最多只有51种颜色
library(scales) # show_col
library(patchwork) # 拼接图
library(ggthemes) # 主题设置
# -- ------------------------------------------------------------------------------------------------------------------
# Serum
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO.RData')
# ggplot2画图
show_col(pal_npg(alpha = 0.8)(10))
my_res.bp <- egoBP_Serum.sim0.6@result
p <- my_res.bp %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.bp$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#F39B7FFF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) +
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Serum GO-BP') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)
my_res.cc <- egoCC_Serum.sim0.6@result
p <- my_res.cc %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.cc$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#4DBBD5FF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) + 
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Serum GO-CC') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)
my_res.mf <- egoMF_Serum.sim0.6@result
p <- my_res.mf %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.mf$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#8491B4FF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) + 
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Serum GO-MF') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO.RData')
# ggplot2画图
show_col(pal_npg(alpha = 0.8)(10))
my_res.bp <- egoBP_Urine.sim0.6@result
p <- my_res.bp %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.bp$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#F39B7FFF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) +
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Urine GO-BP') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.png',
       width = 6, height = 2)
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.pdf',
       width = 6, height = 2)
my_res.cc <- egoCC_Urine.sim0.6@result
p <- my_res.cc %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.cc$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#4DBBD5FF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) + 
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Urine GO-CC') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)
my_res.mf <- egoMF_Urine.sim0.6@result
p <- my_res.mf %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.mf$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#8491B4FF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) + 
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Urine GO-MF') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.png',
       width = 8, height = 2)
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.pdf',
       width = 8, height = 2)
# -- ------------------------------------------------------------------------------------------------------------------
# NIS
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO.RData')
# ggplot2画图
show_col(pal_npg(alpha = 0.8)(10))
my_res.bp <- egoBP_NIS.sim0.6@result
p <- my_res.bp %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.bp$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#F39B7FFF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) +
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('NIS GO-BP') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.png',
       width = 5, height = 2)
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.pdf',
       width = 5, height = 2)
my_res.cc <- egoCC_NIS.sim0.6@result
p <- my_res.cc %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.cc$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#4DBBD5FF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) + 
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('NIS GO-CC') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)
my_res.mf <- egoMF_NIS.sim0.6@result
p <- my_res.mf %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.mf$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#8491B4FF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) + 
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('NIS GO-MF') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.png',
       width = 3, height = 2)
ggsave(p, filename = 'Results/Figure1/15M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.pdf',
       width = 3, height = 2)

# -- ------------------------------------------------------------------------------------------------------------------
# 10M
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data.nor <- read.csv('Results/Figure1/Files/Mouse.10M_Raw.Counts_CPM.quantiles_scaleRow.csv', row.names = 1)
range(data.nor)
data.Serum <- data.nor[9:16]
data.Con <- data.nor[,c(1:8,17:24)]
data.Serum.T.F <- data.Serum >= 1
data.Serum$Count <- rowSums(data.Serum.T.F)
data.Serum <- data.Serum %>% filter(Count >= 4)
data.Con.T.F <- data.Con >= 1
data.Con$Count <- rowSums(data.Con.T.F)
data.Con <- data.Con %>% filter(Count <= 6)
data.Serum$Pro <- rownames(data.Serum)
data.Con$Pro <- rownames(data.Con)
data.Serum.T <- inner_join(data.Serum, data.Con, by = 'Pro')
write_csv(data.Serum.T, 'Results/Figure1/Files/Mouse.10M_Raw.Counts_CPM.quantiles_scaleRow_SerumSpecial.csv')

data.Urine <- data.nor[17:24]
data.Con <- data.nor[1:16]
data.Urine.T.F <- data.Urine >= 1
data.Urine$Count <- rowSums(data.Urine.T.F)
data.Urine <- data.Urine %>% filter(Count >= 4)
data.Con.T.F <- data.Con >= 1
data.Con$Count <- rowSums(data.Con.T.F)
data.Con <- data.Con %>% filter(Count <= 6)
data.Urine$Pro <- rownames(data.Urine)
data.Con$Pro <- rownames(data.Con)
data.Urine.T <- inner_join(data.Urine, data.Con, by = 'Pro')
write_csv(data.Urine.T, 'Results/Figure1/Files/Mouse.10M_Raw.Counts_CPM.quantiles_scaleRow_UrineSpecial.csv')

data.NIS <- data.nor[1:8]
data.Con <- data.nor[9:24]
data.NIS.T.F <- data.NIS >= 1
data.NIS$Count <- rowSums(data.NIS.T.F)
data.NIS <- data.NIS %>% filter(Count >= 4)
data.Con.T.F <- data.Con >= 1
data.Con$Count <- rowSums(data.Con.T.F)
data.Con <- data.Con %>% filter(Count <= 6)
data.NIS$Pro <- rownames(data.NIS)
data.Con$Pro <- rownames(data.Con)
data.NIS.T <- inner_join(data.NIS, data.Con, by = 'Pro')
write_csv(data.NIS.T, 'Results/Figure1/Files/Mouse.10M_Raw.Counts_CPM.quantiles_scaleRow_NISSpecial.csv')

myPro <- c(data.Serum.T$Pro, data.Urine.T$Pro, data.NIS.T$Pro)
myPro %>% unique()
data.nor <- data.nor %>% filter(rownames(data.nor) %in% myPro)
coldata <- names(data.nor) %>% as.data.frame()
names(coldata) <- 'Samples'
coldata$Group <- ifelse(str_detect(coldata$Samples,'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
coldata$Type <- ifelse(str_detect(coldata$Samples,'s'), 'Serum',
                       ifelse(str_detect(coldata$Samples,'u'), 'Urine', 'NIS')) %>% factor(levels = c('Serum', 'Urine', 'NIS'))
rownames(coldata) <- coldata$Samples
coldata <- coldata[-1]
# show_col(pal_igv(alpha = 0.8)(25))
ann_colors = list(Group = c('AD' ='#CC0000FF', 'Con' ='#6699FFFF'),
                  Type = c('Serum' = '#D595A7CC', 'Urine' = '#F0E685CC', 'NIS' = '#6BD76BCC')) #注意ann_colors是列表
pheatmap(data.nor, color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 6,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = T,
         fontsize_number = 4, fontsize_row = 6, fontsize_col = 8, border_color = NA,
         cluster_rows = T, cluster_cols = F, treeheight_row = 40,
         annotation_col = coldata, annotation_colors = ann_colors,
         cutree_rows = 3,
         gaps_col = c(8,16),
         filename = 'Results/Figure1/Mouse.10M_Raw.Counts_CPM.quantiles_pheatmap.TreeRow.RowName_TypeSpecial.png') # 保存，自动调整纸张大小
pheatmap(data.nor, color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 6,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = T,
         fontsize_number = 4, fontsize_row = 6, fontsize_col = 8, border_color = NA,
         cluster_rows = T, cluster_cols = F, treeheight_row = 40,
         annotation_col = coldata, annotation_colors = ann_colors,
         cutree_rows = 3,
         gaps_col = c(8,16),
         filename = 'Results/Figure1/Mouse.10M_Raw.Counts_CPM.quantiles_pheatmap.TreeRow.RowName_TypeSpecial.pdf') # 保存，自动调整纸张大小
# -- ------------------------------------------------------------------------------------------------------------------
# 组织特异性的蛋白进行功能分析
# -- ------------------------------------------------------------------------------------------------------------------
library(clusterProfiler)
library(enrichplot) #用于可视化的包
library(org.Hs.eg.db) #小鼠用org.Mm.eg.db
# -- ------------------------------------------------------------------------------------------------------------------
# Serum
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/Figure1/Files/Mouse.10M_Raw.Counts_CPM.quantiles_scaleRow_SerumSpecial.csv')
gene <- data$Pro
gene <-  str_replace(gene, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2') %>%
  str_replace('PAR2', 'F2RL1')
genelist <- bitr(gene, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
genelist$SYMBOL %>% unique()
myGene <- genelist$ENTREZID %>% unique()
# -- ------------------------------------------------------------------------
# GO富集分析：
# -- ------------------------------------------------------------------------
# GO over-representation analysis
# -- ------------------------------------------------------------------------
# Biological Process
egoBP_Serum <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoBP_Serum <- clusterProfiler::filter(egoBP_Serum, p.adjust < 0.05)
dim(egoBP_Serum@result)
head(egoBP_Serum@result$Description)
write_csv(egoBP_Serum@result, 'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-BP.Sig.csv')
egoBP_Serum.sim0.8 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.8)
egoBP_Serum.sim0.7 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.7)
egoBP_Serum.sim0.6 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.6)
egoBP_Serum.sim0.5 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.5)
egoBP_Serum.sim0.4 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.4)
egoBP_Serum.sim0.3 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.3)
egoBP_Serum.sim0.8@result$Description %>% head(n=10)
egoBP_Serum.sim0.7@result$Description %>% head(n=10)
egoBP_Serum.sim0.6@result$Description %>% head(n=10)
egoBP_Serum.sim0.5@result$Description %>% head(n=10)
egoBP_Serum.sim0.4@result$Description %>% head(n=10)
egoBP_Serum.sim0.3@result$Description %>% head(n=10)
write.csv(egoBP_Serum.sim0.6@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-BP.Sig.sim0.6.csv', row.names = F)
write.csv(egoBP_Serum.sim0.5@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-BP.Sig.sim0.5.csv', row.names = F)
# Cellular Component
egoCC_Serum <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoCC_Serum <- clusterProfiler::filter(egoCC_Serum, p.adjust < 0.05)
dim(egoCC_Serum@result)
head(egoCC_Serum@result$Description)
write_csv(egoCC_Serum@result, 'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-CC.Sig.csv')
egoCC_Serum.sim0.8 <- clusterProfiler::simplify(egoCC_Serum, cutoff = 0.8)
egoCC_Serum.sim0.7 <- clusterProfiler::simplify(egoCC_Serum, cutoff = 0.7)
egoCC_Serum.sim0.6 <- clusterProfiler::simplify(egoCC_Serum, cutoff = 0.6)
egoCC_Serum.sim0.5 <- clusterProfiler::simplify(egoCC_Serum, cutoff = 0.5)
egoCC_Serum.sim0.8@result$Description %>% head(n=10)
egoCC_Serum.sim0.7@result$Description %>% head(n=10)
egoCC_Serum.sim0.6@result$Description %>% head(n=10)
egoCC_Serum.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC_Serum.sim0.6@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-CC.Sig.sim0.6.csv', row.names = F)
write.csv(egoCC_Serum.sim0.5@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-CC.Sig.sim0.5.csv', row.names = F)
# Molecular Function
egoMF_Serum <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoMF_Serum <- clusterProfiler::filter(egoMF_Serum, p.adjust < 0.05)
dim(egoMF_Serum@result)
head(egoMF_Serum@result$Description)
write_csv(egoMF_Serum@result, 'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-MF.Sig.csv')
egoMF_Serum_5 <- gofilter(egoMF_Serum, level = 5)
egoMF_Serum_5@result$Description %>% head(n=10)
write.csv(egoMF_Serum_5@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-MF.Sig.Level5.csv', row.names = F)
egoMF_Serum.sim0.8 <- clusterProfiler::simplify(egoMF_Serum, cutoff = 0.8)
egoMF_Serum.sim0.7 <- clusterProfiler::simplify(egoMF_Serum, cutoff = 0.7)
egoMF_Serum.sim0.6 <- clusterProfiler::simplify(egoMF_Serum, cutoff = 0.6)
egoMF_Serum.sim0.5 <- clusterProfiler::simplify(egoMF_Serum, cutoff = 0.5)
egoMF_Serum.sim0.8@result$Description %>% head(n=10)
egoMF_Serum.sim0.7@result$Description %>% head(n=10)
egoMF_Serum.sim0.6@result$Description %>% head(n=10)
egoMF_Serum.sim0.5@result$Description %>% head(n=10)
write.csv(egoMF_Serum.sim0.6@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-MF.Sig.sim0.6.csv', row.names = F)
write.csv(egoMF_Serum.sim0.5@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-MF.Sig.sim0.5.csv', row.names = F)
# 保存数据，后续作图用
save(data, myGene, 
     egoBP_Serum, egoCC_Serum, egoMF_Serum, 
     egoBP_Serum.sim0.6, egoCC_Serum.sim0.6, egoMF_Serum.sim0.6,
     egoBP_Serum.sim0.5, egoCC_Serum.sim0.5, egoMF_Serum.sim0.5,
     file = "Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO.RData")
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/Figure1/Files/Mouse.10M_Raw.Counts_CPM.quantiles_scaleRow_UrineSpecial.csv')
gene <- data$Pro
gene <-  str_replace(gene, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2') %>%
  str_replace('PAR2', 'F2RL1')
genelist <- bitr(gene, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
genelist$SYMBOL %>% unique()
myGene <- genelist$ENTREZID %>% unique()
# -- ------------------------------------------------------------------------
# GO富集分析：
# -- ------------------------------------------------------------------------
# GO over-representation analysis
# -- ------------------------------------------------------------------------
# Biological Process
egoBP_Urine <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoBP_Urine <- clusterProfiler::filter(egoBP_Urine, p.adjust < 0.05)
dim(egoBP_Urine@result)
head(egoBP_Urine@result$Description)
write_csv(egoBP_Urine@result, 'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-BP.Sig.csv')
egoBP_Urine.sim0.8 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.8)
egoBP_Urine.sim0.7 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.7)
egoBP_Urine.sim0.6 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.6)
egoBP_Urine.sim0.5 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.5)
egoBP_Urine.sim0.4 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.4)
egoBP_Urine.sim0.3 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.3)
egoBP_Urine.sim0.8@result$Description %>% head(n=10)
egoBP_Urine.sim0.7@result$Description %>% head(n=10)
egoBP_Urine.sim0.6@result$Description %>% head(n=10)
egoBP_Urine.sim0.5@result$Description %>% head(n=10)
egoBP_Urine.sim0.4@result$Description %>% head(n=10)
egoBP_Urine.sim0.3@result$Description %>% head(n=10)
write.csv(egoBP_Urine.sim0.6@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-BP.Sig.sim0.6.csv', row.names = F)
write.csv(egoBP_Urine.sim0.5@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-BP.Sig.sim0.5.csv', row.names = F)
# Cellular Component
egoCC_Urine <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoCC_Urine <- clusterProfiler::filter(egoCC_Urine, p.adjust < 0.05)
dim(egoCC_Urine@result)
head(egoCC_Urine@result$Description)
write_csv(egoCC_Urine@result, 'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-CC.Sig.csv')
egoCC_Urine.sim0.8 <- clusterProfiler::simplify(egoCC_Urine, cutoff = 0.8)
egoCC_Urine.sim0.7 <- clusterProfiler::simplify(egoCC_Urine, cutoff = 0.7)
egoCC_Urine.sim0.6 <- clusterProfiler::simplify(egoCC_Urine, cutoff = 0.6)
egoCC_Urine.sim0.5 <- clusterProfiler::simplify(egoCC_Urine, cutoff = 0.5)
egoCC_Urine.sim0.8@result$Description %>% head(n=10)
egoCC_Urine.sim0.7@result$Description %>% head(n=10)
egoCC_Urine.sim0.6@result$Description %>% head(n=10)
egoCC_Urine.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC_Urine.sim0.6@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-CC.Sig.sim0.6.csv', row.names = F)
write.csv(egoCC_Urine.sim0.5@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-CC.Sig.sim0.5.csv', row.names = F)
# Molecular Function
egoMF_Urine <- enrichGO(gene         = myGene,
                        # universe      = gene_All,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        minGSSize     = 5,
                        readable      = TRUE)
egoMF_Urine <- clusterProfiler::filter(egoMF_Urine, p.adjust < 0.05)
dim(egoMF_Urine@result)
head(egoMF_Urine@result$Description)
write_csv(egoMF_Urine@result, 'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-MF.Sig.csv')
egoMF_Urine_5 <- gofilter(egoMF_Urine, level = 5)
egoMF_Urine_5@result$Description %>% head(n=10)
write.csv(egoMF_Urine_5@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-MF.Sig.Level5.csv', row.names = F)
egoMF_Urine.sim0.8 <- clusterProfiler::simplify(egoMF_Urine, cutoff = 0.8)
egoMF_Urine.sim0.7 <- clusterProfiler::simplify(egoMF_Urine, cutoff = 0.7)
egoMF_Urine.sim0.6 <- clusterProfiler::simplify(egoMF_Urine, cutoff = 0.6)
egoMF_Urine.sim0.5 <- clusterProfiler::simplify(egoMF_Urine, cutoff = 0.5)
egoMF_Urine.sim0.8@result$Description %>% head(n=10)
egoMF_Urine.sim0.7@result$Description %>% head(n=10)
egoMF_Urine.sim0.6@result$Description %>% head(n=10)
egoMF_Urine.sim0.5@result$Description %>% head(n=10)
write.csv(egoMF_Urine.sim0.6@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-MF.Sig.sim0.6.csv', row.names = F)
write.csv(egoMF_Urine.sim0.5@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-MF.Sig.sim0.5.csv', row.names = F)
# 保存数据，后续作图用
save(data, myGene, 
     egoBP_Urine, egoCC_Urine, egoMF_Urine, 
     egoBP_Urine.sim0.6, egoCC_Urine.sim0.6, egoMF_Urine.sim0.6,
     egoBP_Urine.sim0.5, egoCC_Urine.sim0.5, egoMF_Urine.sim0.5,
     file = "Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO.RData")
# -- ------------------------------------------------------------------------------------------------------------------
# NIS
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/Figure1/Files/Mouse.10M_Raw.Counts_CPM.quantiles_scaleRow_NISSpecial.csv')
gene <- data$Pro
gene <-  str_replace(gene, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2') %>%
  str_replace('PAR2', 'F2RL1')
genelist <- bitr(gene, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
genelist$SYMBOL %>% unique()
myGene <- genelist$ENTREZID %>% unique()
# -- ------------------------------------------------------------------------
# GO富集分析：
# -- ------------------------------------------------------------------------
# GO over-representation analysis
# -- ------------------------------------------------------------------------
# Biological Process
egoBP_NIS <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoBP_NIS <- clusterProfiler::filter(egoBP_NIS, p.adjust < 0.05)
dim(egoBP_NIS@result)
head(egoBP_NIS@result$Description)
write_csv(egoBP_NIS@result, 'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-BP.Sig.csv')
egoBP_NIS.sim0.8 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.8)
egoBP_NIS.sim0.7 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.7)
egoBP_NIS.sim0.6 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.6)
egoBP_NIS.sim0.5 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.5)
egoBP_NIS.sim0.4 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.4)
egoBP_NIS.sim0.3 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.3)
egoBP_NIS.sim0.8@result$Description %>% head(n=10)
egoBP_NIS.sim0.7@result$Description %>% head(n=10)
egoBP_NIS.sim0.6@result$Description %>% head(n=10)
egoBP_NIS.sim0.5@result$Description %>% head(n=10)
egoBP_NIS.sim0.4@result$Description %>% head(n=10)
egoBP_NIS.sim0.3@result$Description %>% head(n=10)
write.csv(egoBP_NIS.sim0.6@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-BP.Sig.sim0.6.csv', row.names = F)
write.csv(egoBP_NIS.sim0.5@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-BP.Sig.sim0.5.csv', row.names = F)
# Cellular Component
egoCC_NIS <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoCC_NIS <- clusterProfiler::filter(egoCC_NIS, p.adjust < 0.05)
dim(egoCC_NIS@result)
head(egoCC_NIS@result$Description)
write_csv(egoCC_NIS@result, 'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-CC.Sig.csv')
egoCC_NIS.sim0.8 <- clusterProfiler::simplify(egoCC_NIS, cutoff = 0.8)
egoCC_NIS.sim0.7 <- clusterProfiler::simplify(egoCC_NIS, cutoff = 0.7)
egoCC_NIS.sim0.6 <- clusterProfiler::simplify(egoCC_NIS, cutoff = 0.6)
egoCC_NIS.sim0.5 <- clusterProfiler::simplify(egoCC_NIS, cutoff = 0.5)
egoCC_NIS.sim0.8@result$Description %>% head(n=10)
egoCC_NIS.sim0.7@result$Description %>% head(n=10)
egoCC_NIS.sim0.6@result$Description %>% head(n=10)
egoCC_NIS.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC_NIS.sim0.6@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-CC.Sig.sim0.6.csv', row.names = F)
write.csv(egoCC_NIS.sim0.5@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-CC.Sig.sim0.5.csv', row.names = F)
# Molecular Function
egoMF_NIS <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoMF_NIS <- clusterProfiler::filter(egoMF_NIS, p.adjust < 0.05)
dim(egoMF_NIS@result)
head(egoMF_NIS@result$Description)
write_csv(egoMF_NIS@result, 'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-MF.Sig.csv')
egoMF_NIS_5 <- gofilter(egoMF_NIS, level = 5)
egoMF_NIS_5@result$Description %>% head(n=10)
write.csv(egoMF_NIS_5@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-MF.Sig.Level5.csv', row.names = F)
egoMF_NIS.sim0.8 <- clusterProfiler::simplify(egoMF_NIS, cutoff = 0.8)
egoMF_NIS.sim0.7 <- clusterProfiler::simplify(egoMF_NIS, cutoff = 0.7)
egoMF_NIS.sim0.6 <- clusterProfiler::simplify(egoMF_NIS, cutoff = 0.6)
egoMF_NIS.sim0.5 <- clusterProfiler::simplify(egoMF_NIS, cutoff = 0.5)
egoMF_NIS.sim0.8@result$Description %>% head(n=10)
egoMF_NIS.sim0.7@result$Description %>% head(n=10)
egoMF_NIS.sim0.6@result$Description %>% head(n=10)
egoMF_NIS.sim0.5@result$Description %>% head(n=10)
write.csv(egoMF_NIS.sim0.6@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-MF.Sig.sim0.6.csv', row.names = F)
write.csv(egoMF_NIS.sim0.5@result, 
          'Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-MF.Sig.sim0.5.csv', row.names = F)
# 保存数据，后续作图用
save(data, myGene, 
     egoBP_NIS, egoCC_NIS, egoMF_NIS, 
     egoBP_NIS.sim0.6, egoCC_NIS.sim0.6, egoMF_NIS.sim0.6,
     egoBP_NIS.sim0.5, egoCC_NIS.sim0.5, egoMF_NIS.sim0.5,
     file = "Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO.RData")
# -- ------------------------------------------------------------------------
# GO over-representation analysis可视化
# -- ------------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db) #小鼠用org.Mm.eg.db
library(enrichplot) #用于可视化的包
library(ggsci) # 期刊配色，单个离散配色模块中最多只有51种颜色
library(scales) # show_col
library(patchwork) # 拼接图
library(ggthemes) # 主题设置
# -- ------------------------------------------------------------------------------------------------------------------
# Serum
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO.RData')
# ggplot2画图
show_col(pal_npg(alpha = 0.8)(10))
my_res.bp <- egoBP_Serum.sim0.6@result
p <- my_res.bp %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.bp$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#F39B7FFF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) +
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Serum GO-BP') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.png',
       width = 5.6, height = 2)
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.pdf',
       width = 5.6, height = 2)
my_res.cc <- egoCC_Serum.sim0.6@result
p <- my_res.cc %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.cc$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#4DBBD5FF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) + 
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Serum GO-CC') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)
my_res.mf <- egoMF_Serum.sim0.6@result
p <- my_res.mf %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.mf$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#8491B4FF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) + 
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Serum GO-MF') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_SerumSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO.RData')
# ggplot2画图
show_col(pal_npg(alpha = 0.8)(10))
my_res.bp <- egoBP_Urine.sim0.6@result
p <- my_res.bp %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.bp$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#F39B7FFF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) +
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Urine GO-BP') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.png',
       width = 5.6, height = 2)
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.pdf',
       width = 5.6, height = 2)
my_res.cc <- egoCC_Urine@result
p <- my_res.cc %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.cc$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#4DBBD5FF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) + 
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Urine GO-CC') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-CC.Sig Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-CC.Sig Top5.barplot.pdf',
       width = 4, height = 2)
my_res.mf <- egoMF_Urine@result
p <- my_res.mf %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.mf$p.adjust[4]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#8491B4FF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) + 
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Urine GO-MF') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-MF.Sig.barplot.png',
       width = 3.6, height = 1.8)
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO-MF.Sig.barplot.pdf',
       width = 3.6, height = 1.8)
# -- ------------------------------------------------------------------------------------------------------------------
# NIS
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure1/Files/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO.RData')
# ggplot2画图
show_col(pal_npg(alpha = 0.8)(10))
my_res.bp <- egoBP_NIS.sim0.6@result
p <- my_res.bp %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.bp$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#F39B7FFF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) +
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('NIS GO-BP') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.png',
       width = 5.2, height = 2)
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.pdf',
       width = 5.2, height = 2)
my_res.cc <- egoCC_NIS.sim0.6@result
p <- my_res.cc %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.cc$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#4DBBD5FF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) + 
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('NIS GO-CC') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)
my_res.mf <- egoMF_NIS.sim0.6@result
p <- my_res.mf %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.mf$p.adjust[5]) %>% 
  mutate(logP = -log(p.adjust, 10)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#8491B4FF') +
  geom_vline(xintercept = -log(0.05, 10),
             color = "white",
             linetype = 2,
             size = 0.5) + 
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('NIS GO-MF') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure1/10M_Raw.Counts_CPM.quantiles_NISSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)


# -- ------------------------------------------------------------------------------------------------------------------
# 3.标准化的表达值PCA分析
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse) # 处理数据神包
library(FactoMineR) # PCA分析
library(factoextra) # PCA分析可视化
# -- ------------------------------------------------------------------------------------------------------------------
# 15M
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/Mouse.15M_Raw.Counts_CPM.quantiles.csv', row.names = 1)
data$var <- apply(data, 1, var) # 对每一行求方差
range(data$var)
data <- data %>% filter(var > 0) %>% dplyr::select(!var) # 将表达量全相等的Pro去除
data.pca <- data %>% t() %>% as.data.frame() # PCA分析数据为行Sample，列是变量(维度，成分)
data.pca$Group <- ifelse(str_detect(rownames(data.pca),'D'), 'AD', 'Con') %>% factor(levels = c('AD','Con'))
data.pca$Type <- ifelse(str_detect(rownames(data.pca),'s'), 'Serum',
                        ifelse(str_detect(rownames(data.pca),'u'), 'Urine', 'NIS')) %>% 
  factor(levels = c('Serum', 'Urine', 'NIS'))
data.pca$Color <- str_c(data.pca$Type, data.pca$Group, sep = '.') %>% 
  factor(levels = c('Serum.AD','Serum.Con','Urine.AD','Urine.Con','NIS.AD','NIS.Con'))
data.pca <- data.pca %>% relocate(Group, Type, Color)
pca <- PCA(data.pca[,-c(1:3)], graph = F)
p <- fviz_pca_ind(pca,
                  mean.point = F, # 去除分组的中心点，否则每个群中间会有一个比较大的点
                  # label = "ind", # 展示每一个样本的标签
                  geom.ind = "point", # 样本采用点来展示
                  repel = T, # 避免重叠
                  pointsize = 3,      # 点的大小
                  pointshape = 21,    # 点的形状
                  fill.ind = data.pca$Color, # 根据样本类型来着色
                  alpha.ind = 0.8,
                  palette = c('#AE1F63CC','#D595A7CC','#CC9900CC','#F0E685CC','#6BD76BCC','#749B58CC'), # 3个组设置3种颜色
                  title = '15M',
                  addEllipses = F, # 添加边界线，默认为椭圆
                  ellipse.level = 0.95) + theme_bw() + theme(legend.title=element_blank())
ggsave(p, filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_PCA.png', width = 4.2, height = 3)
ggsave(p, filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_PCA.pdf', width = 4.2, height = 3)
# -- ------------------------------------------------------------------------------------------------------------------
# 10M
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/Mouse.10M_Raw.Counts_CPM.quantiles.csv', row.names = 1)
data$var <- apply(data, 1, var) # 对每一行求方差
range(data$var)
data <- data %>% filter(var > 0) %>% dplyr::select(!var) # 将表达量全相等的Pro去除
data.pca <- data %>% t() %>% as.data.frame() # PCA分析数据为行Sample，列是变量(维度，成分)
data.pca$Group <- ifelse(str_detect(rownames(data.pca),'D'), 'AD', 'Con') %>% factor(levels = c('AD','Con'))
data.pca$Type <- ifelse(str_detect(rownames(data.pca),'s'), 'Serum',
                        ifelse(str_detect(rownames(data.pca),'u'), 'Urine', 'NIS')) %>% 
  factor(levels = c('Serum', 'Urine', 'NIS'))
data.pca$Color <- str_c(data.pca$Type, data.pca$Group, sep = '.') %>% 
  factor(levels = c('Serum.AD','Serum.Con','Urine.AD','Urine.Con','NIS.AD','NIS.Con'))
data.pca <- data.pca %>% relocate(Group, Type, Color)
pca <- PCA(data.pca[,-c(1:3)], graph = F)
p <- fviz_pca_ind(pca,
                  mean.point = F, # 去除分组的中心点，否则每个群中间会有一个比较大的点
                  # label = "ind", # 展示每一个样本的标签
                  geom.ind = "point", # 样本采用点来展示
                  repel = T, # 避免重叠
                  pointsize = 3,      # 点的大小
                  pointshape = 21,    # 点的形状
                  fill.ind = data.pca$Color, # 根据样本类型来着色
                  alpha.ind = 0.8,
                  palette = c('#AE1F63CC','#D595A7CC','#CC9900CC','#F0E685CC','#6BD76BCC','#749B58CC'), # 3个组设置3种颜色
                  title = '10M',
                  addEllipses = F, # 添加边界线，默认为椭圆
                  ellipse.level = 0.95) + theme_bw() + theme(legend.title=element_blank())
ggsave(p, filename = 'Results/Figure1/Mouse.10M_Raw.Counts_CPM.quantiles_PCA.png', width = 4.2, height = 3)
ggsave(p, filename = 'Results/Figure1/Mouse.10M_Raw.Counts_CPM.quantiles_PCA.pdf', width = 4.2, height = 3)
# -- ------------------------------------------------------------------------------------------------------------------
# 3M
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/Mouse_Raw.Counts.csv', row.names = 1)
coldata <- names(data) %>% as.data.frame()
names(coldata) <- 'Samples'
coldata$Group <- ifelse(str_detect(coldata$Samples,'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
coldata$Type <- ifelse(str_detect(coldata$Samples,'s'), 'Serum',
                       ifelse(str_detect(coldata$Samples,'u'), 'Urine', 'NIS')) %>% factor(levels = c('Serum', 'Urine', 'NIS'))
ID <- str_split(coldata$Samples, '_') %>% as.data.frame() %>% t() %>% as.data.frame()
coldata$Month <- ifelse(str_detect(ID$V1,'3'), '3M',
                        ifelse(str_detect(ID$V1,'10'), '10M','15M')) %>% factor(levels = c('3M', '10M', '15M'))
rownames(coldata) <- coldata$Samples
coldata <- coldata[-1]
# 3M
coldata.sub <- coldata %>% filter(Month == '3M') %>% dplyr::select(!Month)
dt <- data %>% dplyr::select(rownames(coldata.sub))
dt$mean <- rowMeans(dt)
range(dt$mean)
dt <- dt %>% filter(mean > 0) %>% dplyr::select(!mean) # 将表达量为0的Pro去除
dt.nor <- dt %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(dt.nor) <- rownames(dt)
names(dt.nor) <- names(dt)
write.csv(dt.nor, 'Results/Mouse.3M_Raw.Counts_CPM.quantiles.csv')
data <- read.csv('Results/Mouse.3M_Raw.Counts_CPM.quantiles.csv', row.names = 1)
data$var <- apply(data, 1, var) # 对每一行求方差
range(data$var)
data <- data %>% filter(var > 0) %>% dplyr::select(!var) # 将表达量全相等的Pro去除
data.pca <- data %>% t() %>% as.data.frame() # PCA分析数据为行Sample，列是变量(维度，成分)
data.pca$Group <- ifelse(str_detect(rownames(data.pca),'D'), 'AD', 'Con') %>% factor(levels = c('AD','Con'))
data.pca$Type <- ifelse(str_detect(rownames(data.pca),'s'), 'Serum',
                        ifelse(str_detect(rownames(data.pca),'u'), 'Urine', 'NIS')) %>% 
  factor(levels = c('Serum', 'Urine', 'NIS'))
data.pca$Color <- str_c(data.pca$Type, data.pca$Group, sep = '.') %>% 
  factor(levels = c('Serum.AD','Serum.Con','Urine.AD','Urine.Con','NIS.AD','NIS.Con'))
data.pca <- data.pca %>% relocate(Group, Type, Color)
pca <- PCA(data.pca[,-c(1:3)], graph = F)
p <- fviz_pca_ind(pca,
                  mean.point = F, # 去除分组的中心点，否则每个群中间会有一个比较大的点
                  # label = "ind", # 展示每一个样本的标签
                  geom.ind = "point", # 样本采用点来展示
                  repel = T, # 避免重叠
                  pointsize = 3,      # 点的大小
                  pointshape = 21,    # 点的形状
                  fill.ind = data.pca$Color, # 根据样本类型来着色
                  alpha.ind = 0.8,
                  palette = c('#AE1F63CC','#D595A7CC','#CC9900CC','#F0E685CC','#6BD76BCC','#749B58CC'), # 3个组设置3种颜色
                  title = '3M',
                  addEllipses = F, # 添加边界线，默认为椭圆
                  ellipse.level = 0.95) + theme_bw() + theme(legend.title=element_blank())
ggsave(p, filename = 'Results/Figure1/Mouse.3M_Raw.Counts_CPM.quantiles_PCA.png', width = 4.2, height = 3)
ggsave(p, filename = 'Results/Figure1/Mouse.3M_Raw.Counts_CPM.quantiles_PCA.pdf', width = 4.2, height = 3)
