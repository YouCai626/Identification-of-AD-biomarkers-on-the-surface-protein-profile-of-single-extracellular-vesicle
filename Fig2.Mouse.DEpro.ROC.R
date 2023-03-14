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
# 1.将每一种样本每一个时间点的Raw.Counts单独拧出来，因为要单独分析
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse) # 包含了ggplot2等包
library(edgeR) # 用于cpm标准化
library(preprocessCore) # normalize.quantiles函数进行标准化
rm(list=ls())
gc()
data <- read.csv('Results/Mouse_Raw.Counts.csv', row.names = 1)
sampleID <- read_csv('Mouse.PT_SampleID_SampleName.csv', show_col_types = F)
# MouseSerum
sub.ID <- sampleID %>% filter(Source == 'Serum')
sub.data <- data %>% dplyr::select(sub.ID$SampleName)
sub.data$mean <- rowMeans(sub.data)
sub.data <- sub.data %>% filter(mean > 0) %>% dplyr::select(!mean) # 将表达量为0的Pro去除
write_csv(sub.ID, 'Results/MouseSerum_Samples.Annotation.csv')
write.csv(sub.data, 'Results/MouseSerum_Raw.Counts.csv')
sub.data.nor <- sub.data %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(sub.data.nor) <- rownames(sub.data)
names(sub.data.nor) <- names(sub.data)
write.csv(sub.data.nor, 'Results/MouseSerum_Raw.Counts_CPM.quantiles.csv')
# MouseUrine
sub.ID <- sampleID %>% filter(Source == 'Urine')
sub.data <- data %>% dplyr::select(sub.ID$SampleName)
sub.data$mean <- rowMeans(sub.data)
sub.data <- sub.data %>% filter(mean > 0) %>% dplyr::select(!mean) # 将表达量为0的Pro去除
write_csv(sub.ID, 'Results/MouseUrine_Samples.Annotation.csv')
write.csv(sub.data, 'Results/MouseUrine_Raw.Counts.csv')
sub.data.nor <- sub.data %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(sub.data.nor) <- rownames(sub.data)
names(sub.data.nor) <- names(sub.data)
write.csv(sub.data.nor, 'Results/MouseUrine_Raw.Counts_CPM.quantiles.csv')
# MouseNIS
sub.ID <- sampleID %>% filter(Source == 'NIS')
sub.data <- data %>% dplyr::select(sub.ID$SampleName)
sub.data$mean <- rowMeans(sub.data)
sub.data <- sub.data %>% filter(mean > 0) %>% dplyr::select(!mean) # 将表达量为0的Pro去除
write_csv(sub.ID, 'Results/MouseNIS_Samples.Annotation.csv')
write.csv(sub.data, 'Results/MouseNIS_Raw.Counts.csv')
sub.data.nor <- sub.data %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(sub.data.nor) <- rownames(sub.data)
names(sub.data.nor) <- names(sub.data)
write.csv(sub.data.nor, 'Results/MouseNIS_Raw.Counts_CPM.quantiles.csv')

# -- ------------------------------------------------------------------------------------------------------------------
# 2.每一种体液每个时间点的差异表达蛋白
# -- ------------------------------------------------------------------------------------------------------------------
# 因为我们的数据已经经过标准化处理，所以用limma来进行差异表达分析
library(tidyverse) # 包含了ggplot2等包
library(limma)
# -- ------------------------------------------------------------------------------------------------------------------
# Serum
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
split <- c('3s','10s','15s')
month <- c('3M','10M','15M')
data <- read.csv('Results/MouseSerum_Raw.Counts_CPM.quantiles.csv', row.names = 1)
for (j in 1:3) {
  data.sub <- data %>% dplyr::select(contains(split[j]))
  data.sub.log <- log(data.sub + 0.01, 2)
  group <- ifelse(str_detect(colnames(data.sub), 'D'), 'AD', 'Con') %>% factor()
  design <- model.matrix(~0 + group)
  colnames(design) <- c("AD", "Con")
  rownames(design) <- colnames(data.sub)
  contrast.matrix <- makeContrasts(AD - Con, levels = design)
  fit <- lmFit(data.sub.log, design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  DEpro <- topTable(fit, number = Inf, adjust.method = "fdr")
  DEpro <- DEpro %>% arrange(P.Value)
  fn <- str_c('Results/Figure2/Files/MouseSerum.', month[j], '_ADvsCon_DEpro_Limma.csv')
  write.csv(DEpro, fn)
}
# 用wilcox.test进行差异表达分析
for (j in 1:3) {
  data.sub <- data %>% dplyr::select(contains(split[j]))
  conditions <- ifelse(str_detect(colnames(data.sub), 'D'), 'AD', 'Con') %>% factor()
  pvalues <- sapply(1:nrow(data.sub),function(i){
    temp <- cbind.data.frame(gene = as.numeric(t(data.sub[i,])), conditions)
    p = wilcox.test(gene~conditions, temp)$p.value
    return(p)
  })
  fdr <- p.adjust(pvalues, method = "fdr")
  conditionsLevel <- levels(conditions)
  data.AD <- data.sub[,c(which(conditions==conditionsLevel[1]))]
  data.Con <- data.sub[,c(which(conditions==conditionsLevel[2]))]
  foldChanges <- log2(rowMeans(data.AD)/rowMeans(data.Con))
  outRst <- data.frame(log2FC = foldChanges, pValues = pvalues, FDR = fdr)
  outRst <- na.omit(outRst) %>% arrange(pValues)
  fn <- str_c('Results/Figure2/Files/MouseSerum.', month[j], '_ADvsCon_DEpro_wilcox.test.csv')
  write.csv(outRst, fn)
}
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
split <- c('3u','10u','15u')
month <- c('3M','10M','15M')
data <- read.csv('Results/MouseUrine_Raw.Counts_CPM.quantiles.csv', row.names = 1)
for (j in 1:3) {
  data.sub <- data %>% dplyr::select(contains(split[j]))
  data.sub.log <- log(data.sub + 0.01, 2)
  group <- ifelse(str_detect(colnames(data.sub), 'D'), 'AD', 'Con') %>% factor()
  design <- model.matrix(~0 + group)
  colnames(design) <- c("AD", "Con")
  rownames(design) <- colnames(data.sub)
  contrast.matrix <- makeContrasts(AD - Con, levels = design)
  fit <- lmFit(data.sub.log, design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  DEpro <- topTable(fit, number = Inf, adjust.method = "fdr")
  DEpro <- DEpro %>% arrange(P.Value)
  fn <- str_c('Results/Figure2/Files/MouseUrine.', month[j], '_ADvsCon_DEpro_Limma.csv')
  write.csv(DEpro, fn)
}
# 用wilcox.test进行差异表达分析
for (j in 1:3) {
  data.sub <- data %>% dplyr::select(contains(split[j]))
  conditions <- ifelse(str_detect(colnames(data.sub), 'D'), 'AD', 'Con') %>% factor()
  pvalues <- sapply(1:nrow(data.sub),function(i){
    temp <- cbind.data.frame(gene = as.numeric(t(data.sub[i,])), conditions)
    p = wilcox.test(gene~conditions, temp)$p.value
    return(p)
  })
  fdr <- p.adjust(pvalues, method = "fdr")
  conditionsLevel <- levels(conditions)
  data.AD <- data.sub[,c(which(conditions==conditionsLevel[1]))]
  data.Con <- data.sub[,c(which(conditions==conditionsLevel[2]))]
  foldChanges <- log2(rowMeans(data.AD)/rowMeans(data.Con))
  outRst <- data.frame(log2FC = foldChanges, pValues = pvalues, FDR = fdr)
  outRst <- na.omit(outRst) %>% arrange(pValues)
  fn <- str_c('Results/Figure2/Files/MouseUrine.', month[j], '_ADvsCon_DEpro_wilcox.test.csv')
  write.csv(outRst, fn)
}

# -- ------------------------------------------------------------------------------------------------------------------
# NIS
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
split <- c('3n','10n','15n')
month <- c('3M','10M','15M')
data <- read.csv('Results/MouseNIS_Raw.Counts_CPM.quantiles.csv', row.names = 1)
for (j in 1:3) {
  data.sub <- data %>% dplyr::select(contains(split[j]))
  data.sub.log <- log(data.sub + 0.01, 2)
  group <- ifelse(str_detect(colnames(data.sub), 'D'), 'AD', 'Con') %>% factor()
  design <- model.matrix(~0 + group)
  colnames(design) <- c("AD", "Con")
  rownames(design) <- colnames(data.sub)
  contrast.matrix <- makeContrasts(AD - Con, levels = design)
  fit <- lmFit(data.sub.log, design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  DEpro <- topTable(fit, number = Inf, adjust.method = "fdr")
  DEpro <- DEpro %>% arrange(P.Value)
  fn <- str_c('Results/Figure2/Files/MouseNIS.', month[j], '_ADvsCon_DEpro_Limma.csv')
  write.csv(DEpro, fn)
}
# 用wilcox.test进行差异表达分析
for (j in 1:3) {
  data.sub <- data %>% dplyr::select(contains(split[j]))
  conditions <- ifelse(str_detect(colnames(data.sub), 'D'), 'AD', 'Con') %>% factor()
  pvalues <- sapply(1:nrow(data.sub),function(i){
    temp <- cbind.data.frame(gene = as.numeric(t(data.sub[i,])), conditions)
    p = wilcox.test(gene~conditions, temp)$p.value
    return(p)
  })
  fdr <- p.adjust(pvalues, method = "fdr")
  conditionsLevel <- levels(conditions)
  data.AD <- data.sub[,c(which(conditions==conditionsLevel[1]))]
  data.Con <- data.sub[,c(which(conditions==conditionsLevel[2]))]
  foldChanges <- log2(rowMeans(data.AD)/rowMeans(data.Con))
  outRst <- data.frame(log2FC = foldChanges, pValues = pvalues, FDR = fdr)
  outRst <- na.omit(outRst) %>% arrange(pValues)
  fn <- str_c('Results/Figure2/Files/MouseNIS.', month[j], '_ADvsCon_DEpro_wilcox.test.csv')
  write.csv(outRst, fn)
}
# -- ------------------------------------------------------------------------------------------------------------------
# 3.将阳性的DEpro合并在一起
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
# Serum
filenames <- Sys.glob('Results/Figure2/Files/MouseSerum.*_ADvsCon_DEpro_Limma.csv')
data <- read.csv(filenames[1])
data <- data %>% filter(P.Value < 0.05)
month <- filenames %>% str_remove('Results/Figure2/Files/MouseSerum.') %>%
  str_remove('_ADvsCon_DEpro_Limma.csv')
data$Month <- month[1]
for (i in 2:3) {
  dt <- read.csv(filenames[i])
  dt <- dt %>% filter(P.Value < 0.05)
  dt$Month <- month[i]
  data <- rbind(data, dt)
}
write.csv(data, 'Results/Figure2/Files/MouseSerum.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
# Urine
filenames <- Sys.glob('Results/Figure2/Files/MouseUrine.*_ADvsCon_DEpro_Limma.csv')
data <- read.csv(filenames[1])
data <- data %>% filter(P.Value < 0.05)
month <- filenames %>% str_remove('Results/Figure2/Files/MouseUrine.') %>%
  str_remove('_ADvsCon_DEpro_Limma.csv')
data$Month <- month[1]
for (i in 2:3) {
  dt <- read.csv(filenames[i])
  dt <- dt %>% filter(P.Value < 0.05)
  dt$Month <- month[i]
  data <- rbind(data, dt)
}
write.csv(data, 'Results/Figure2/Files/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
# NIS
filenames <- Sys.glob('Results/Figure2/Files/MouseNIS.*_ADvsCon_DEpro_Limma.csv')
data <- read.csv(filenames[1])
data <- data %>% filter(P.Value < 0.05)
month <- filenames %>% str_remove('Results/Figure2/Files/MouseNIS.') %>%
  str_remove('_ADvsCon_DEpro_Limma.csv')
data$Month <- month[1]
for (i in 2:3) {
  dt <- read.csv(filenames[i])
  dt <- dt %>% filter(P.Value < 0.05)
  dt$Month <- month[i]
  data <- rbind(data, dt)
}
write.csv(data, 'Results/Figure2/Files/MouseNIS.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')

rm(list=ls())
gc()
# Serum
filenames <- Sys.glob('Results/Figure2/Files/MouseSerum.*_ADvsCon_DEpro_wilcox.test.csv')
data <- read.csv(filenames[1])
data <- data %>% filter(pValues < 0.05)
month <- filenames %>% str_remove('Results/Figure2/Files/MouseSerum.') %>%
  str_remove('_ADvsCon_DEpro_wilcox.test.csv')
data$Month <- month[1]
for (i in 2:3) {
  dt <- read.csv(filenames[i])
  dt <- dt %>% filter(pValues < 0.05)
  dt$Month <- month[i]
  data <- rbind(data, dt)
}
write.csv(data, 'Results/Figure2/Files/MouseSerum.3-10-15M.ADvsCon_DEpro.Pos_wilcox.test.csv')
# Urine
filenames <- Sys.glob('Results/Figure2/Files/MouseUrine.*_ADvsCon_DEpro_wilcox.test.csv')
data <- read.csv(filenames[1])
data <- data %>% filter(pValues < 0.05)
month <- filenames %>% str_remove('Results/Figure2/Files/MouseUrine.') %>%
  str_remove('_ADvsCon_DEpro_wilcox.test.csv')
data$Month <- month[1]
for (i in 2:3) {
  dt <- read.csv(filenames[i])
  dt <- dt %>% filter(pValues < 0.05)
  dt$Month <- month[i]
  data <- rbind(data, dt)
}
write.csv(data, 'Results/Figure2/Files/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_wilcox.test.csv')
# NIS
filenames <- Sys.glob('Results/Figure2/Files/MouseNIS.*_ADvsCon_DEpro_wilcox.test.csv')
data <- read.csv(filenames[1])
data <- data %>% filter(pValues < 0.05)
month <- filenames %>% str_remove('Results/Figure2/Files/MouseNIS.') %>%
  str_remove('_ADvsCon_DEpro_wilcox.test.csv')
data$Month <- month[1]
for (i in 2:3) {
  dt <- read.csv(filenames[i])
  dt <- dt %>% filter(pValues < 0.05)
  dt$Month <- month[i]
  data <- rbind(data, dt)
}
write.csv(data, 'Results/Figure2/Files/MouseNIS.3-10-15M.ADvsCon_DEpro.Pos_wilcox.test.csv')

# -- ------------------------------------------------------------------------------------------------------------------
# 20230224Tang说将3-10-15M的差异表达基因合并在一起看3种体液有没有Overlap
# 将特异性的和Overlap的分别进行功能富集分析
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data.Serum <- read.csv('Results/Figure2/Files/MouseSerum.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
data.Serum <- data.Serum %>% filter(AveExpr >= 9)
data.Urine <- read.csv('Results/Figure2/Files/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
data.Urine <- data.Urine %>% filter(AveExpr >= 9)
data.NIS <- read.csv('Results/Figure2/Files/MouseNIS.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
data.NIS <- data.NIS %>% filter(AveExpr >= 9)
# 可视化常用包
library(tidyverse) # 包含了ggplot2等包
library(ggsci) # 期刊配色，单个离散配色模块中最多只有51种颜色
library(scales) # show_col(c('#CC0000FF','#6699FFFF'))
library(patchwork) # 拼接图
library(ggbreak) # 坐标轴截断
library(ggthemes) # 主题设置
library(stringr) # 文本处理包，包含在tidyverse中
library(ggsignif) # 计算差异显著性并添加显著性线段
library(ggrepel)

library(VennDiagram)
# 'Serum' = '#D595A7CC', 'Urine' = '#F0E685CC', 'NIS' = '#6BD76BCC'
venn_list <- list(Serum = data.Serum$X, Urine = data.Urine$X,
                  NIS = data.NIS$X)
p <- venn.diagram(x = venn_list, filename = NULL, output=TRUE, resolution = 600, scaled = F,
                  col = "white", alpha = 0.4,
                  fill = c('red','yellow','green'), cex = 2.5,
                  fontfamily = "serif", fontface = "bold",
                  cat.col = 'black', cat.cex = 2.5, cat.default.pos = "outer",
                  cat.fontfamily = "serif", cat.fontface = "bold",
                  cat.dist = c(0.01, 0.01, 0.01),
                  cat.pos = c(0,0,180))
p <- venn.diagram(x = venn_list, filename = NULL, output=TRUE, resolution = 600, scaled = F,
                 col = "white", alpha = 0.8,
                 fill = c('#D595A7FF','#F0E685FF','#6BD76BFF'), cex = 2.5,
                 fontfamily = "serif", fontface = "bold",
                 cat.col = 'black', cat.cex = 2.5, cat.default.pos = "outer",
                 cat.fontfamily = "serif", cat.fontface = "bold",
                 cat.dist = c(0.01, 0.01, 0.01),
                 cat.pos = c(0,0,180))
png("Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma.Overlap.png", height = 600, width = 600)
grid.draw(p)
dev.off()
pdf("Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma.Overlap.pdf", height = 6, width = 6)
grid.draw(p)
dev.off()
inter <- get.venn.partitions(venn_list) %>% as.data.frame()
write.csv(as.matrix(inter), 'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma.Overlap.csv')
serum <- inter$..values..[7] %>% unlist()
urine <- inter$..values..[6] %>% unlist()
nis <- inter$..values..[4] %>% unlist()
over <- inter$..values..[1] %>% unlist()
# -- ------------------------------------------------------------------------------------------------------------------
# 功能注释
library(clusterProfiler)
library(enrichplot) #用于可视化的包
library(org.Hs.eg.db) #小鼠用org.Mm.eg.db
# -- ------------------------------------------------------------------------------------------------------------------
# Serum
# -- ------------------------------------------------------------------------------------------------------------------
serum <-  str_replace(serum, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2') %>%
  str_replace('PAR2', 'F2RL1') %>% str_replace('ITGA4B7', 'ITGA4')
genelist <- bitr(serum, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
genelist$SYMBOL %>% unique()
myGene <- genelist$ENTREZID %>% unique()
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
write_csv(egoBP_Serum@result, 'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-BP.Sig.csv')
egoBP_Serum.sim0.6 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.6)
egoBP_Serum.sim0.5 <- clusterProfiler::simplify(egoBP_Serum, cutoff = 0.5)
egoBP_Serum.sim0.6@result$Description %>% head(n=10)
egoBP_Serum.sim0.5@result$Description %>% head(n=10)
write.csv(egoBP_Serum.sim0.6@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-BP.Sig.sim0.6.csv', row.names = F)
write.csv(egoBP_Serum.sim0.5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-BP.Sig.sim0.5.csv', row.names = F)
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
write_csv(egoCC_Serum@result, 'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-CC.Sig.csv')
egoCC_Serum.sim0.6 <- clusterProfiler::simplify(egoCC_Serum, cutoff = 0.6)
egoCC_Serum.sim0.5 <- clusterProfiler::simplify(egoCC_Serum, cutoff = 0.5)
egoCC_Serum.sim0.6@result$Description %>% head(n=10)
egoCC_Serum.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC_Serum.sim0.6@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-CC.Sig.sim0.6.csv', row.names = F)
write.csv(egoCC_Serum.sim0.5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-CC.Sig.sim0.5.csv', row.names = F)
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
write_csv(egoMF_Serum@result, 'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-MF.Sig.csv')
egoMF_Serum_5 <- gofilter(egoMF_Serum, level = 5)
egoMF_Serum_5@result$Description %>% head(n=10)
write.csv(egoMF_Serum_5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-MF.Sig.Level5.csv', row.names = F)
egoMF_Serum.sim0.8 <- clusterProfiler::simplify(egoMF_Serum, cutoff = 0.8)
egoMF_Serum.sim0.7 <- clusterProfiler::simplify(egoMF_Serum, cutoff = 0.7)
egoMF_Serum.sim0.6 <- clusterProfiler::simplify(egoMF_Serum, cutoff = 0.6)
egoMF_Serum.sim0.5 <- clusterProfiler::simplify(egoMF_Serum, cutoff = 0.5)
egoMF_Serum.sim0.8@result$Description %>% head(n=10)
egoMF_Serum.sim0.7@result$Description %>% head(n=10)
egoMF_Serum.sim0.6@result$Description %>% head(n=10)
egoMF_Serum.sim0.5@result$Description %>% head(n=10)
write.csv(egoMF_Serum.sim0.7@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-MF.Sig.sim0.7.csv', row.names = F)
write.csv(egoMF_Serum.sim0.6@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-MF.Sig.sim0.6.csv', row.names = F)
# 保存数据，后续作图用
save(serum, egoBP_Serum, egoCC_Serum, egoMF_Serum, 
     egoBP_Serum.sim0.6, egoCC_Serum.sim0.6, egoMF_Serum.sim0.6,
     egoBP_Serum.sim0.5, egoCC_Serum.sim0.5, egoMF_Serum.sim0.7,
     file = "Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO.RData")
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
urine <-  str_replace(urine, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2') %>%
  str_replace('PAR2', 'F2RL1') %>% str_replace('ITGA4B7', 'ITGA4')
genelist <- bitr(urine, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
genelist$SYMBOL %>% unique()
myGene <- genelist$ENTREZID %>% unique()
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
write_csv(egoBP_Urine@result, 'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-BP.Sig.csv')
egoBP_Urine.sim0.6 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.6)
egoBP_Urine.sim0.5 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.5)
egoBP_Urine.sim0.6@result$Description %>% head(n=10)
egoBP_Urine.sim0.5@result$Description %>% head(n=10)
write.csv(egoBP_Urine.sim0.6@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-BP.Sig.sim0.6.csv', row.names = F)
write.csv(egoBP_Urine.sim0.5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-BP.Sig.sim0.5.csv', row.names = F)
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
write_csv(egoCC_Urine@result, 'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-CC.Sig.csv')
egoCC_Urine.sim0.6 <- clusterProfiler::simplify(egoCC_Urine, cutoff = 0.6)
egoCC_Urine.sim0.5 <- clusterProfiler::simplify(egoCC_Urine, cutoff = 0.5)
egoCC_Urine.sim0.6@result$Description %>% head(n=10)
egoCC_Urine.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC_Urine.sim0.6@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-CC.Sig.sim0.6.csv', row.names = F)
write.csv(egoCC_Urine.sim0.5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-CC.Sig.sim0.5.csv', row.names = F)
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
write_csv(egoMF_Urine@result, 'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-MF.Sig.csv')
egoMF_Urine_5 <- gofilter(egoMF_Urine, level = 5)
egoMF_Urine_5@result$Description %>% head(n=10)
write.csv(egoMF_Urine_5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-MF.Sig.Level5.csv', row.names = F)
egoMF_Urine.sim0.6 <- clusterProfiler::simplify(egoMF_Urine, cutoff = 0.6)
egoMF_Urine.sim0.5 <- clusterProfiler::simplify(egoMF_Urine, cutoff = 0.5)
egoMF_Urine.sim0.6@result$Description %>% head(n=10)
egoMF_Urine.sim0.5@result$Description %>% head(n=10)
write.csv(egoMF_Urine.sim0.5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-MF.Sig.sim0.5.csv', row.names = F)
write.csv(egoMF_Urine.sim0.6@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-MF.Sig.sim0.6.csv', row.names = F)
# 保存数据，后续作图用
save(urine, egoBP_Urine, egoCC_Urine, egoMF_Urine, 
     egoBP_Urine.sim0.6, egoCC_Urine.sim0.6, egoMF_Urine.sim0.6,
     egoBP_Urine.sim0.5, egoCC_Urine.sim0.5, egoMF_Urine.sim0.5,
     file = "Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO.RData")
# -- ------------------------------------------------------------------------------------------------------------------
# NIS
# -- ------------------------------------------------------------------------------------------------------------------
nis <-  str_replace(nis, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2') %>%
  str_replace('PAR2', 'F2RL1') %>% str_replace('ITGA4B7', 'ITGA4') %>% str_replace('ITGAVB5', 'ITGAV')
genelist <- bitr(nis, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
genelist$SYMBOL %>% unique()
myGene <- genelist$ENTREZID %>% unique()
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
write_csv(egoBP_NIS@result, 'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-BP.Sig.csv')
egoBP_NIS.sim0.6 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.6)
egoBP_NIS.sim0.5 <- clusterProfiler::simplify(egoBP_NIS, cutoff = 0.5)
egoBP_NIS.sim0.6@result$Description %>% head(n=10)
egoBP_NIS.sim0.5@result$Description %>% head(n=10)
write.csv(egoBP_NIS.sim0.6@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-BP.Sig.sim0.6.csv', row.names = F)
write.csv(egoBP_NIS.sim0.5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-BP.Sig.sim0.5.csv', row.names = F)
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
write_csv(egoCC_NIS@result, 'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-CC.Sig.csv')
egoCC_NIS.sim0.6 <- clusterProfiler::simplify(egoCC_NIS, cutoff = 0.6)
egoCC_NIS.sim0.5 <- clusterProfiler::simplify(egoCC_NIS, cutoff = 0.5)
egoCC_NIS.sim0.6@result$Description %>% head(n=10)
egoCC_NIS.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC_NIS.sim0.6@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-CC.Sig.sim0.6.csv', row.names = F)
write.csv(egoCC_NIS.sim0.5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-CC.Sig.sim0.5.csv', row.names = F)
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
write_csv(egoMF_NIS@result, 'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-MF.Sig.csv')
egoMF_NIS_5 <- gofilter(egoMF_NIS, level = 5)
egoMF_NIS_5@result$Description %>% head(n=10)
write.csv(egoMF_NIS_5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-MF.Sig.Level5.csv', row.names = F)
egoMF_NIS.sim0.6 <- clusterProfiler::simplify(egoMF_NIS, cutoff = 0.6)
egoMF_NIS.sim0.5 <- clusterProfiler::simplify(egoMF_NIS, cutoff = 0.5)
egoMF_NIS.sim0.6@result$Description %>% head(n=10)
egoMF_NIS.sim0.5@result$Description %>% head(n=10)
write.csv(egoMF_NIS.sim0.5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-MF.Sig.sim0.5.csv', row.names = F)
write.csv(egoMF_NIS.sim0.6@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-MF.Sig.sim0.6.csv', row.names = F)
# 保存数据，后续作图用
save(nis, egoBP_NIS, egoCC_NIS, egoMF_NIS, 
     egoBP_NIS.sim0.6, egoCC_NIS.sim0.6, egoMF_NIS.sim0.6,
     egoBP_NIS.sim0.5, egoCC_NIS.sim0.5, egoMF_NIS.sim0.5,
     file = "Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO.RData")
# -- ------------------------------------------------------------------------------------------------------------------
# Overlap
# -- ------------------------------------------------------------------------------------------------------------------
over <-  str_replace(over, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2') %>%
  str_replace('PAR2', 'F2RL1') %>% str_replace('ITGA4B7', 'ITGA4') %>% str_replace('ITGAVB5', 'ITGAV')
genelist <- bitr(over, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
genelist$SYMBOL %>% unique()
myGene <- genelist$ENTREZID %>% unique()
# Biological Process
egoBP_over <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoBP_over <- clusterProfiler::filter(egoBP_over, p.adjust < 0.05)
dim(egoBP_over@result)
head(egoBP_over@result$Description)
write_csv(egoBP_over@result, 'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-BP.Sig.csv')
egoBP_over.sim0.6 <- clusterProfiler::simplify(egoBP_over, cutoff = 0.6)
egoBP_over.sim0.5 <- clusterProfiler::simplify(egoBP_over, cutoff = 0.5)
egoBP_over.sim0.6@result$Description %>% head(n=10)
egoBP_over.sim0.5@result$Description %>% head(n=10)
write.csv(egoBP_over.sim0.6@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-BP.Sig.sim0.6.csv', row.names = F)
write.csv(egoBP_over.sim0.5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-BP.Sig.sim0.5.csv', row.names = F)
# Cellular Component
egoCC_over <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoCC_over <- clusterProfiler::filter(egoCC_over, p.adjust < 0.05)
dim(egoCC_over@result)
head(egoCC_over@result$Description)
write_csv(egoCC_over@result, 'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-CC.Sig.csv')
egoCC_over.sim0.6 <- clusterProfiler::simplify(egoCC_over, cutoff = 0.6)
egoCC_over.sim0.5 <- clusterProfiler::simplify(egoCC_over, cutoff = 0.5)
egoCC_over.sim0.6@result$Description %>% head(n=10)
egoCC_over.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC_over.sim0.6@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-CC.Sig.sim0.6.csv', row.names = F)
write.csv(egoCC_over.sim0.5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-CC.Sig.sim0.5.csv', row.names = F)
# Molecular Function
egoMF_over <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoMF_over <- clusterProfiler::filter(egoMF_over, p.adjust < 0.05)
dim(egoMF_over@result)
head(egoMF_over@result$Description)
write_csv(egoMF_over@result, 'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-MF.Sig.csv')
egoMF_over_5 <- gofilter(egoMF_over, level = 5)
egoMF_over_5@result$Description %>% head(n=10)
write.csv(egoMF_over_5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-MF.Sig.Level5.csv', row.names = F)
egoMF_over.sim0.6 <- clusterProfiler::simplify(egoMF_over, cutoff = 0.6)
egoMF_over.sim0.5 <- clusterProfiler::simplify(egoMF_over, cutoff = 0.5)
egoMF_over.sim0.6@result$Description %>% head(n=10)
egoMF_over.sim0.5@result$Description %>% head(n=10)
write.csv(egoMF_over.sim0.5@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-MF.Sig.sim0.5.csv', row.names = F)
write.csv(egoMF_over.sim0.6@result, 
          'Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-MF.Sig.sim0.6.csv', row.names = F)
# 保存数据，后续作图用
save(over, egoBP_over, egoCC_over, egoMF_over, 
     egoBP_over.sim0.6, 
     egoBP_over.sim0.5, 
     file = "Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO.RData")
# -- ------------------------------------------------------------------------
# GO over-representation analysis可视化
# -- ------------------------------------------------------------------------
# Serum
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO.RData')
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
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.png',
       width = 6, height = 2)
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.pdf',
       width = 6, height = 2)
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
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.pdf',
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
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_SerumSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)
# -- ------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO.RData')
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
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.png',
       width = 5, height = 2)
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.pdf',
       width = 5, height = 2)
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
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)
my_res.mf <- egoMF_Urine.sim0.6@result %>% arrange(p.adjust)
p <- my_res.mf[1:5,] %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.mf$p.adjust[5]) %>% 
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
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.png',
       width = 3.6, height = 2)
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_UrineSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.pdf',
       width = 3.6, height = 2)
# -- ------------------------------------------------------------------------
# NIS
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO.RData')
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
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.png',
       width = 6, height = 2)
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-BP.Sig.sim0.6 Top5.barplot.pdf',
       width = 6, height = 2)
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
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.png',
       width = 5, height = 2)
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-CC.Sig.sim0.6 Top5.barplot.pdf',
       width = 5, height = 2)
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
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_NISSpecial_enrichGO-MF.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)
# -- ------------------------------------------------------------------------
# Overlap
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure2/Files/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO.RData')
# ggplot2画图
show_col(pal_npg(alpha = 0.8)(10))
my_res.bp <- egoBP_over.sim0.6@result
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
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Overlap GO-BP') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-BP.Sig.sim0.6 Top5.barplot.png',
       width = 6, height = 2)
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-BP.Sig.sim0.6 Top5.barplot.pdf',
       width = 6, height = 2)
my_res.cc <- egoCC_over@result
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
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Overlap GO-CC') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-CC.Sig Top5.barplot.png',
       width = 5, height = 2)
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-CC.Sig Top5.barplot.pdf',
       width = 5, height = 2)
my_res.mf <- egoMF_over@result
p <- my_res.mf %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.mf$p.adjust[3]) %>% 
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
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('Overlap GO-MF') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-MF.Sig.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure2/Mouse.3-10-15M.ADvsCon_DEpro.Pos_Limma_Overlap_enrichGO-MF.Sig.barplot.pdf',
       width = 4, height = 2)


# -- ------------------------------------------------------------------------------------------------------------------
# 4.可视化
# -- ------------------------------------------------------------------------------------------------------------------
# 差异表达基因的数目统计
library(ggpubr)
library(ggsignif)
# -- ------------------------------------------------------------------------------------------------------------------
rm(list = ls())
gc()
select = dplyr::select
mycol.month <- c('#F39B7FFF','#4DBBD5FF','#DC0000FF')
show_col(mycol.month)
# Serum
deg_list <- read.csv('Results/Figure2/Files/MouseSerum.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
names(deg_list)[1] <- 'Pro'
df <- deg_list %>% 
  mutate(change = case_when(logFC > 0 & P.Value < 0.05~ "Up",
                            logFC < 0 & P.Value < 0.05 ~ "Down",
                            T ~ "No"))
df %>% 
  group_by(Month, change) %>%
  summarise(deg_nmber = n()) %>%
  filter(change != "No") %>%
  ungroup() %>%
  mutate(Month = factor(Month, levels = c("3M","10M","15M")),
         change = factor(change, levels = c("Up","Down"))) %>%
  ggplot(aes(x = Month ,y = deg_nmber ,fill = Month)) + 
  geom_bar(stat= "identity",width = 0.8) + 
  scale_fill_manual(values = mycol.month) +
  ggthemes::theme_few() +
  xlab("Month") +
  ylab("Number of DEpro")+
  theme(axis.text = element_text(color = "black"))+
  facet_wrap(change~.,nrow = 2) +
  coord_flip() +
  guides(fill = "none")
ggsave("Results/Figure2/MouseSerum.DEpro.Pos_Limma.Number.png",width = 3.6,height = 4)
ggsave("Results/Figure2/MouseSerum.DEpro.Pos_Limma.Number.pdf",width = 3.6,height = 4)
# Urine
deg_list <- read.csv('Results/Figure2/Files/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
names(deg_list)[1] <- 'Pro'
df <- deg_list %>% 
  mutate(change = case_when(logFC > 0 & P.Value < 0.05~ "Up",
                            logFC < 0 & P.Value < 0.05 ~ "Down",
                            T ~ "No"))
df %>% 
  group_by(Month, change) %>%
  summarise(deg_nmber = n()) %>%
  filter(change != "No") %>%
  ungroup() %>%
  mutate(Month = factor(Month, levels = c("3M","10M","15M")),
         change = factor(change, levels = c("Up","Down"))) %>%
  ggplot(aes(x = Month ,y = deg_nmber ,fill = Month)) + 
  geom_bar(stat= "identity",width = 0.8) + 
  scale_fill_manual(values = mycol.month) +
  ggthemes::theme_few() +
  xlab("Month") +
  ylab("Number of DEpro")+
  theme(axis.text = element_text(color = "black"))+
  facet_wrap(change~.,nrow = 2) +
  coord_flip() +
  guides(fill = "none")
ggsave("Results/Figure2/MouseUrine.DEpro.Pos_Limma.Number.png",width = 3.6,height = 4)
ggsave("Results/Figure2/MouseUrine.DEpro.Pos_Limma.Number.pdf",width = 3.6,height = 4)
# NIS
deg_list <- read.csv('Results/Figure2/Files/MouseNIS.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
names(deg_list)[1] <- 'Pro'
df <- deg_list %>% 
  mutate(change = case_when(logFC > 0 & P.Value < 0.05~ "Up",
                            logFC < 0 & P.Value < 0.05 ~ "Down",
                            T ~ "No"))
df %>% 
  group_by(Month, change) %>%
  summarise(deg_nmber = n()) %>%
  filter(change != "No") %>%
  ungroup() %>%
  mutate(Month = factor(Month, levels = c("3M","10M","15M")),
         change = factor(change, levels = c("Up","Down"))) %>%
  ggplot(aes(x = Month ,y = deg_nmber ,fill = Month)) + 
  geom_bar(stat= "identity",width = 0.8) + 
  scale_fill_manual(values = mycol.month) +
  ggthemes::theme_few() +
  xlab("Month") +
  ylab("Number of DEpro")+
  theme(axis.text = element_text(color = "black"))+
  facet_wrap(change~.,nrow = 2) +
  coord_flip() +
  guides(fill = "none")
ggsave("Results/Figure2/MouseNIS.DEpro.Pos_Limma.Number.png",width = 3.6,height = 4)
ggsave("Results/Figure2/MouseNIS.DEpro.Pos_Limma.Number.pdf",width = 3.6,height = 4)
# -- ------------------------------------------------------------------------------------------------------------------
# plot valcano----
# -- ------------------------------------------------------------------------------------------------------------------
# Serum
# -- ------------------------------------------------------------------------------------------------------------------
rm(list = ls())
gc()
select = dplyr::select
data <- read.csv('Results/Figure2/Files/MouseSerum.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
data <- data %>% filter(X != 'IgA') %>% filter(X != 'IgM')
data$Month <- factor(data$Month, levels = c('3M','10M','15M'))
data$Sig <- ifelse(data$logFC > 0, 'Up', 'Down')
data$Sig <- data$Sig %>% factor(levels = c('Up','Down'))
p <- ggplot(data, aes(x = Month, y = logFC, color = Sig)) + theme_bw() +
  geom_jitter(alpha = 0.8, size = 2, width = 0.4) + xlab('') + 
  scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
  theme(axis.text.x = element_text(colour = 'black', size = rel(1.5)),
        axis.text.y = element_text(colour = 'black', size = rel(1.5)),
        axis.title.y = element_text(colour = 'black', size = rel(1.6)))
p  
#根据图p中log2FC区间确定背景柱长度：
tmp <- data %>%
  group_by(Month) %>%
  slice_max(logFC)
bar_up <- data.frame(x = tmp$Month,
                     y = tmp$logFC + 0.1)
tmp <- data %>%
  group_by(Month) %>%
  slice_min(logFC)
bar_bottom <- data.frame(x = tmp$Month,
                         y = tmp$logFC-0.1)
#绘制背景柱：
p1 <- ggplot()+
  geom_col(data = bar_up,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",
           alpha = 0.6,
           width = 0.8) +
  geom_col(data = bar_bottom,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",
           alpha = 0.6,
           width = 0.8) +
  geom_jitter(data = data, 
              mapping = aes(x = Month, y = logFC, color = Sig),
              alpha = 0.8, 
              size = 2, 
              width = 0.4) + 
  theme_bw() + xlab('') + ggtitle('MouseSerum') +
  scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
  theme(axis.text.x = element_text(colour = 'black', size = rel(1.5)),
        axis.text.y = element_text(colour = 'black', size = rel(1.5)),
        axis.title.y = element_text(colour = 'black', size = rel(1.6)))
p1
#添加X轴的cluster色块标签：
df.box <- data.frame(x=c('3M','10M','15M'),
                     y=0,
                     label=c('3M','10M','15M'))
mycol.month <- c('#F39B7FFF','#4DBBD5FF','#DC0000FF')
p2 <- p1 + 
  geom_tile(data = df.box,
            aes(x=x,y=y),
            height = 0.16,
            width = 0.8,
            color = "black",
            fill = mycol.month,
            alpha = 0.8,
            show.legend = F) +
  xlab("Months")+
  ylab("log2FoldChange")+
  geom_text(data = df.box,
            aes(x = x,y = y,label=label),
            size = 4,
            color = "white")
p2
p3  <- p2 +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 12,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 0.5),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10,color = "black",face = "bold"),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0)#,legend.text = element_text(size = 15)
  )
p3
ggsave('Results/Figure2/MouseSerum.3-10-15M.ADvsCon_DEpro.Pos_Limma.png', p3, width = 4, height = 8)
ggsave('Results/Figure2/MouseSerum.3-10-15M.ADvsCon_DEpro.Pos_Limma.pdf', p3, width = 4, height = 8)
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list = ls())
gc()
select = dplyr::select
data <- read.csv('Results/Figure2/Files/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
data <- data %>% filter(AveExpr >= 9)
data$Month <- factor(data$Month, levels = c('3M','10M','15M'))
data$Sig <- ifelse(data$logFC > 0, 'Up', 'Down')
data$Sig <- data$Sig %>% factor(levels = c('Up','Down'))
p <- ggplot(data, aes(x = Month, y = logFC, color = Sig)) + theme_bw() +
  geom_jitter(alpha = 0.8, size = 2, width = 0.4) + xlab('') + 
  scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
  theme(axis.text.x = element_text(colour = 'black', size = rel(1.5)),
        axis.text.y = element_text(colour = 'black', size = rel(1.5)),
        axis.title.y = element_text(colour = 'black', size = rel(1.6)))
p  
#根据图p中log2FC区间确定背景柱长度：
tmp <- data %>%
  group_by(Month) %>%
  slice_max(logFC)
bar_up <- data.frame(x = tmp$Month,
                     y = tmp$logFC + 0.1)
tmp <- data %>%
  group_by(Month) %>%
  slice_min(logFC)
bar_bottom <- data.frame(x = tmp$Month,
                         y = tmp$logFC-0.1)
#绘制背景柱：
p1 <- ggplot()+
  geom_col(data = bar_up,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",
           alpha = 0.6,
           width = 0.8) +
  geom_col(data = bar_bottom,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",
           alpha = 0.6,
           width = 0.8) +
  geom_jitter(data = data, 
              mapping = aes(x = Month, y = logFC, color = Sig),
              alpha = 0.8, 
              size = 2, 
              width = 0.4) + 
  theme_bw() + xlab('') + ggtitle('MouseUrine') +
  scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
  theme(axis.text.x = element_text(colour = 'black', size = rel(1.5)),
        axis.text.y = element_text(colour = 'black', size = rel(1.5)),
        axis.title.y = element_text(colour = 'black', size = rel(1.6)))
p1
#添加X轴的cluster色块标签：
df.box <- data.frame(x=c('3M','10M','15M'),
                     y=0,
                     label=c('3M','10M','15M'))
mycol.month <- c('#F39B7FFF','#4DBBD5FF','#DC0000FF')
p2 <- p1 + 
  geom_tile(data = df.box,
            aes(x=x,y=y),
            height = 0.16,
            width = 0.8,
            color = "black",
            fill = mycol.month,
            alpha = 0.8,
            show.legend = F) +
  xlab("Months")+
  ylab("log2FoldChange")+
  geom_text(data = df.box,
            aes(x = x,y = y,label=label),
            size = 4,
            color = "white")
p2
p3  <- p2 +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 12,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 0.5),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10,color = "black",face = "bold"),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0)#,legend.text = element_text(size = 15)
  )
p3
ggsave('Results/Figure2/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma.png', p3, width = 4, height = 8)
ggsave('Results/Figure2/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma.pdf', p3, width = 4, height = 8)
# -- ------------------------------------------------------------------------------------------------------------------
# 20230303画有趋势线的图
# -- ------------------------------------------------------------------------------------------------------------------
rm(list = ls())
gc()
select = dplyr::select
data <- read.csv('Results/Figure2/Files/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
data <- data %>% filter(AveExpr >= 9)
data$Month <- factor(data$Month, levels = c('3M','10M','15M'))
data$Time <- ifelse(data$Month == '3M', 3, ifelse(data$Month == '10M', 10, 15))
data$Sig <- ifelse(data$logFC > 0, 'Up', 'Down')
data$Sig <- data$Sig %>% factor(levels = c('Up','Down'))
# 添加趋势线
data1 <- data %>% filter(Sig == 'Up')
ggplot(data1, aes(x = Time, y = logFC, color = Sig)) + theme_bw() +
  geom_jitter(alpha = 0.7, size = 2, color = '#CC0000FF') + xlab('') + 
  theme(axis.text.x = element_text(colour = 'black', size = rel(1.5)),
        axis.text.y = element_text(colour = 'black', size = rel(1.5)),
        axis.title.y = element_text(colour = 'black', size = rel(1.4))) + 
  geom_hline(aes(yintercept = 0), colour = "grey", linetype = 1, size = 6, alpha = 0.8) +
  geom_vline(aes(xintercept = 6.5), colour = "black", linetype = 2, size = 1, alpha = 0.6) + 
  geom_vline(aes(xintercept = 12.5), colour = "black", linetype = 2, size = 1, alpha = 0.6) +
  # geom_text_repel(aes(label = Label), colour = 'black', size = 3) + 
  stat_smooth(method='lm', formula = y~x, colour='black', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(method = "spearman", digits = 3, size = 4.5, label.y = 1, label.x = 1, colour = 'black') + 
  scale_x_continuous(breaks = c(3,10,15), n.breaks = 3)
ggsave('Results/Figure2/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma_Up 2.png', width = 4, height = 3)
ggsave('Results/Figure2/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma_Up 2.pdf', width = 4, height = 3)
data1 <- data %>% filter(Sig == 'Down')
ggplot(data1, aes(x = Time, y = logFC)) + theme_bw() +
  geom_jitter(alpha = 0.7, size = 2, color = '#6699FFFF') + xlab('') + 
  theme(axis.text.x = element_text(colour = 'black', size = rel(1.5)),
        axis.text.y = element_text(colour = 'black', size = rel(1.5)),
        axis.title.y = element_text(colour = 'black', size = rel(1.4))) + 
  geom_hline(aes(yintercept = 0), colour = "grey", linetype = 1, size = 6, alpha = 0.8) +
  geom_vline(aes(xintercept = 6.5), colour = "black", linetype = 2, size = 1, alpha = 0.6) + 
  geom_vline(aes(xintercept = 12.5), colour = "black", linetype = 2, size = 1, alpha = 0.6) +
  # geom_text_repel(aes(label = Label), colour = 'black', size = 3) + 
  stat_smooth(method='lm', formula = y~x, colour='black', linetype = 2, size = 0.5, alpha = 0.2) +
  stat_cor(method = "spearman", digits = 3, size = 4.5, label.y = -0.85, label.x = 1, colour = 'black') + 
  scale_x_continuous(breaks = c(3,10,15), n.breaks = 3)
ggsave('Results/Figure2/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma_Down 2.png', width = 4, height = 3)
ggsave('Results/Figure2/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma_Down 2.pdf', width = 4, height = 3)
# -- ------------------------------------------------------------------------------------------------------------------
# 20230224Tang说将Urine的3-10-15M的差异表达基因分别进行功能富集分析
# -- ------------------------------------------------------------------------------------------------------------------
rm(list = ls())
gc()
select = dplyr::select
data <- read.csv('Results/Figure2/Files/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
data <- data %>% filter(AveExpr >= 9)
M3 <- data %>% filter(Month == '3M')
M10 <- data %>% filter(Month == '10M')
M15 <- data %>% filter(Month == '15M')
# -- ------------------------------------------------------------------------------------------------------------------
# 15M
# -- ------------------------------------------------------------------------------------------------------------------
gene <-  str_replace(M15$X, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2') %>%
  str_replace('PAR2', 'F2RL1')
genelist <- bitr(gene, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
genelist$SYMBOL %>% unique()
myGene <- genelist$ENTREZID %>% unique()
# Biological Process
egoBP_15M <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoBP_15M <- clusterProfiler::filter(egoBP_15M, p.adjust < 0.05)
dim(egoBP_15M@result)
head(egoBP_15M@result$Description)
write_csv(egoBP_15M@result, 'Results/Figure2/Files/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.csv')
egoBP_15M.sim0.6 <- clusterProfiler::simplify(egoBP_15M, cutoff = 0.6)
egoBP_15M.sim0.5 <- clusterProfiler::simplify(egoBP_15M, cutoff = 0.5)
egoBP_15M.sim0.6@result$Description %>% head(n=10)
egoBP_15M.sim0.5@result$Description %>% head(n=10)
write.csv(egoBP_15M.sim0.6@result, 
          'Results/Figure2/Files/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.sim0.6.csv', row.names = F)
write.csv(egoBP_15M.sim0.5@result, 
          'Results/Figure2/Files/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.sim0.5.csv', row.names = F)
# Cellular Component
egoCC_15M <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoCC_15M <- clusterProfiler::filter(egoCC_15M, p.adjust < 0.05)
dim(egoCC_15M@result)
head(egoCC_15M@result$Description)
write_csv(egoCC_15M@result, 'Results/Figure2/Files/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.csv')
egoCC_15M.sim0.6 <- clusterProfiler::simplify(egoCC_15M, cutoff = 0.6)
egoCC_15M.sim0.5 <- clusterProfiler::simplify(egoCC_15M, cutoff = 0.5)
egoCC_15M.sim0.6@result$Description %>% head(n=10)
egoCC_15M.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC_15M.sim0.6@result, 
          'Results/Figure2/Files/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.sim0.6.csv', row.names = F)
write.csv(egoCC_15M.sim0.5@result, 
          'Results/Figure2/Files/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.sim0.5.csv', row.names = F)
# Molecular Function
egoMF_15M <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoMF_15M <- clusterProfiler::filter(egoMF_15M, p.adjust < 0.05)
dim(egoMF_15M@result)
head(egoMF_15M@result$Description)
write_csv(egoMF_15M@result, 'Results/Figure2/Files/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.csv')
egoMF_15M.sim0.6 <- clusterProfiler::simplify(egoMF_15M, cutoff = 0.6)
egoMF_15M.sim0.5 <- clusterProfiler::simplify(egoMF_15M, cutoff = 0.5)
egoMF_15M.sim0.6@result$Description %>% head(n=10)
egoMF_15M.sim0.5@result$Description %>% head(n=10)
write.csv(egoMF_15M.sim0.6@result, 
          'Results/Figure2/Files/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.sim0.6.csv', row.names = F)
write.csv(egoMF_15M.sim0.5@result, 
          'Results/Figure2/Files/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.sim0.5.csv', row.names = F)
# 保存数据，后续作图用
save(data, myGene, 
     egoBP_15M, egoCC_15M, egoMF_15M, 
     egoBP_15M.sim0.6, egoCC_15M.sim0.6, egoMF_15M.sim0.6,
     egoBP_15M.sim0.5, egoCC_15M.sim0.5, egoMF_15M.sim0.5,
     file = "Results/Figure2/Files/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO.RData")
# -- ------------------------------------------------------------------------------------------------------------------
# 10M
# -- ------------------------------------------------------------------------------------------------------------------
gene <-  str_replace(M10$X, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2') %>%
  str_replace('PAR2', 'F2RL1')
genelist <- bitr(gene, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
genelist$SYMBOL %>% unique()
myGene <- genelist$ENTREZID %>% unique()
# Biological Process
egoBP_10M <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoBP_10M <- clusterProfiler::filter(egoBP_10M, p.adjust < 0.05)
dim(egoBP_10M@result)
head(egoBP_10M@result$Description)
write_csv(egoBP_10M@result, 'Results/Figure2/Files/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.csv')
egoBP_10M.sim0.6 <- clusterProfiler::simplify(egoBP_10M, cutoff = 0.6)
egoBP_10M.sim0.5 <- clusterProfiler::simplify(egoBP_10M, cutoff = 0.5)
egoBP_10M.sim0.6@result$Description %>% head(n=10)
egoBP_10M.sim0.5@result$Description %>% head(n=10)
write.csv(egoBP_10M.sim0.6@result, 
          'Results/Figure2/Files/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.sim0.6.csv', row.names = F)
write.csv(egoBP_10M.sim0.5@result, 
          'Results/Figure2/Files/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.sim0.5.csv', row.names = F)
# Cellular Component
egoCC_10M <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoCC_10M <- clusterProfiler::filter(egoCC_10M, p.adjust < 0.05)
dim(egoCC_10M@result)
head(egoCC_10M@result$Description)
write_csv(egoCC_10M@result, 'Results/Figure2/Files/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.csv')
egoCC_10M.sim0.6 <- clusterProfiler::simplify(egoCC_10M, cutoff = 0.6)
egoCC_10M.sim0.5 <- clusterProfiler::simplify(egoCC_10M, cutoff = 0.5)
egoCC_10M.sim0.6@result$Description %>% head(n=10)
egoCC_10M.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC_10M.sim0.6@result, 
          'Results/Figure2/Files/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.sim0.6.csv', row.names = F)
write.csv(egoCC_10M.sim0.5@result, 
          'Results/Figure2/Files/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.sim0.5.csv', row.names = F)
# Molecular Function
egoMF_10M <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoMF_10M <- clusterProfiler::filter(egoMF_10M, p.adjust < 0.05)
dim(egoMF_10M@result)
head(egoMF_10M@result$Description)
write_csv(egoMF_10M@result, 'Results/Figure2/Files/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.csv')
egoMF_10M.sim0.6 <- clusterProfiler::simplify(egoMF_10M, cutoff = 0.6)
egoMF_10M.sim0.5 <- clusterProfiler::simplify(egoMF_10M, cutoff = 0.5)
egoMF_10M.sim0.6@result$Description %>% head(n=10)
egoMF_10M.sim0.5@result$Description %>% head(n=10)
write.csv(egoMF_10M.sim0.6@result, 
          'Results/Figure2/Files/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.sim0.6.csv', row.names = F)
write.csv(egoMF_10M.sim0.5@result, 
          'Results/Figure2/Files/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.sim0.5.csv', row.names = F)
# 保存数据，后续作图用
save(data, myGene, 
     egoBP_10M, egoCC_10M, egoMF_10M, 
     egoBP_10M.sim0.6, egoCC_10M.sim0.6, egoMF_10M.sim0.6,
     egoBP_10M.sim0.5, egoCC_10M.sim0.5, egoMF_10M.sim0.5,
     file = "Results/Figure2/Files/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO.RData")
# -- ------------------------------------------------------------------------------------------------------------------
# 3M
# -- ------------------------------------------------------------------------------------------------------------------
gene <-  str_replace(M3$X, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2') %>%
  str_replace('PAR2', 'F2RL1')
genelist <- bitr(gene, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
genelist$SYMBOL %>% unique()
myGene <- genelist$ENTREZID %>% unique()
# Biological Process
egoBP_3M <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoBP_3M <- clusterProfiler::filter(egoBP_3M, p.adjust < 0.05)
dim(egoBP_3M@result)
head(egoBP_3M@result$Description)
write_csv(egoBP_3M@result, 'Results/Figure2/Files/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.csv')
egoBP_3M.sim0.6 <- clusterProfiler::simplify(egoBP_3M, cutoff = 0.6)
egoBP_3M.sim0.5 <- clusterProfiler::simplify(egoBP_3M, cutoff = 0.5)
egoBP_3M.sim0.6@result$Description %>% head(n=10)
egoBP_3M.sim0.5@result$Description %>% head(n=10)
write.csv(egoBP_3M.sim0.6@result, 
          'Results/Figure2/Files/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.sim0.6.csv', row.names = F)
write.csv(egoBP_3M.sim0.5@result, 
          'Results/Figure2/Files/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.sim0.5.csv', row.names = F)
# Cellular Component
egoCC_3M <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoCC_3M <- clusterProfiler::filter(egoCC_3M, p.adjust < 0.05)
dim(egoCC_3M@result)
head(egoCC_3M@result$Description)
write_csv(egoCC_3M@result, 'Results/Figure2/Files/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.csv')
egoCC_3M.sim0.6 <- clusterProfiler::simplify(egoCC_3M, cutoff = 0.6)
egoCC_3M.sim0.5 <- clusterProfiler::simplify(egoCC_3M, cutoff = 0.5)
egoCC_3M.sim0.6@result$Description %>% head(n=10)
egoCC_3M.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC_3M.sim0.6@result, 
          'Results/Figure2/Files/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.sim0.6.csv', row.names = F)
write.csv(egoCC_3M.sim0.5@result, 
          'Results/Figure2/Files/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.sim0.5.csv', row.names = F)
# Molecular Function
egoMF_3M <- enrichGO(gene         = myGene,
                      # universe      = gene_All,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 5,
                      readable      = TRUE)
egoMF_3M <- clusterProfiler::filter(egoMF_3M, p.adjust < 0.05)
dim(egoMF_3M@result)
head(egoMF_3M@result$Description)
write_csv(egoMF_3M@result, 'Results/Figure2/Files/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.csv')
egoMF_3M.sim0.6 <- clusterProfiler::simplify(egoMF_3M, cutoff = 0.6)
egoMF_3M.sim0.5 <- clusterProfiler::simplify(egoMF_3M, cutoff = 0.5)
egoMF_3M.sim0.6@result$Description %>% head(n=10)
egoMF_3M.sim0.5@result$Description %>% head(n=10)
write.csv(egoMF_3M.sim0.6@result, 
          'Results/Figure2/Files/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.sim0.6.csv', row.names = F)
write.csv(egoMF_3M.sim0.5@result, 
          'Results/Figure2/Files/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.sim0.5.csv', row.names = F)
# 保存数据，后续作图用
save(data, myGene, 
     egoBP_3M, egoCC_3M, egoMF_3M, 
     egoBP_3M.sim0.6, egoCC_3M.sim0.6, egoMF_3M.sim0.6,
     egoBP_3M.sim0.5, egoCC_3M.sim0.5, egoMF_3M.sim0.5,
     file = "Results/Figure2/Files/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO.RData")
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
# 15M
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure2/Files/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO.RData')
# ggplot2画图
show_col(pal_npg(alpha = 0.8)(10))
my_res.bp <- egoBP_15M.sim0.6@result
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
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('15M GO-BP') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure2/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.sim0.6 Top5.barplot.png',
       width = 5.6, height = 2)
ggsave(p, filename = 'Results/Figure2/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.sim0.6 Top5.barplot.pdf',
       width = 5.6, height = 2)
my_res.cc <- egoCC_15M.sim0.6@result
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
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('15M GO-CC') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure2/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure2/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)
my_res.mf <- egoMF_15M.sim0.6@result
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
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('15M GO-MF') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure2/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.sim0.6 Top5.barplot.png',
       width = 5.2, height = 2)
ggsave(p, filename = 'Results/Figure2/MouseUrine.15M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.sim0.6 Top5.barplot.pdf',
       width = 5.2, height = 2)
# -- ------------------------------------------------------------------------------------------------------------------
# 10M
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure2/Files/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO.RData')
# ggplot2画图
show_col(pal_npg(alpha = 0.8)(10))
my_res.bp <- egoBP_10M.sim0.6@result %>% arrange(p.adjust)
p <- my_res.bp[c(1,4,5,6,8,9),] %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.bp$p.adjust[9]) %>% 
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
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('10M GO-BP') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure2/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.sim0.6 Tang.barplot.png',
       width = 6, height = 2)
ggsave(p, filename = 'Results/Figure2/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.sim0.6 Tang.barplot.pdf',
       width = 6, height = 2)
my_res.cc <- egoCC_10M.sim0.6@result
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
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('10M GO-CC') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure2/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.sim0.6 Top5.barplot.png',
       width = 3.6, height = 2)
ggsave(p, filename = 'Results/Figure2/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.sim0.6 Top5.barplot.pdf',
       width = 3.6, height = 2)
my_res.mf <- egoMF_10M.sim0.6@result
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
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('10M GO-MF') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure2/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.sim0.6 Top5.barplot.png',
       width = 3.6, height = 2)
ggsave(p, filename = 'Results/Figure2/MouseUrine.10M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.sim0.6 Top5.barplot.pdf',
       width = 3.6, height = 2)
# -- ------------------------------------------------------------------------------------------------------------------
# 3M
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
load('Results/Figure2/Files/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO.RData')
# ggplot2画图
show_col(pal_npg(alpha = 0.8)(10))
my_res.bp <- egoBP_3M.sim0.6@result %>% arrange(p.adjust)
p <- my_res.bp[1:5,] %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.bp$p.adjust[5]) %>% 
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
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('3M GO-BP') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure2/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.sim0.6 Top5.barplot.png',
       width = 6, height = 2)
ggsave(p, filename = 'Results/Figure2/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-BP.Sig.sim0.6 Top5.barplot.pdf',
       width = 6, height = 2)
my_res.cc <- egoCC_3M.sim0.6@result %>% arrange(p.adjust)
p <- my_res.cc[1:5,] %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.cc$p.adjust[5]) %>% 
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
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('3M GO-CC') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure2/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.sim0.6 Top5.barplot.png',
       width = 5, height = 2)
ggsave(p, filename = 'Results/Figure2/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-CC.Sig.sim0.6 Top5.barplot.pdf',
       width = 5, height = 2)
my_res.mf <- egoMF_3M.sim0.6@result %>% arrange(p.adjust)
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
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('3M GO-MF') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'Results/Figure2/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.sim0.6 Top5.barplot.png',
       width = 4, height = 2)
ggsave(p, filename = 'Results/Figure2/MouseUrine.3M.ADvsCon_DEpro.Pos_Limma_enrichGO-MF.Sig.sim0.6 Top5.barplot.pdf',
       width = 4, height = 2)


# -- ------------------------------------------------------------------------------------------------------------------
# NIS
# -- ------------------------------------------------------------------------------------------------------------------
rm(list = ls())
gc()
select = dplyr::select
data <- read.csv('Results/Figure2/Files/MouseNIS.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
data <- data %>% filter(X != 'IgA') %>% filter(X != 'IgM')
data$Month <- factor(data$Month, levels = c('3M','10M','15M'))
data$Sig <- ifelse(data$logFC > 0, 'Up', 'Down')
data$Sig <- data$Sig %>% factor(levels = c('Up','Down'))
p <- ggplot(data, aes(x = Month, y = logFC, color = Sig)) + theme_bw() +
  geom_jitter(alpha = 0.8, size = 2, width = 0.4) + xlab('') + 
  scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
  theme(axis.text.x = element_text(colour = 'black', size = rel(1.5)),
        axis.text.y = element_text(colour = 'black', size = rel(1.5)),
        axis.title.y = element_text(colour = 'black', size = rel(1.6)))
p  
#根据图p中log2FC区间确定背景柱长度：
tmp <- data %>%
  group_by(Month) %>%
  slice_max(logFC)
bar_up <- data.frame(x = tmp$Month,
                     y = tmp$logFC + 0.1)
tmp <- data %>%
  group_by(Month) %>%
  slice_min(logFC)
bar_bottom <- data.frame(x = tmp$Month,
                         y = tmp$logFC-0.1)
#绘制背景柱：
p1 <- ggplot()+
  geom_col(data = bar_up,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",
           alpha = 0.6,
           width = 0.8) +
  geom_col(data = bar_bottom,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",
           alpha = 0.6,
           width = 0.8) +
  geom_jitter(data = data, 
              mapping = aes(x = Month, y = logFC, color = Sig),
              alpha = 0.8, 
              size = 2, 
              width = 0.4) + 
  theme_bw() + xlab('') + ggtitle('MouseNIS') +
  scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
  theme(axis.text.x = element_text(colour = 'black', size = rel(1.5)),
        axis.text.y = element_text(colour = 'black', size = rel(1.5)),
        axis.title.y = element_text(colour = 'black', size = rel(1.6)))
p1
#添加X轴的cluster色块标签：
df.box <- data.frame(x=c('3M','10M','15M'),
                     y=0,
                     label=c('3M','10M','15M'))
mycol.month <- c('#F39B7FFF','#4DBBD5FF','#DC0000FF')
p2 <- p1 + 
  geom_tile(data = df.box,
            aes(x=x,y=y),
            height = 0.16,
            width = 0.8,
            color = "black",
            fill = mycol.month,
            alpha = 0.8,
            show.legend = F) +
  xlab("Months")+
  ylab("log2FoldChange")+
  geom_text(data = df.box,
            aes(x = x,y = y,label=label),
            size = 4,
            color = "white")
p2
p3  <- p2 +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 12,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 0.5),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10,color = "black",face = "bold"),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0)#,legend.text = element_text(size = 15)
  )
p3
ggsave('Results/Figure2/MouseNIS.3-10-15M.ADvsCon_DEpro.Pos_Limma.png', p3, width = 4, height = 8)
ggsave('Results/Figure2/MouseNIS.3-10-15M.ADvsCon_DEpro.Pos_Limma.pdf', p3, width = 4, height = 8)

# -- ------------------------------------------------------------------------------------------------------------------
# 5.将每一种样本的不同时间点的AD和Con合在一起比较，只考虑AD vs Con，进行ROC分析
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

# 可视化常用包
library(tidyverse) # 包含了ggplot2等包
library(ggsci) # 期刊配色，单个离散配色模块中最多只有51种颜色
library(scales) # show_col
library(patchwork) # 拼接图
library(ggbreak) # 坐标轴截断
library(ggsignif) # 计算差异显著性并添加显著性线段
library(ggthemes) # 主题设置
library(stringr) # 文本处理包，包含在tidyverse中
library(ggrepel)
library(pROC)

# -- ------------------------------------------------------------------------------------------------------------------
# MouseSerum
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/MouseSerum_Raw.Counts_CPM.quantiles.csv', row.names = 1) #行为基因，列为样本
roc.data <- data %>% t() %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'AD'), 'AD', 'Con') %>% factor()
roc.data <- relocate(roc.data, Group)
auc <- c()
z <- c()
pvalue <- c()
for (j in 1:(ncol(roc.data)-1)){
  rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('Con','AD'), direction = '<', data = roc.data)
  auc[j] <- rocPlot$auc
  se <- sqrt(pROC::var(rocPlot))
  b <- auc[j] - 0.5
  z[j] <- (b / se)
  pvalue[j] <- 2 * pt(-abs(z[j]), df=Inf)
  if (pvalue[j] < 0.05 & auc[j] >= 0.7) {
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    filename5 <- paste('Results/Figure2/MouseSerum.', proteinname, '.png', sep = '')
    png(filename5, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)  
    dev.off()
    filename6 <- paste('Results/Figure2/MouseSerum.', proteinname, '.pdf', sep = '')
    pdf(filename6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)  
    dev.off()
  }
  if (pvalue[j] < 0.05 & auc[j] <= 0.3) {
    rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('Con','AD'), direction = '>', data = roc.data)
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    fn6 <- paste('Results/Figure2/MouseSerum.', proteinname, '.png', sep = '')
    png(fn6, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "yellow",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure2/MouseSerum.', proteinname, '.pdf', sep = '')
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
auc_pvalue <- auc_pvalue %>% arrange(pvalue)
write.csv(auc_pvalue, file = 'Results/Figure2/Files/MouseSerum_AUC_pvalue.csv')
auc_pvalue <- filter(auc_pvalue, pvalue < 0.05)
write.csv(auc_pvalue, file = 'Results/Figure2/Files/MouseSerum_AUC_pvalue.Pos.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# MouseUrine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/MouseUrine_Raw.Counts_CPM.quantiles.csv', row.names = 1) #行为基因，列为样本
roc.data <- data %>% t() %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'AD'), 'AD', 'Con') %>% factor()
roc.data <- relocate(roc.data, Group)
auc <- c()
z <- c()
pvalue <- c()
for (j in 1:(ncol(roc.data)-1)){
  rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('Con','AD'), direction = '<', data = roc.data)
  auc[j] <- rocPlot$auc
  se <- sqrt(pROC::var(rocPlot))
  b <- auc[j] - 0.5
  z[j] <- (b / se)
  pvalue[j] <- 2 * pt(-abs(z[j]), df=Inf)
  if (pvalue[j] < 0.05 & auc[j] >= 0.7) {
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    filename5 <- paste('Results/Figure2/MouseUrine.', proteinname, '.png', sep = '')
    png(filename5, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)  
    dev.off()
    filename6 <- paste('Results/Figure2/MouseUrine.', proteinname, '.pdf', sep = '')
    pdf(filename6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)  
    dev.off()
  }
  if (pvalue[j] < 0.05 & auc[j] <= 0.3) {
    rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('Con','AD'), direction = '>', data = roc.data)
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    fn6 <- paste('Results/Figure2/MouseUrine.', proteinname, '.png', sep = '')
    png(fn6, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "yellow",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure2/MouseUrine.', proteinname, '.pdf', sep = '')
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
auc_pvalue <- auc_pvalue %>% arrange(pvalue)
write.csv(auc_pvalue, file = 'Results/Figure2/Files/MouseUrine_AUC_pvalue.csv')
auc_pvalue <- filter(auc_pvalue, pvalue < 0.05)
write.csv(auc_pvalue, file = 'Results/Figure2/Files/MouseUrine_AUC_pvalue.Pos.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# MouseNIS
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/MouseNIS_Raw.Counts_CPM.quantiles.csv', row.names = 1) #行为基因，列为样本
roc.data <- data %>% t() %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'AD'), 'AD', 'Con') %>% factor()
roc.data <- relocate(roc.data, Group)
auc <- c()
z <- c()
pvalue <- c()
for (j in 1:(ncol(roc.data)-1)){
  rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('Con','AD'), direction = '<', data = roc.data)
  auc[j] <- rocPlot$auc
  se <- sqrt(pROC::var(rocPlot))
  b <- auc[j] - 0.5
  z[j] <- (b / se)
  pvalue[j] <- 2 * pt(-abs(z[j]), df=Inf)
  if (pvalue[j] < 0.05 & auc[j] >= 0.7) {
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    filename5 <- paste('Results/Figure2/MouseNIS.', proteinname, '.png', sep = '')
    png(filename5, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)  
    dev.off()
    filename6 <- paste('Results/Figure2/MouseNIS.', proteinname, '.pdf', sep = '')
    pdf(filename6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)  
    dev.off()
  }
  if (pvalue[j] < 0.05 & auc[j] <= 0.3) {
    rocPlot <- roc(roc.data[,1] ~ roc.data[,(j + 1)], levels = c('Con','AD'), direction = '>', data = roc.data)
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    fn6 <- paste('Results/Figure2/MouseNIS.', proteinname, '.png', sep = '')
    png(fn6, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "yellow",
             print.auc.x=0.7, print.auc.y=0.1,# AUC的坐标为（x，y）
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure2/MouseNIS.', proteinname, '.pdf', sep = '')
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
auc_pvalue <- auc_pvalue %>% arrange(pvalue)
write.csv(auc_pvalue, file = 'Results/Figure2/Files/MouseNIS_AUC_pvalue.csv')
auc_pvalue <- filter(auc_pvalue, pvalue < 0.05)
write.csv(auc_pvalue, file = 'Results/Figure2/Files/MouseNIS_AUC_pvalue.Pos.csv')

# -- ------------------------------------------------------------------------------------------------------------------
# 6.将每一种样本的不同时间点的AD和Con合在一起比较，只考虑AD vs Con，进行ROC分析，将阳性的Pro进行箱线图展示
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

# 可视化常用包
library(tidyverse) # 包含了ggplot2等包
library(ggsci) # 期刊配色，单个离散配色模块中最多只有51种颜色
library(scales) # show_col
library(patchwork) # 拼接图
library(ggbreak) # 坐标轴截断
library(ggsignif) # 计算差异显著性并添加显著性线段
library(ggthemes) # 主题设置
library(stringr) # 文本处理包，包含在tidyverse中
library(ggrepel)
# -- ------------------------------------------------------------------------------------------------------------------
# MouseSerum
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/MouseSerum_Raw.Counts_CPM.quantiles.csv', row.names = 1) #行为基因，列为样本
data <- data %>% t() %>% as.data.frame()
data$Group <- ifelse(str_detect(rownames(data), 'AD'), 'AD', 'Con') %>% factor()
data <- relocate(data, Group)
SigPro <- read.csv('Results/Figure2/Files/MouseSerum_AUC_pvalue.Pos.csv', row.names = 1)
data.sub <- data %>% dplyr::select(Group, rownames(SigPro))
for (i in 2:ncol(data.sub)) {
  myTitle <- gsub("/", "_", colnames(data.sub)[i])
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(data.sub, aes(x = Group, fill = Group)) + theme_bw() +
    # geom_violin(aes_(y = as.name(colnames(data.sub)[i])), scale = 'width', width = 0.8, alpha = 0.8) +
    geom_boxplot(aes_(y = as.name(colnames(data.sub)[i])), width = 0.8, alpha = 0.8) +
    geom_jitter(aes_(y = as.name(colnames(data.sub)[i])), width = 0.3, alpha = 0.6, size = 2) + 
    geom_signif(aes_(y = as.name(colnames(data.sub)[i])), 
                comparisons = list(c('AD', 'Con')), 
                map_signif_level = F, #step_increase = 0.01, 
                tip_length = 0.01, test = "wilcox.test", textsize = 2) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure2/MouseSerum.', myTitle, '.wilcox.test.png')
  ggsave(p2, filename = fn2, width = 3, height = 3)
  fn2 <- str_c('Results/Figure2/MouseSerum.', myTitle, '.wilcox.test.pdf')
  ggsave(p2, filename = fn2, width = 3, height = 3)
}
# -- ------------------------------------------------------------------------------------------------------------------
# MouseUrine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/MouseUrine_Raw.Counts_CPM.quantiles.csv', row.names = 1) #行为基因，列为样本
data <- data %>% t() %>% as.data.frame()
data$Group <- ifelse(str_detect(rownames(data), 'AD'), 'AD', 'Con') %>% factor()
data <- relocate(data, Group)
SigPro <- read.csv('Results/Figure2/Files/MouseUrine_AUC_pvalue.Pos.csv', row.names = 1)
data.sub <- data %>% dplyr::select(Group, rownames(SigPro))
for (i in 2:ncol(data.sub)) {
  myTitle <- gsub("/", "_", colnames(data.sub)[i])
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(data.sub, aes(x = Group, fill = Group)) + theme_bw() +
    # geom_violin(aes_(y = as.name(colnames(data.sub)[i])), scale = 'width', width = 0.8, alpha = 0.8) +
    geom_boxplot(aes_(y = as.name(colnames(data.sub)[i])), width = 0.8, alpha = 0.8) +
    geom_jitter(aes_(y = as.name(colnames(data.sub)[i])), width = 0.3, alpha = 0.6, size = 2) + 
    geom_signif(aes_(y = as.name(colnames(data.sub)[i])), 
                comparisons = list(c('AD', 'Con')), 
                map_signif_level = F, #step_increase = 0.01, 
                tip_length = 0.01, test = "wilcox.test", textsize = 2) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure2/MouseUrine.', myTitle, '.wilcox.test.png')
  ggsave(p2, filename = fn2, width = 3, height = 3)
  fn2 <- str_c('Results/Figure2/MouseUrine.', myTitle, '.wilcox.test.pdf')
  ggsave(p2, filename = fn2, width = 3, height = 3)
}
# -- ------------------------------------------------------------------------------------------------------------------
# MouseNIS
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data <- read.csv('Results/MouseNIS_Raw.Counts_CPM.quantiles.csv', row.names = 1) #行为基因，列为样本
data <- data %>% t() %>% as.data.frame()
data$Group <- ifelse(str_detect(rownames(data), 'AD'), 'AD', 'Con') %>% factor()
data <- relocate(data, Group)
SigPro <- read.csv('Results/Figure2/Files/MouseNIS_AUC_pvalue.Pos.csv', row.names = 1)
data.sub <- data %>% dplyr::select(Group, rownames(SigPro))
for (i in 2:ncol(data.sub)) {
  myTitle <- gsub("/", "_", colnames(data.sub)[i])
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(data.sub, aes(x = Group, fill = Group)) + theme_bw() +
    # geom_violin(aes_(y = as.name(colnames(data.sub)[i])), scale = 'width', width = 0.8, alpha = 0.8) +
    geom_boxplot(aes_(y = as.name(colnames(data.sub)[i])), width = 0.8, alpha = 0.8) +
    geom_jitter(aes_(y = as.name(colnames(data.sub)[i])), width = 0.3, alpha = 0.6, size = 2) + 
    geom_signif(aes_(y = as.name(colnames(data.sub)[i])), 
                comparisons = list(c('AD', 'Con')), 
                map_signif_level = F, #step_increase = 0.01, 
                tip_length = 0.01, test = "wilcox.test", textsize = 2) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure2/MouseNIS.', myTitle, '.wilcox.test.png')
  ggsave(p2, filename = fn2, width = 3, height = 3)
  fn2 <- str_c('Results/Figure2/MouseNIS.', myTitle, '.wilcox.test.pdf')
  ggsave(p2, filename = fn2, width = 3, height = 3)
}

# -- ------------------------------------------------------------------------------------------------------------------
# 7.将每一种样本的不同时间点的AD和Con合在一起比较，只考虑AD vs Con，进行ROC分析，2-3个蛋白的联合ROC
# -- ------------------------------------------------------------------------------------------------------------------
# MouseSerum
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data.nor <- read.csv('Results/MouseSerum_Raw.Counts_CPM.quantiles.csv', row.names = 1)
roc.data <- t(data.nor) %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
roc.data <- relocate(roc.data, Group)
data.pos <- read.csv('Results/Figure2/Files/MouseSerum_AUC_pvalue.Pos.csv', row.names = 1)
data.pos$auc <- ifelse(data.pos$auc > 0.5, data.pos$auc, 1-data.pos$auc)
data.pos <- data.pos %>% arrange(desc(auc))
roc.data <- roc.data %>% dplyr::select(Group, rownames(data.pos)[1:5])
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
png('Results/Figure2/MouseSerum_Raw.Counts_CPM.quantiles_ROC_2-5.png', height = 400, width = 400)
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
pdf('Results/Figure2/MouseSerum_Raw.Counts_CPM.quantiles_ROC_2-5.pdf', height = 6, width = 6)
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
# MouseUrine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data.nor <- read.csv('Results/MouseUrine_Raw.Counts_CPM.quantiles.csv', row.names = 1)
roc.data <- t(data.nor) %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
roc.data <- relocate(roc.data, Group)
data.pos <- read.csv('Results/Figure2/Files/MouseUrine_AUC_pvalue.Pos.csv', row.names = 1)
data.pos$auc <- ifelse(data.pos$auc > 0.5, data.pos$auc, 1-data.pos$auc)
data.pos <- data.pos %>% arrange(desc(auc))
roc.data <- roc.data %>% dplyr::select(Group, rownames(data.pos)[1:5])
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
png('Results/Figure2/MouseUrine_Raw.Counts_CPM.quantiles_ROC_2-5.png', height = 400, width = 400)
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
pdf('Results/Figure2/MouseUrine_Raw.Counts_CPM.quantiles_ROC_2-5.pdf', height = 6, width = 6)
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
# MouseNIS
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data.nor <- read.csv('Results/MouseNIS_Raw.Counts_CPM.quantiles.csv', row.names = 1)
roc.data <- t(data.nor) %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
roc.data <- relocate(roc.data, Group)
data.pos <- read.csv('Results/Figure2/Files/MouseNIS_AUC_pvalue.Pos.csv', row.names = 1)
data.pos$auc <- ifelse(data.pos$auc > 0.5, data.pos$auc, 1-data.pos$auc)
data.pos <- data.pos %>% arrange(desc(auc))
roc.data <- roc.data %>% dplyr::select(Group, rownames(data.pos)[1:5])
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
png('Results/Figure2/MouseNIS_Raw.Counts_CPM.quantiles_ROC_2-5.png', height = 400, width = 400)
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
pdf('Results/Figure2/MouseNIS_Raw.Counts_CPM.quantiles_ROC_2-5.pdf', height = 6, width = 6)
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

