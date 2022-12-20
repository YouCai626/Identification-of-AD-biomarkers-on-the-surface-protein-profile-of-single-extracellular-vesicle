getwd()
setwd("D:/2.AD.PBA.CAM.Article.V1.Proj")
# 设置最大的内存
if(.Platform$OS.type == "windows") withAutoprint({
  memory.size()
  memory.size(TRUE)
  memory.limit()
})
# options(future.globals.maxSize = 20*1024^3) #设置为20G

# -- ------------------------------------------------------------------------------------------------------------------
# 可视化常用包
library(tidyverse) # 包含了ggplot2等包
library(ggsci) # 期刊配色，单个离散配色模块中最多只有51种颜色
library(scales) # show_col
library(patchwork) # 拼接图
library(ggbreak) # 坐标轴截断
library(ggsignif) # 计算差异显著性并添加显著性线段
library(ggthemes) # 主题设置
library(stringr) # 文本处理包，包含在tidyverse中

library(readxl) #读取EXCEL文件
library(clusterProfiler)
library(enrichplot) #用于可视化的包
library(org.Hs.eg.db) #小鼠用org.Mm.eg.db
# library(GOSemSim) # GO相似度分析
library(simplifyEnrichment) # GO富集结果整体可视化(有相似性展示)BiocManager::install("simplifyEnrichment")
library(magick) # 让simplifyEnrichment的图更好的栅格化 install.packages('magick')
library(ComplexHeatmap)

# -- ------------------------------------------------------------------------------------------------------------------
# Fig1.CAM panel的功能注释 + 每一种体液每一个时间点的差异表达Pro的功能注释
# -- ------------------------------------------------------------------------------------------------------------------
# CAM panel的功能注释
rm(list=ls())
gc()
dt <- read.csv('CAM_Annotation/CAM_PBA_Panel.csv')
gene <- dt$Protein %>% unique()
gene <-  str_replace(gene, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2')
#translate gene name type
genelist <- bitr(gene, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
write_csv(as.data.frame(genelist),"CAM_Annotation/CAM_PBA_Panel_SYMBOL-ENSEMBL-ENTREZID.csv")
# -- ------------------------------------------------------------------------------------------------------------------
# GO classification 先对gene进行GO分类
# - -----------------------------------------------------------------------
# Level 1 provides the highest list coverage with the least amount of term specificity. 
# With each increasing level coverage decreases while specificity increases so that 
# level 5 provides the least amount of coverage with the highest term specificity.
ggoBP <- groupGO(gene     = unique(genelist$ENTREZID),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "BP",
                 level    = 5,
                 readable = TRUE)
ggoBP <- filter(ggoBP, Count > 0) #去除没有意义的条目
ggoBP <- arrange(ggoBP, desc(Count)) #逆序排列
write_csv(as.data.frame(ggoBP), "CAM_Annotation/CAM_PBA_Panel_groupGO_5_BP.csv")
ggoCC <- groupGO(gene     = unique(genelist$ENTREZID),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "CC",
                 level    = 5,
                 readable = TRUE)
ggoCC <- filter(ggoCC, Count > 0)
ggoCC <- arrange(ggoCC, desc(Count)) #逆序排列
write_csv(as.data.frame(ggoCC), "CAM_Annotation/CAM_PBA_Panel_groupGO_5_CC.csv")
ggoMF <- groupGO(gene     = unique(genelist$ENTREZID),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "MF",
                 level    = 5,
                 readable = TRUE)
ggoMF <- filter(ggoMF, Count > 0)
ggoMF <- arrange(ggoMF, desc(Count)) #逆序排列
write_csv(as.data.frame(ggoMF), "CAM_Annotation/CAM_PBA_Panel_groupGO_5_MF.csv")
# 保存数据，后续作图用
save(genelist, ggoBP, ggoCC, ggoMF, 
     file = "CAM_Annotation/CAM_PBA_Panel_GO Classification5.RData")

# GO Classification4
ggoBP <- groupGO(gene     = unique(genelist$ENTREZID),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "BP",
                 level    = 4,
                 readable = TRUE)
ggoBP <- filter(ggoBP, Count > 0) #去除没有意义的条目
ggoBP <- arrange(ggoBP, desc(Count)) #逆序排列
write_csv(as.data.frame(ggoBP), "CAM_Annotation/CAM_PBA_Panel_groupGO_4_BP.csv")
ggoCC <- groupGO(gene     = unique(genelist$ENTREZID),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "CC",
                 level    = 4,
                 readable = TRUE)
ggoCC <- filter(ggoCC, Count > 0)
ggoCC <- arrange(ggoCC, desc(Count)) #逆序排列
write_csv(as.data.frame(ggoCC), "CAM_Annotation/CAM_PBA_Panel_groupGO_4_CC.csv")
ggoMF <- groupGO(gene     = unique(genelist$ENTREZID),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "MF",
                 level    = 4,
                 readable = TRUE)
ggoMF <- filter(ggoMF, Count > 0)
ggoMF <- arrange(ggoMF, desc(Count)) #逆序排列
write_csv(as.data.frame(ggoMF), "CAM_Annotation/CAM_PBA_Panel_groupGO_4_MF.csv")
# 保存数据，后续作图用
save(genelist, ggoBP, ggoCC, ggoMF, 
     file = "CAM_Annotation/CAM_PBA_Panel_GO Classification4.RData")

# GO Classification3
ggoBP <- groupGO(gene     = unique(genelist$ENTREZID),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "BP",
                 level    = 3,
                 readable = TRUE)
ggoBP <- filter(ggoBP, Count > 0) #去除没有意义的条目
ggoBP <- arrange(ggoBP, desc(Count)) #逆序排列
write_csv(as.data.frame(ggoBP), "CAM_Annotation/CAM_PBA_Panel_groupGO_3_BP.csv")
ggoCC <- groupGO(gene     = unique(genelist$ENTREZID),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "CC",
                 level    = 3,
                 readable = TRUE)
ggoCC <- filter(ggoCC, Count > 0)
ggoCC <- arrange(ggoCC, desc(Count)) #逆序排列
write_csv(as.data.frame(ggoCC), "CAM_Annotation/CAM_PBA_Panel_groupGO_3_CC.csv")
ggoMF <- groupGO(gene     = unique(genelist$ENTREZID),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "MF",
                 level    = 3,
                 readable = TRUE)
ggoMF <- filter(ggoMF, Count > 0)
ggoMF <- arrange(ggoMF, desc(Count)) #逆序排列
write_csv(as.data.frame(ggoMF), "CAM_Annotation/CAM_PBA_Panel_groupGO_3_MF.csv")
# 保存数据，后续作图用
save(genelist, ggoBP, ggoCC, ggoMF, 
     file = "CAM_Annotation/CAM_PBA_Panel_GO Classification3.RData")

# GO Classification2
ggoBP <- groupGO(gene     = unique(genelist$ENTREZID),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "BP",
                 level    = 2,
                 readable = TRUE)
ggoBP <- filter(ggoBP, Count > 0) #去除没有意义的条目
ggoBP <- arrange(ggoBP, desc(Count)) #逆序排列
write_csv(as.data.frame(ggoBP), "CAM_Annotation/CAM_PBA_Panel_groupGO_2_BP.csv")
ggoCC <- groupGO(gene     = unique(genelist$ENTREZID),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "CC",
                 level    = 2,
                 readable = TRUE)
ggoCC <- filter(ggoCC, Count > 0)
ggoCC <- arrange(ggoCC, desc(Count)) #逆序排列
write_csv(as.data.frame(ggoCC), "CAM_Annotation/CAM_PBA_Panel_groupGO_2_CC.csv")
ggoMF <- groupGO(gene     = unique(genelist$ENTREZID),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "MF",
                 level    = 2,
                 readable = TRUE)
ggoMF <- filter(ggoMF, Count > 0)
ggoMF <- arrange(ggoMF, desc(Count)) #逆序排列
write_csv(as.data.frame(ggoMF), "CAM_Annotation/CAM_PBA_Panel_groupGO_2_MF.csv")
# 保存数据，后续作图用
save(genelist, ggoBP, ggoCC, ggoMF, 
     file = "CAM_Annotation/CAM_PBA_Panel_GO Classification2.RData")

# -- ------------------------------------------------------------------------
# GO over-representation analysis
# - -----------------------------------------------------------------------
egoBP <- enrichGO(gene          = unique(genelist$ENTREZID),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  readable      = TRUE)
egoBP <- clusterProfiler::filter(egoBP, p.adjust < 0.05)
head(egoBP@result)
dim(egoBP@result)
write.csv(egoBP@result, 'CAM_Annotation/CAM_PBA_Panel_enrichGO-BP.Sig.csv', row.names = F)
# filter GO enriched result at specific level
# GO分层是有向无环图：有建议说取单纯的某一个层会有问题：有的通路，级层为3都还很冗余，而其他的通路却在3的时候已经是最细节的分层
egoBP5 <- gofilter(egoBP, level = 5)
egoBP4 <- gofilter(egoBP, level = 4)
egoBP3 <- gofilter(egoBP, level = 3)
egoBP2 <- gofilter(egoBP, level = 2)
egoBP5@result$Description %>% head(n=10)
egoBP4@result$Description %>% head(n=10)
egoBP3@result$Description %>% head(n=10)
egoBP2@result$Description %>% head(n=10)
write.csv(egoBP3@result, 'CAM_Annotation/CAM_PBA_Panel_enrichGO-BP.Sig level3.csv', row.names = F)
# simplify 对GO富集分析结果进行精简
egoBP.sim0.8 <- clusterProfiler::simplify(egoBP, cutoff = 0.8)
egoBP.sim0.7 <- clusterProfiler::simplify(egoBP, cutoff = 0.7)
egoBP.sim0.6 <- clusterProfiler::simplify(egoBP, cutoff = 0.6)
egoBP.sim0.5 <- clusterProfiler::simplify(egoBP, cutoff = 0.5)
egoBP.sim0.4 <- clusterProfiler::simplify(egoBP, cutoff = 0.4)
egoBP.sim0.3 <- clusterProfiler::simplify(egoBP, cutoff = 0.3)
egoBP.sim0.8@result$Description %>% head(n=10)
egoBP.sim0.7@result$Description %>% head(n=10)
egoBP.sim0.6@result$Description %>% head(n=10)
egoBP.sim0.5@result$Description %>% head(n=10)
egoBP.sim0.4@result$Description %>% head(n=10)
egoBP.sim0.3@result$Description %>% head(n=10)
write.csv(egoBP.sim0.4@result, 'CAM_Annotation/CAM_PBA_Panel_enrichGO-BP.Sig.sim0.4.csv', row.names = F)

egoCC <- enrichGO(gene          = unique(genelist$ENTREZID),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  readable      = TRUE)
egoCC <- clusterProfiler::filter(egoCC, p.adjust < 0.05)
head(egoCC@result)
dim(egoCC@result)
write.csv(egoCC@result, 'CAM_Annotation/CAM_PBA_Panel_enrichGO-CC.Sig.csv', row.names = F)
egoCC5 <- gofilter(egoCC, level = 5)
egoCC4 <- gofilter(egoCC, level = 4)
egoCC3 <- gofilter(egoCC, level = 3)
egoCC5@result$Description %>% head(n=10)
egoCC4@result$Description %>% head(n=10)
egoCC3@result$Description %>% head(n=10)
write.csv(egoCC3@result, 'CAM_Annotation/CAM_PBA_Panel_enrichGO-CC.Sig level3.csv', row.names = F)
egoCC.sim0.8 <- clusterProfiler::simplify(egoCC, cutoff = 0.8)
egoCC.sim0.7 <- clusterProfiler::simplify(egoCC, cutoff = 0.7)
egoCC.sim0.6 <- clusterProfiler::simplify(egoCC, cutoff = 0.6)
egoCC.sim0.5 <- clusterProfiler::simplify(egoCC, cutoff = 0.5)
egoCC.sim0.8@result$Description %>% head(n=10)
egoCC.sim0.7@result$Description %>% head(n=10)
egoCC.sim0.6@result$Description %>% head(n=10)
egoCC.sim0.5@result$Description %>% head(n=10)
write.csv(egoCC.sim0.5@result, 'CAM_Annotation/CAM_PBA_Panel_enrichGO-CC.Sig.sim0.5.csv', row.names = F)

egoMF <- enrichGO(gene          = unique(genelist$ENTREZID),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  readable      = TRUE)
egoMF <- clusterProfiler::filter(egoMF, p.adjust < 0.05)
egoMF@result$Description %>% head(n=10)
dim(egoMF@result)
write.csv(egoMF@result, 'CAM_Annotation/CAM_PBA_Panel_enrichGO-MF.Sig.csv', row.names = F)
egoMF5 <- gofilter(egoMF, level = 5)
egoMF4 <- gofilter(egoMF, level = 4)
egoMF3 <- gofilter(egoMF, level = 3)
egoMF5@result$Description %>% head(n=10)
egoMF4@result$Description %>% head(n=10)
egoMF3@result$Description %>% head(n=10)
write.csv(egoMF5@result, 'CAM_Annotation/CAM_PBA_Panel_enrichGO-MF.Sig level5.csv', row.names = F)
egoMF.sim0.8 <- clusterProfiler::simplify(egoMF, cutoff = 0.8)
egoMF.sim0.7 <- clusterProfiler::simplify(egoMF, cutoff = 0.7)
egoMF.sim0.6 <- clusterProfiler::simplify(egoMF, cutoff = 0.6)
egoMF.sim0.5 <- clusterProfiler::simplify(egoMF, cutoff = 0.5)
egoMF.sim0.4 <- clusterProfiler::simplify(egoMF, cutoff = 0.4)
egoMF.sim0.8@result$Description %>% head(n=10)
egoMF.sim0.7@result$Description %>% head(n=10)
egoMF.sim0.6@result$Description %>% head(n=10)
egoMF.sim0.5@result$Description %>% head(n=10)
egoMF.sim0.4@result$Description %>% head(n=10)
write.csv(egoMF.sim0.6@result, 'CAM_Annotation/CAM_PBA_Panel_enrichGO-MF.Sig.sim0.5.csv', row.names = F)
# 保存数据，后续作图用
save(genelist, egoBP, egoCC, egoMF, 
     egoBP.sim0.4, egoCC.sim0.5, egoMF.sim0.5,
     egoBP3, egoCC3, egoMF5,
     file = "CAM_Annotation/CAM_PBA_Panel_enrichGO.RData")

# -- ------------------------------------------------------------------------
# KEGG pathway over-representation analysis
# -- ------------------------------------------------------------------------
options(clusterProfiler.download.method = "wininet")
kk <- enrichKEGG(gene         = unique(genelist$ENTREZID),
                 organism     = 'hsa',
                 pvalueCutoff = 1,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
kk <- clusterProfiler::filter(kk, pvalue < 0.05)
head(kk@result)
write.csv(kk@result, 'CAM_Annotation/CAM_PBA_Panel_enrichKEGG.csv', row.names = F)
# 保存数据，后续作图用
save(genelist, kk,
     file = "CAM_Annotation/CAM_PBA_Panel_enrichKEGG.RData")

# -- ------------------------------------------------------------------------
# Reactome pathway over-representation analysis
# -- ------------------------------------------------------------------------
# BiocManager::install('ReactomePA', force = TRUE)
# BiocManager::install('reactome.db') # 总报错
# install.packages("D:/Downloads/RGtk2_2.20.36.3.tar.gz", repos = NULL, type = "source") # 本地安装
library(ReactomePA)
rp <- enrichPathway(gene         = unique(genelist$ENTREZID),
                    organism     = 'human',
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)
rp <- setReadable(rp, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
rp <- clusterProfiler::filter(rp, pvalue < 0.05)
head(rp@result)
write.csv(rp@result, 'CAM_Annotation/CAM_PBA_Panel_enrichReactomePathway.csv')
# 保存数据，后续作图用
save(genelist, rp,
     file = "CAM_Annotation/CAM_PBA_Panel_enrichReactomePathway.RData")

# -- ------------------------------------------------------------------------
# Disease over-representation analysis
# -- ------------------------------------------------------------------------
library(DOSE)
do <- enrichDO(gene          = unique(genelist$ENTREZID),
               ont           = "DO",
               pvalueCutoff  = 1,
               pAdjustMethod = "BH",
               minGSSize     = 5,
               qvalueCutoff  = 1,
               readable      = T)
do <- clusterProfiler::filter(do, pvalue < 0.05)
head(do@result)
write.csv(do@result, 'CAM_Annotation/CAM_PBA_Panel_enrichDO.csv', row.names = F)
# 保存数据，后续作图用
save(genelist, do,
     file = "CAM_Annotation/CAM_PBA_Panel_enrichDO.RData")

# -- ------------------------------------------------------------------------
# WikiPathways over-representation analysis
# -- ------------------------------------------------------------------------
wp <- enrichWP(gene          = unique(genelist$ENTREZID),
               organism      = "Homo sapiens",
               pvalueCutoff  = 1,
               qvalueCutoff  = 1
               ) 
wp <- setReadable(wp, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
wp <- clusterProfiler::filter(wp, pvalue < 0.05)
head(wp@result)
write.csv(wp@result, 'CAM_Annotation/CAM_PBA_Panel_enrichWikiPathways.csv', row.names = F)
# 保存数据，后续作图用
save(genelist, wp,
     file = "CAM_Annotation/CAM_PBA_Panel_enrichWP.RData")



# -- ------------------------------------------------------------------------
# GO Classification可视化
# -- ------------------------------------------------------------------------
rm(list=ls())
gc()
load('CAM_Annotation/CAM_PBA_Panel_GO Classification3.RData')
show_col(pal_aaas(alpha = 0.7)(10))
p <- barplot(ggoBP, showCategory = 20, label_format = 60) + ggtitle("GO-BP Classification3 Top20")
ggsave('CAM_Annotation/CAM_PBA_Panel_GO-BP Classification3 Top20.barplot.png', width = 7, height = 6)
p <- barplot(ggoCC, showCategory = 20, label_format = 60) + ggtitle("GO-CC Classification3 Top20")
ggsave('CAM_Annotation/CAM_PBA_Panel_GO-CC Classification3 Top20.barplot.png', width = 7, height = 6)
p <- barplot(ggoMF, showCategory = 20, label_format = 60) + ggtitle("GO-MF Classification3 Top20")
ggsave('CAM_Annotation/CAM_PBA_Panel_GO-MF Classification3 Top20.barplot.png', width = 7, height = 6)
# 自己用ggplot2画图
ggoBP_res <- ggoBP@result
p <- ggoBP_res %>% arrange(Count) %>% filter(Count > 60) %>% 
  mutate(Description = factor(Description, levels = Description)) %>%
  ggplot(aes(x = Count , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#008280B2') + 
  geom_vline(xintercept = 50,
             color = "white",
             linetype = 2, 
             size = 0.6) + 
  xlab('') + ylab('') + 
  theme(axis.text.x = element_text(color = 'black', size = rel(1.5)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'CAM_Annotation/CAM_PBA_Panel_GO-BP Classification3 Count60.barplot.png', width = 6, height = 6)

# 可进一步筛选想关注的通路，比如：神经/炎症/免疫/细胞死亡等于AD相关的通路
c1 <- str_subset(ggoBP@result$Description, '(N|n)eu')
c2 <- str_subset(ggoBP@result$Description, 'synapse')
c3 <- str_subset(ggoBP@result$Description, '(I|i)nflamm')
c4 <- str_subset(ggoBP@result$Description, '(I|i)mmu')
c5 <- str_subset(ggoBP@result$Description, '(D|d)eath')
c6 <- str_subset(ggoBP@result$Description, 'osis')
c7 <- str_subset(ggoBP@result$Description, 'commu')
c8 <- str_subset(ggoBP@result$Description, 'signal')
c9 <- str_subset(ggoBP@result$Description, 'localiza')
c10 <- str_subset(ggoBP@result$Description, 'trans')
c11 <- str_subset(ggoBP@result$Description, 'acti')
c12 <- str_subset(ggoBP@result$Description, 'prolife')
c13 <- str_subset(ggoBP@result$Description, 'cycle')
c14 <- str_subset(ggoBP@result$Description, 'vesi')
c <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14)
ggoBP1 <- filter(ggoBP, Description %in% c)
write_csv(as.data.frame(ggoBP1), "CAM_Annotation/CAM_PBA_Panel_groupGO3_BP-AD相关.csv")
ggoBP1_res <- ggoBP1@result
p <- ggoBP1_res %>% arrange(Count) %>% filter(Count > 38) %>% 
  mutate(Description = factor(Description, levels = Description)) %>%
  ggplot(aes(x = Count , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#00828099') + 
  geom_vline(xintercept = 30,
             color = "white",
             linetype = 2, 
             size = 0.6) + 
  xlab('') + ylab('') + 
  # ggtitle('GO-BP Classification') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.4)),
        axis.text.y = element_text(color = 'black', size = rel(1.4)))
ggsave(p, filename = 'CAM_Annotation/CAM_PBA_Panel_GO-BP Classification3 AD相关.barplot.png', width = 9, height = 3.5)

# Tang说用人为手法挑选
ggoBP1 <- filter(ggoBP, Count >= 3)
c <- ggoBP1@result$Description
c1 <- str_subset(c, 'regulation')
c <- setdiff(c, c1)
c2 <- str_subset(c, 'response')
c2 <- c2[-c(1,3:5,10:12)]
c <- setdiff(c, c2)
c3 <- str_subset(c, 'develop')
c <- setdiff(c, c3)
c4 <- str_subset(c, 'multi')
c <- setdiff(c, c4)
c5 <- str_subset(c, 'based')
c <- setdiff(c, c5)
c6 <- str_subset(c, 'cycle')
c <- setdiff(c, c6[-2])
c7 <- str_subset(c, 'process')
c <- setdiff(c, c7[-c(2,7,14)])
c8 <- str_subset(c, 'anato')
c <- setdiff(c, c8)
ggoBP1 <- filter(ggoBP1, Description %in% c)
write_csv(as.data.frame(ggoBP1), "CAM_Annotation/CAM_PBA_Panel_groupGO_3_BP.Less.csv")

# 自己用ggplot2画图
ggoBP1_res <- ggoBP1@result
p <- ggoBP1_res %>% arrange(Count) %>% 
  mutate(Description = factor(Description, levels = Description)) %>%
  ggplot(aes(x = Count , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#00828099') + 
  geom_vline(xintercept = 50,
             color = "white",
             linetype = 2, 
             size = 0.6) + 
  xlab('') + ylab('') + ggtitle('GO_Biological Process') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.5)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'CAM_Annotation/CAM_PBA_Panel_GO-BP Classification3.Less.barplot.png', width = 7, height = 9)

ggoBP1_res <- read.csv('CAM_Annotation/CAM_PBA_Panel_groupGO_3_BP.Less.Tang.csv')
p <- ggoBP1_res %>% arrange(Count) %>% 
  mutate(Description = factor(Description, levels = Description)) %>%
  ggplot(aes(x = Count , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#00828099') + 
  geom_vline(xintercept = 10,
             color = "white",
             linetype = 2, 
             size = 0.6) + 
  xlab('') + ylab('') + ggtitle('GO_Biological Process') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.5)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'CAM_Annotation/CAM_PBA_Panel_GO-BP Classification3.Less.Tang.barplot.png', width = 7, height = 8)

# -- ------------------------------------------------------------------------
# GO over-representation analysis可视化
# -- ------------------------------------------------------------------------
rm(list=ls())
gc()
load('CAM_Annotation/CAM_PBA_Panel_enrichGO.RData')
# -- ------------------------------------------------------------------------
# 用simplifyEnrichment全局展示富集结果
# -- ------------------------------------------------------------------------
go_id = egoBP@result$ID
mat = GO_similarity(go_id, db = 'org.Hs.eg.db', ont = 'BP')
png('CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP.Sig simplifyEnrichment.compare.png',
    width = 10, height = 6, units = 'in', bg = "white", res = 300)
set.seed(2022)
compare_clustering_methods(mat, verbose = F)
dev.off()
png('CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP.Sig simplifyEnrichment.compare.heatmap.png',
    width = 16, height = 12, units = 'in', bg = "white", res = 300)
set.seed(2022)
compare_clustering_methods(mat, plot_type = "heatmap", verbose = F)
dev.off()
pdf('CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP.Sig simplifyEnrichment.compare.pdf', width = 10, height = 6)
set.seed(2022)
compare_clustering_methods(mat, verbose = F)
dev.off()
pdf('CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP.Sig simplifyEnrichment.compare.heatmap.pdf', width = 16, height = 12)
set.seed(2022)
compare_clustering_methods(mat, plot_type = "heatmap", verbose = F)
dev.off()

all_clustering_methods()
set.seed(2022)
cl = binary_cut(mat)
png('CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP simplifyEnrichment.binary_cut.png',
    width = 10, height = 5, units = 'in', bg = "white", res = 300)
ht_clusters(mat, cl, word_cloud_grob_param = list(max_width = 100),
            exclude_words = c('regulation', 'positive', 'negative','response', 'process','cell','pathway','signaling'),
            order_by_size = TRUE,
            max_words = 9,
            fontsize_range = c(8, 16)
)
dev.off()
pdf('CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP simplifyEnrichment.binary_cut.pdf', width = 10, height = 5)
ht_clusters(mat, cl, word_cloud_grob_param = list(max_width = 100),
            exclude_words = c('regulation', 'positive', 'negative','response', 'process','cell','pathway','signaling'),
            order_by_size = TRUE,
            max_words = 9,
            fontsize_range = c(8, 16)
)
dev.off()
?ht_clusters
cl = cluster_by_kmeans(mat) # 除了binary_cut直接用之外，其他的是cluster_by_*

all_clustering_methods()
png('CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP simplifyEnrichment.binary_cut.simplifyGO.png',
    width = 10, height = 5, units = 'in', bg = "white", res = 300)
set.seed(2022)
df = simplifyGO(mat, method = "binary_cut", word_cloud_grob_param = list(max_width = 100),
                exclude_words = c('regulation', 'positive', 'negative','response', 'process','cell','pathway','signaling'),
                order_by_size = TRUE,
                max_words = 9,
                fontsize_range = c(8, 16),
                column_title = NULL,
                verbose = F)
dev.off()
head(df)
write.csv(df, file = 'CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP simplifyEnrichment.binary_cut.csv', row.names = F)
png('CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP simplifyEnrichment.walktrap.png',
    width = 8, height = 4, units = 'in', bg = "white", res = 300)
set.seed(2022)
df = simplifyGO(mat, method = "walktrap", word_cloud_grob_param = list(max_width = 80),
                exclude_words = c('regulation', 'positive', 'negative','response', 'process','cell','pathway','signaling'),
                order_by_size = TRUE,
                max_words = 8,
                fontsize_range = c(8, 16),
                column_title = NULL,
                verbose = F)
dev.off()
head(df)
write.csv(df, file = 'CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP simplifyEnrichment.walktrap.csv', row.names = F)
pdf('CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP simplifyEnrichment.walktrap.pdf', width = 8, height = 4)
set.seed(2022)
df = simplifyGO(mat, method = "walktrap", word_cloud_grob_param = list(max_width = 80),
                exclude_words = c('regulation', 'positive', 'negative','response', 'process','cell','pathway','signaling'),
                order_by_size = TRUE,
                max_words = 8,
                fontsize_range = c(8, 16),
                column_title = NULL,
                verbose = F)
dev.off()
png('CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP simplifyEnrichment.louvain.png',
    width = 8, height = 4, units = 'in', bg = "white", res = 300)
set.seed(2022)
df = simplifyGO(mat, method = "louvain", word_cloud_grob_param = list(max_width = 80),
                exclude_words = c('regulation', 'positive', 'negative','response', 'process','cell','pathway','signaling'),
                order_by_size = TRUE,
                max_words = 9,
                fontsize_range = c(8, 16),
                column_title = NULL,
                verbose = F)
dev.off()
head(df)
write.csv(df, file = 'CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP simplifyEnrichment.louvain.csv', row.names = F)
pdf('CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP simplifyEnrichment.louvain.pdf', width = 8, height = 4)
set.seed(2022)
df = simplifyGO(mat, method = "louvain", word_cloud_grob_param = list(max_width = 80),
                exclude_words = c('regulation', 'positive', 'negative','response', 'process','cell','pathway','signaling'),
                order_by_size = TRUE,
                max_words = 9,
                fontsize_range = c(8, 16),
                column_title = NULL,
                verbose = F)
dev.off()
# -- ------------------------------------------------------------------------
# 自己用ggplot2画图
show_col(pal_npg(alpha = 0.8)(10))
my_res.bp <- egoBP.sim0.4@result
# p <- my_res.bp %>% arrange(pvalue) %>% filter(pvalue <= my_res.bp$pvalue[10]) %>% arrange(Count) %>% 
#   mutate(Description = factor(Description, levels = Description)) %>%
#   ggplot(aes(x = Count , y = Description)) + theme_bw() +
#   geom_bar(stat = "identity", 
#            width = 0.8, 
#            position = position_dodge(0.8),
#            fill = '#F39B7FFF') + 
#   geom_vline(xintercept = 10,
#              color = "white",
#              linetype = 2, 
#              size = 0.6) + 
#   xlab('') + ylab('') + 
#   # ggtitle('GO-BP Classification') +
#   theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
#         axis.text.y = element_text(color = 'black', size = rel(1.2)))
p <- my_res.bp %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.bp$p.adjust[10]) %>% 
  mutate(logP = -log(p.adjust)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#F39B7FFF') +  
  xlab('') + ylab('') + 
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP.Sig.sim0.4 Top10.barplot.png', width = 5, height = 2.5)
ggsave(p, filename = 'CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP.Sig.sim0.4 Top10.barplot.pdf', width = 5, height = 2.5)
my_res.cc <- egoCC.sim0.5@result
p <- my_res.mf %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.mf$p.adjust[10]) %>% 
  mutate(logP = -log(p.adjust)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#4DBBD5FF') +  
  xlab('') + ylab('') + 
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-CC.Sig.sim0.5 Top10.barplot.png', width = 3.7, height = 2.5)
ggsave(p, filename = 'CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-CC.Sig.sim0.5 Top10.barplot.pdf', width = 3.7, height = 2.5)
my_res.mf <- egoMF5@result
p <- my_res.mf %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.mf$p.adjust[10]) %>% 
  mutate(logP = -log(p.adjust)) %>% arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  ggplot(aes(x = logP , y = Description)) + theme_bw() +
  geom_bar(stat = "identity", 
           width = 0.8, 
           position = position_dodge(0.8),
           fill = '#8491B4FF') + 
  xlab('') + ylab('') + 
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-MF.Sig level5 Top10.barplot.png', width = 3.13, height = 2.5)
ggsave(p, filename = 'CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-MF.Sig level5 Top10.barplot.pdf', width = 3.13, height = 2.5)

# -- ------------------------------------------------------------------------
# 自带函数展示结果
p <- barplot(egoBP, showCategory = 10, label_format = 40) + ggtitle("BP") + xlab('') + 
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichGO-BP Top10.barplot.png', width = 4, height = 4)
p <- barplot(egoBP.sim0.7, showCategory = 10, label_format = 40) + ggtitle("BP") + xlab('') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichGO-BP.sim0.7 Top10.barplot.png', width = 4, height = 4)
p <- barplot(egoBP.sim0.6, showCategory = 10, label_format = 40) + ggtitle("BP") + xlab('') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichGO-BP.sim0.6 Top10.barplot.png', width = 4, height = 4)
p <- barplot(egoBP.sim0.5, showCategory = 10, label_format = 40) + ggtitle("BP") + xlab('') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichGO-BP.sim0.5 Top10.barplot.png', width = 4, height = 4)
p <- barplot(egoBP.sim0.4, showCategory = 10, label_format = 40) + ggtitle("BP") + xlab('') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichGO-BP.sim0.4 Top10.barplot.png', width = 4, height = 4)
p <- barplot(egoCC, showCategory = 10, label_format = 40) + ggtitle("CC") + xlab('') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichGO-CC Top10.barplot.png', width = 4, height = 4)
p <- barplot(egoCC.sim0.7, showCategory = 10, label_format = 40) + ggtitle("CC") + xlab('') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichGO-CC.sim0.7 Top10.barplot.png', width = 4, height = 4)
p <- barplot(egoCC.sim0.6, showCategory = 10, label_format = 40) + ggtitle("CC") + xlab('') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichGO-CC.sim0.6 Top10.barplot.png', width = 4, height = 4)
p <- barplot(egoMF, showCategory = 10, label_format = 40) + ggtitle("MF") + xlab('') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichGO-MF Top10.barplot.png', width = 4, height = 4)
p <- barplot(egoMF.sim0.7, showCategory = 10, label_format = 40) + ggtitle("MF") + xlab('') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichGO-MF.sim0.7 Top10.barplot.png', width = 4, height = 4)
p <- barplot(egoMF.sim0.6, showCategory = 10, label_format = 40) + ggtitle("MF") + xlab('') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichGO-MF.sim0.6 Top10.barplot.png', width = 4, height = 4)

# -- ------------------------------------------------------------------------
# 按照自己的思路, 筛选AD相关的通路
c1 <- str_subset(egoBP.sim0.7@result$Description, 'amyloid')
c2 <- str_subset(egoBP.sim0.7@result$Description, '(T|t)au')
c3 <- str_subset(egoBP.sim0.7@result$Description, '(N|n)euroinflamma')
c4 <- str_subset(egoBP.sim0.7@result$Description, 'synap')
c5 <- str_subset(egoBP.sim0.7@result$Description, 'transmi')
c6 <- str_subset(egoBP.sim0.7@result$Description, 'vesi')
c7 <- str_subset(egoBP.sim0.7@result$Description, '(D|d)eath')
c8 <- str_subset(egoBP.sim0.7@result$Description, 'glia')
c9 <- str_subset(egoBP.sim0.7@result$Description, '(A|a)stro')
c10 <- str_subset(egoBP.sim0.7@result$Description, 'osis')
c11 <- str_subset(egoBP.sim0.7@result$Description, 'totic')
c12 <- str_subset(egoBP.sim0.7@result$Description, '(I|i)mmu')
c <- c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) %>% unique()
c13 <- str_subset(c, 'regulation')
c <- setdiff(c, c13)
egoBP.sim0.7.sub <- filter(egoBP.sim0.7, Description %in% c)
write_csv(as.data.frame(egoBP.sim0.7.sub), "CAM_Annotation/CAM_PBA_Panel_enrichGO-BP.sim0.7_AD相关.csv")
p <- barplot(egoBP.sim0.7.sub, color = 'pvalue', showCategory = 20, label_format = 50) + ggtitle("GO-BP") + xlab('') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        # legend.position = 'none'
        )
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichGO-BP.sim0.7_AD相关.barplot.png', width = 6, height = 8)
p <- dotplot(egoBP.sim0.7.sub, color = 'pvalue', showCategory = 20, label_format = 50) + ggtitle("GO-BP") + xlab('') +
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        # legend.position = 'none'
  )
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichGO-BP.sim0.7_AD相关.dotplot.png', width = 6, height = 8)

# -- ------------------------------------------------------------------------
# KEGG over-representation analysis可视化
# -- ------------------------------------------------------------------------
rm(list=ls())
gc()
load('CAM_Annotation/CAM_PBA_Panel_enrichKEGG.RData')
p <- barplot(kk, showCategory = 10, label_format = 40) + ggtitle("KEGG") + xlab('') + 
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichKEGG Top10.barplot.png', width = 4, height = 4)

# -- ------------------------------------------------------------------------
# Disease over-representation analysis可视化
# -- ------------------------------------------------------------------------
rm(list=ls())
gc()
load('CAM_Annotation/CAM_PBA_Panel_enrichDO.RData')
p <- barplot(do, showCategory = 10, label_format = 40) + ggtitle("enrichDisease") + xlab('') + 
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichDO Top10.barplot.png', width = 4, height = 4)

# -- ------------------------------------------------------------------------
# ReactomePathway over-representation analysis可视化
# -- ------------------------------------------------------------------------
rm(list=ls())
gc()
load('CAM_Annotation/CAM_PBA_Panel_enrichReactomePathway.RData')
p <- barplot(rp, showCategory = 10, label_format = 40) + ggtitle("ReactomePathway") + xlab('') + 
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichReactomePathway Top10.barplot.png', width = 4, height = 4)

# -- ------------------------------------------------------------------------
# WikiPathway over-representation analysis可视化
# -- ------------------------------------------------------------------------
rm(list=ls())
gc()
load('CAM_Annotation/CAM_PBA_Panel_enrichWP.RData')
p <- barplot(wp, showCategory = 10, label_format = 40) + ggtitle("WikiPathway") + xlab('') + 
  theme(axis.text.x = element_text(color = 'black', size = rel(1.2)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)),
        legend.position = 'none')
ggsave('CAM_Annotation/CAM_PBA_Panel_enrichWikiPathways Top10.barplot.png', width = 4, height = 4)
