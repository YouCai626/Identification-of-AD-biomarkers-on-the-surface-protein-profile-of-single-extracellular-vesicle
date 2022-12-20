library(tidyverse) # 包含了ggplot2等包
library(ggsci) # 期刊配色，单个离散配色模块中最多只有51种颜色
library(scales) # show_col
library(patchwork) # 拼接图
library(ggbreak) # 坐标轴截断
library(ggsignif) # 计算差异显著性并添加显著性线段
library(ggthemes) # 主题设置

library(readxl) # 读取EXCEL文件
library(clusterProfiler)
library(enrichplot) # 用于可视化的包
library(org.Hs.eg.db) # 小鼠用org.Mm.eg.db

# - -----------------------------------------------------------------------
# CAM panel的功能注释
# - -----------------------------------------------------------------------
rm(list=ls())
gc()
dt <- read.csv('CAM_Annotation/CAM_PBA_Panel.csv')
gene <- dt$Protein %>% unique()
gene <-  str_replace(gene, 'CD26', 'DPP4') %>% str_replace('Thy1', 'THY1') %>% str_replace('ITGA2b', 'ITGA2B') %>% 
  str_replace('TIM3', 'HAVCR2') %>% str_replace('CD20', 'MS4A1') %>% str_replace('ULBP2/5/6', 'ULBP2')
# translate gene name type
genelist <- bitr(gene, fromType = "SYMBOL", 
                 toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
head(genelist)
write_csv(as.data.frame(genelist),"CAM_Annotation/CAM_PBA_Panel_SYMBOL-ENSEMBL-ENTREZID.csv")
# - -----------------------------------------------------------------------
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
# simplify # 对GO富集分析结果进行精简
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
# GO over-representation analysis可视化
# -- ------------------------------------------------------------------------
rm(list=ls())
gc()
load('CAM_Annotation/CAM_PBA_Panel_enrichGO.RData')
show_col(pal_npg(alpha = 0.8)(10))
my_res.bp <- egoBP.sim0.4@result
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






