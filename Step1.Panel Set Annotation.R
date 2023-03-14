library(tidyverse) # includes packages such as ggplot2
library(ggsci) # journal color schemes, up to 51 colors in a single discrete palette
library(scales) # show_col
library(patchwork) # for creating composite plots
library(ggbreak) # for breaking the axis
library(ggsignif) # for calculating and adding significance brackets
library(ggthemes) # for theme settings
library(readxl) # for reading Excel files
library(clusterProfiler)
library(enrichplot) # for visualization
library(org.Hs.eg.db) # for mouse, use org.Mm.eg.db

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
# Streamline GO enrichment analysis results with simplify - remove redundancy
egoBP.sim0.6 <- clusterProfiler::simplify(egoBP, cutoff = 0.6)
egoBP.sim0.5 <- clusterProfiler::simplify(egoBP, cutoff = 0.5)
egoBP.sim0.6@result$Description %>% head(n=10)
egoBP.sim0.5@result$Description %>% head(n=10)

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
egoCC.sim0.6 <- clusterProfiler::simplify(egoCC, cutoff = 0.6)
egoCC.sim0.5 <- clusterProfiler::simplify(egoCC, cutoff = 0.5)
egoCC.sim0.6@result$Description %>% head(n=10)
egoCC.sim0.5@result$Description %>% head(n=10)

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
egoMF.sim0.6 <- clusterProfiler::simplify(egoMF, cutoff = 0.6)
egoMF.sim0.5 <- clusterProfiler::simplify(egoMF, cutoff = 0.5)
egoMF.sim0.6@result$Description %>% head(n=10)
egoMF.sim0.5@result$Description %>% head(n=10)
# Save the data
save(genelist, egoBP, egoCC, egoMF, 
     egoBP.sim0.5, egoCC.sim0.5, egoMF.sim0.5,
     egoBP.sim0.6, egoCC.sim0.6, egoMF.sim0.6,
     file = "CAM_Annotation/CAM_PBA_Panel_enrichGO.RData")

# -- ------------------------------------------------------------------------
# Visualization of GO over-representation analysis
# -- ------------------------------------------------------------------------
rm(list=ls())
gc()
load('CAM_Annotation/CAM_PBA_Panel_enrichGO.RData')
show_col(pal_npg(alpha = 0.8)(10))
my_res.bp <- egoBP.sim0.5@result
p <- my_res.bp %>% arrange(p.adjust) %>% filter(p.adjust <= my_res.bp$p.adjust[15]) %>% 
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
  xlab('-Log10(p.adjust)') + ylab('') + ggtitle('GO-BP') +
  theme(axis.text.x = element_text(color = 'black', size = rel(0.9)),
        axis.text.y = element_text(color = 'black', size = rel(1.2)))
ggsave(p, filename = 'CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP.Sig.sim0.5 Top15.barplot.png', width = 6, height = 3)
ggsave(p, filename = 'CAM_Annotation/Plots/CAM_PBA_Panel_enrichGO-BP.Sig.sim0.5 Top15.barplot.pdf', width = 6, height = 3)

# Visualization of CC and MF groups is the same as BP

