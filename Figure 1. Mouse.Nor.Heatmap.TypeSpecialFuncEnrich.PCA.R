# Set the maximum memory
if(.Platform$OS.type == "windows") withAutoprint({
  memory.size()
  memory.size(TRUE)
  memory.limit()
})

# -- ------------------------------------------------------------------------------------------------------------------
# 1.Normalization of Mouse's expression values
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse) # Processing data
library(edgeR) # cpm standardization
library(preprocessCore) # normalize.quantiles function for quantile normalization

rm(list=ls())
gc()
data <- read.csv('Results/Mouse_Raw.Counts.csv',row.names = 1)
data$mean <- rowMeans(data)
data <- data %>% filter(mean > 0) %>% dplyr::select(!mean) # Remove the Pro with expression 0
data.nor <- data %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(data.nor) <- rownames(data)
names(data.nor) <- names(data)
write.csv(data.nor, 'Results/Mouse_Raw.Counts_CPM.quantiles.csv')

# -- ------------------------------------------------------------------------------------------------------------------
# 2.draw heatmap
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse) 
library(ggsci) 
library(scales)
library(pheatmap)

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
dt <- dt %>% filter(mean > 0) %>% dplyr::select(!mean) # Remove the Pro with expression 0
dt.nor <- dt %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(dt.nor) <- rownames(dt)
names(dt.nor) <- names(dt)
write.csv(dt.nor, 'Results/Mouse.15M_Raw.Counts_CPM.quantiles.csv')
dt.nor$var <- apply(dt.nor, 1, var) # Find the variance for each row
range(dt.nor$var)
dt.nor <- dt.nor %>% filter(var > 0) %>% dplyr::select(!var) # Remove the Pro with all equal expressions
ann_colors.sub = list(Group = c('AD' ='#CC0000FF', 'Con' ='#6699FFFF'),
                      Type = c('Serum' = '#D595A7CC', 'Urine' = '#F0E685CC', 'NIS' = '#6BD76BCC')) 
pheatmap(dt.nor, scale = 'row', color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 2.5,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = F,
         fontsize_number = 4, fontsize_row = 2.5, fontsize_col = 10,
         cluster_rows = T, cluster_cols = T, treeheight_row = 30, treeheight_col = 20,
         annotation_col = coldata.sub, annotation_colors = ann_colors.sub,
         filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_pheatmap.png')
pheatmap(dt.nor, scale = 'row', color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 2.5,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = F,
         fontsize_number = 4, fontsize_row = 2.5, fontsize_col = 10,
         cluster_rows = T, cluster_cols = T, treeheight_row = 30, treeheight_col = 20,
         annotation_col = coldata.sub, annotation_colors = ann_colors.sub,
         filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_pheatmap.pdf')
# Get the expression matrix after scale
dt.nor_scaled <- t(scale(t(dt.nor))) %>% as.data.frame()
head(dt.nor_scaled)
write.csv(dt.nor_scaled, 'Results/Figure1/Files/Mouse.15M_Raw.Counts_CPM.quantiles_scaleRow.csv')

# -- ------------------------------------------------------------------------------------------------------------------
# 3.Define Type Special according to the criterion that there are more than 5 expressions with values greater than 1 in the group after scale
#   and less than 6 expressions with values greater than 1 outside the group.
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse) 

# 15M
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
ann_colors = list(Group = c('AD' ='#CC0000FF', 'Con' ='#6699FFFF'),
                  Type = c('Serum' = '#D595A7CC', 'Urine' = '#F0E685CC', 'NIS' = '#6BD76BCC'))
pheatmap(data.nor, color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 6,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = T,
         fontsize_number = 4, fontsize_row = 6, fontsize_col = 8, border_color = NA,
         cluster_rows = T, cluster_cols = F, treeheight_row = 40,
         annotation_col = coldata, annotation_colors = ann_colors,
         cutree_rows = 3,
         gaps_col = c(8,16),
         filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_pheatmap_TypeSpecial.png')
pheatmap(data.nor, color = pal_gsea(palette = c("default"), n = 60)(60), cellwidth = 10, cellheight = 6,
         display_numbers = F, number_color = 'black', number_format = '%.3f', show_rownames = T,
         fontsize_number = 4, fontsize_row = 6, fontsize_col = 8, border_color = NA,
         cluster_rows = T, cluster_cols = F, treeheight_row = 40,
         annotation_col = coldata, annotation_colors = ann_colors,
         cutree_rows = 3,
         gaps_col = c(8,16),
         filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_pheatmap_TypeSpecial.pdf')
# -- ------------------------------------------------------------------------------------------------------------------
# 4.Tissue-specific proteins for functional analysis
# -- ------------------------------------------------------------------------------------------------------------------
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggsci)
library(scales)

# Urine
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
# GO over-representation analysis
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
egoBP_Urine.sim0.6 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.6)
egoBP_Urine.sim0.5 <- clusterProfiler::simplify(egoBP_Urine, cutoff = 0.5)
egoBP_Urine.sim0.6@result$Description %>% head(n=10)
egoBP_Urine.sim0.5@result$Description %>% head(n=10)
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
egoCC_Urine.sim0.6 <- clusterProfiler::simplify(egoCC_Urine, cutoff = 0.6)
egoCC_Urine.sim0.5 <- clusterProfiler::simplify(egoCC_Urine, cutoff = 0.5)
egoCC_Urine.sim0.6@result$Description %>% head(n=10)
egoCC_Urine.sim0.5@result$Description %>% head(n=10)
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
egoMF_Urine.sim0.6 <- clusterProfiler::simplify(egoMF_Urine, cutoff = 0.6)
egoMF_Urine.sim0.5 <- clusterProfiler::simplify(egoMF_Urine, cutoff = 0.5)
egoMF_Urine.sim0.6@result$Description %>% head(n=10)
egoMF_Urine.sim0.5@result$Description %>% head(n=10)
# Save the data
save(data, myGene, 
     egoBP_Urine, egoCC_Urine, egoMF_Urine, 
     egoBP_Urine.sim0.6, egoCC_Urine.sim0.6, egoMF_Urine.sim0.6,
     egoBP_Urine.sim0.5, egoCC_Urine.sim0.5, egoMF_Urine.sim0.5,
     file = "Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO.RData")
# Visualization of GO over-representation analysis
rm(list=ls())
gc()
load('Results/Figure1/Files/15M_Raw.Counts_CPM.quantiles_UrineSpecial_enrichGO.RData')
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
# 5.PCA analysis
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(FactoMineR) # PCA analysis
library(factoextra) # Visualization of PCA analysis
library(ggthemes)

# 15M
rm(list=ls())
gc()
data <- read.csv('Results/Mouse.15M_Raw.Counts_CPM.quantiles.csv', row.names = 1)
data$var <- apply(data, 1, var)
range(data$var)
data <- data %>% filter(var > 0) %>% dplyr::select(!var)
data.pca <- data %>% t() %>% as.data.frame() # When PCA analyzes data, rows are samples and columns are variables (dimensions, components).
data.pca$Group <- ifelse(str_detect(rownames(data.pca),'D'), 'AD', 'Con') %>% factor(levels = c('AD','Con'))
data.pca$Type <- ifelse(str_detect(rownames(data.pca),'s'), 'Serum',
                        ifelse(str_detect(rownames(data.pca),'u'), 'Urine', 'NIS')) %>% 
  factor(levels = c('Serum', 'Urine', 'NIS'))
data.pca$Color <- str_c(data.pca$Type, data.pca$Group, sep = '.') %>% 
  factor(levels = c('Serum.AD','Serum.Con','Urine.AD','Urine.Con','NIS.AD','NIS.Con'))
data.pca <- data.pca %>% relocate(Group, Type, Color)
pca <- PCA(data.pca[,-c(1:3)], graph = F)
p <- fviz_pca_ind(pca,
                  mean.point = F, 
                  geom.ind = "point", 
                  repel = T, 
                  pointsize = 3, 
                  pointshape = 21,
                  fill.ind = data.pca$Color, 
                  alpha.ind = 0.8,
                  palette = c('#AE1F63CC','#D595A7CC','#CC9900CC','#F0E685CC','#6BD76BCC','#749B58CC'),
                  title = '15M') + theme_bw() + theme(legend.title=element_blank())
ggsave(p, filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_PCA.png', width = 4.2, height = 3)
ggsave(p, filename = 'Results/Figure1/Mouse.15M_Raw.Counts_CPM.quantiles_PCA.pdf', width = 4.2, height = 3)

# The analysis for the other groups was similar
