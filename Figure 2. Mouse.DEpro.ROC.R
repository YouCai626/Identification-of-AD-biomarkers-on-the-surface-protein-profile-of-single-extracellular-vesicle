# Set the maximum memory
if(.Platform$OS.type == "windows") withAutoprint({
  memory.size()
  memory.size(TRUE)
  memory.limit()
})

# -- ------------------------------------------------------------------------------------------------------------------
# 1. Wring out Raw.Counts for each sample at each time point separately, because to analyze separately
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(edgeR)
library(preprocessCore)

rm(list=ls())
gc()
data <- read.csv('Results/Mouse_Raw.Counts.csv', row.names = 1)
sampleID <- read_csv('Mouse_SampleID_SampleName.csv', show_col_types = F)
# MouseUrine
sub.ID <- sampleID %>% filter(Source == 'Urine')
sub.data <- data %>% dplyr::select(sub.ID$SampleName)
sub.data$mean <- rowMeans(sub.data)
sub.data <- sub.data %>% filter(mean > 0) %>% dplyr::select(!mean)
write_csv(sub.ID, 'Results/MouseUrine_Samples.Annotation.csv')
write.csv(sub.data, 'Results/MouseUrine_Raw.Counts.csv')
sub.data.nor <- sub.data %>% edgeR::cpm() %>% normalize.quantiles() %>% round(., 3) %>% as.data.frame()
rownames(sub.data.nor) <- rownames(sub.data)
names(sub.data.nor) <- names(sub.data)
write.csv(sub.data.nor, 'Results/MouseUrine_Raw.Counts_CPM.quantiles.csv')

# -- ------------------------------------------------------------------------------------------------------------------
# 2. Differentially expressed proteins at each time point for each body fluid
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse) 
library(limma)

# Urine
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

# -- ------------------------------------------------------------------------------------------------------------------
# 3. Combining the positive DEpro
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
rm(list=ls())
gc()
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

# -- ------------------------------------------------------------------------------------------------------------------
# 4. Overlap
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data.Serum <- read.csv('Results/Figure2/Files/MouseSerum.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
data.Serum <- data.Serum %>% filter(AveExpr >= 9)
data.Urine <- read.csv('Results/Figure2/Files/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
data.Urine <- data.Urine %>% filter(AveExpr >= 9)
data.NIS <- read.csv('Results/Figure2/Files/MouseNIS.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
data.NIS <- data.NIS %>% filter(AveExpr >= 9)

library(tidyverse) 
library(ggsci) 
library(scales) 
library(ggthemes) 
library(ggsignif) 
library(ggrepel)
library(VennDiagram)

venn_list <- list(Serum = data.Serum$X, Urine = data.Urine$X,
                  NIS = data.NIS$X)
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

# The functional analysis and its visualization are similar to the previous sections

# -- ------------------------------------------------------------------------------------------------------------------
# 5. Visualization of differentially expressed proteins
# -- ------------------------------------------------------------------------------------------------------------------
# Draw a scatter chart with trend lines
# Urine
rm(list = ls())
gc()
select = dplyr::select
data <- read.csv('Results/Figure2/Files/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma.csv')
data <- data %>% filter(AveExpr >= 9)
data$Month <- factor(data$Month, levels = c('3M','10M','15M'))
data$Time <- ifelse(data$Month == '3M', 3, ifelse(data$Month == '10M', 10, 15))
data$Sig <- ifelse(data$logFC > 0, 'Up', 'Down')
data$Sig <- data$Sig %>% factor(levels = c('Up','Down'))
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
ggsave('Results/Figure2/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma_Up.png', width = 4, height = 3)
ggsave('Results/Figure2/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma_Up.pdf', width = 4, height = 3)
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
ggsave('Results/Figure2/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma_Down.png', width = 4, height = 3)
ggsave('Results/Figure2/MouseUrine.3-10-15M.ADvsCon_DEpro.Pos_Limma_Down.pdf', width = 4, height = 3)

# Functional enrichment analysis of the differentially expressed proteins of 3-10-15M of Urine respectively

# -- ------------------------------------------------------------------------------------------------------------------
# 6. AD and Con at different time points for each sample were compared together, considering only AD vs Con, for ROC analysis
# -- ------------------------------------------------------------------------------------------------------------------
library(pROC)

# MouseUrine
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
             print.auc.x=0.7, print.auc.y=0.1,
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)  
    dev.off()
    filename6 <- paste('Results/Figure2/MouseUrine.', proteinname, '.pdf', sep = '')
    pdf(filename6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,
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
             print.auc.x=0.7, print.auc.y=0.1,
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure2/MouseUrine.', proteinname, '.pdf', sep = '')
    pdf(fn6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "yellow",
             print.auc.x=0.7, print.auc.y=0.1,
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
# 7. Box plot display of positive Pro
# -- ------------------------------------------------------------------------------------------------------------------
library(ggsignif)
# MouseUrine
rm(list=ls())
gc()
data <- read.csv('Results/MouseUrine_Raw.Counts_CPM.quantiles.csv', row.names = 1)
data <- data %>% t() %>% as.data.frame()
data$Group <- ifelse(str_detect(rownames(data), 'AD'), 'AD', 'Con') %>% factor()
data <- relocate(data, Group)
SigPro <- read.csv('Results/Figure2/Files/MouseUrine_AUC_pvalue.Pos.csv', row.names = 1)
data.sub <- data %>% dplyr::select(Group, rownames(SigPro))
for (i in 2:ncol(data.sub)) {
  myTitle <- gsub("/", "_", colnames(data.sub)[i])
  # you need to pass the result of the expression with the aes_() function
  p2 <- ggplot(data.sub, aes(x = Group, fill = Group)) + theme_bw() +
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
# 8.Combined ROC of 2-4 proteins
# -- ------------------------------------------------------------------------------------------------------------------
library(glmnet)
library(pROC)

# MouseUrine
rm(list=ls())
gc()
data.nor <- read.csv('Results/MouseUrine_Raw.Counts_CPM.quantiles.csv', row.names = 1)
roc.data <- t(data.nor) %>% as.data.frame()
roc.data$Group <- ifelse(str_detect(rownames(roc.data), 'D'), 'AD', 'Con') %>% factor(levels = c('AD', 'Con'))
roc.data <- relocate(roc.data, Group)
data.pos <- read.csv('Results/Figure2/Files/MouseUrine_AUC_pvalue.Pos.csv', row.names = 1)
data.pos$auc <- ifelse(data.pos$auc > 0.5, data.pos$auc, 1-data.pos$auc)
data.pos <- data.pos %>% arrange(desc(auc))
roc.data <- roc.data %>% dplyr::select(Group, rownames(data.pos)[1:4])
f2 <- glm(Group ~ roc.data[,2] + roc.data[,3], data = roc.data, family = binomial())
roc.data$Pro2 <- predict(f2, newdata = roc.data, type = "response")
f3 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4], data = roc.data, family = binomial())
roc.data$Pro3 <- predict(f3, newdata = roc.data, type = "response")
f4 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4] + roc.data[,5], data = roc.data, family = binomial())
roc.data$Pro4 <- predict(f4, newdata = roc.data, type = "response")
roc.pro2 <- roc(roc.data$Group, roc.data$Pro2)
roc.pro3 <- roc(roc.data$Group, roc.data$Pro3)
roc.pro4 <- roc(roc.data$Group, roc.data$Pro4)
show_col(pal_npg(alpha = 0.8)(10))
png('Results/Figure2/MouseUrine_Raw.Counts_CPM.quantiles_ROC_2-4.png', height = 400, width = 400)
plot.roc(roc.pro2, col = "#3C5488CC", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#3C5488CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.35,
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro3, col = "#F39B7FCC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F39B7FCC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.3,
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro4, col = "#4DBBD5CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#4DBBD5CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.25,
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
pdf('Results/Figure2/MouseUrine_Raw.Counts_CPM.quantiles_ROC_2-4.pdf', height = 6, width = 6)
plot.roc(roc.pro2, col = "#3C5488CC", lwd = 3, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#3C5488CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.35,
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro3, col = "#F39B7FCC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#F39B7FCC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.3,
         cex.axis=1.4, cex.lab=2, cex.main=2)
plot.roc(roc.pro4, col = "#4DBBD5CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#4DBBD5CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.25,
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()

# -- ------------------------------------------------------------------------------------------------------------------
# 9.Pseudo processing of samples
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
# 10.Correlation of pseudo samples with raw
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
# 11.Machine learning
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
# Normalize the data and fill in the missing values, this can be done with the preProcess function
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

# -- ------------------------------------------------------------------------------------------------------------------
# 12.Correlation of protein used for modeling with mouse Aβ concentration
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

# Similar codes have been omitted
