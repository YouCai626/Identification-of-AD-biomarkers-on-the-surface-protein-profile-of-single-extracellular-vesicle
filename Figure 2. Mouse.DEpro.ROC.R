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
sampleID <- read_csv('Mouse.PT_SampleID_SampleName.csv', show_col_types = F)
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

# Similar codes have been omitted
