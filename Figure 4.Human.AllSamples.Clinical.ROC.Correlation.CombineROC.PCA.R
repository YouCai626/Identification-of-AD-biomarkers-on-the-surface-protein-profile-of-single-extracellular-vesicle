# Set the maximum memory
if(.Platform$OS.type == "windows") withAutoprint({
  memory.size()
  memory.size(TRUE)
  memory.limit()
})

# Commonly used packages for visualization
library(tidyverse) 
library(ggsci) 
library(scales) 
library(ggthemes) 
library(patchwork) 
library(ggsignif) 
library(ggpubr) 
library(ggrepel)

# -- ------------------------------------------------------------------------------------------------------------------
# 1.Visualization of clinical information from Human samples：Age；Gender；MMSE；MoCA-B
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
  theme(strip.text.x = element_text(size = 14))
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
# 2.Get the raw counts for each sample
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID$EVcounts <- c()
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('RawData/',sampleID[i,1],'.total_ev_protein.csv')
  data <- read.csv(fn1, header = T)
  wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
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
write_csv(sampleID, 'Human_SampleID_SampleName EV counts.csv')

# -- ------------------------------------------------------------------------------------------------------------------
# 3.Constructed as an expression matrix
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
sampleID <- read.csv('Human_SampleID_SampleName.csv')
fn <- str_c('pseudoData/Raw/', sampleID[1,1], '.Raw.csv')
Raw.Counts <- read_csv(fn, show_col_types = F)
for (i in 2:nrow(sampleID)) {
  fn <- str_c('pseudoData/Raw/', sampleID[i,1], '.Raw.csv')
  data <- read_csv(fn,show_col_types = F)
  Raw.Counts <- full_join(Raw.Counts, data, by = 'Pro', fill = 0)
}
Raw.Counts <- relocate(Raw.Counts, Pro)
Raw.Counts[is.na(Raw.Counts)] <- 0  # Replace NA with 0
write_csv(Raw.Counts, 'Results/Human_Raw.Counts.csv')
# Separate samples of each body fluid
rm(list=ls())
gc()
data <- read.csv('Results/Human_Raw.Counts.csv', row.names = 1)
sampleID <- read_csv('Human_SampleID_SampleName.csv', show_col_types = F)
# Urine
sub.ID <- sampleID %>% filter(Source == 'Urine')
sub.data <- data %>% dplyr::select(sub.ID$SampleName)
sub.data$mean <- rowMeans(sub.data)
sub.data <- sub.data %>% filter(mean > 0) %>% dplyr::select(!mean) # Remove the Pro with expression 0
write_csv(sub.ID, 'Results/HumanUrine_Samples.Annotation.csv')
write.csv(sub.data, 'Results/HumanUrine_Raw.Counts.csv')

# Other groups is similar to the urine group

# -- ------------------------------------------------------------------------------------------------------------------
# 4.ROC analysis + joint ROC for multiple Pro
# -- ------------------------------------------------------------------------------------------------------------------
library(edgeR) 
library(preprocessCore) 
library(pROC) 
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
  if (pvalue[j] < 0.05 & auc[j] >= 0.7) {
    proteinname <- gsub("/", "_", colnames(roc.data)[j + 1])
    fn6 <- paste('Results/Figure4/HumanUrine_Raw.Counts_cpm.quantiles_', proteinname, '.png', sep = '')
    png(fn6, width = 400, height = 400)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure4/HumanUrine_Raw.Counts_cpm.quantiles_', proteinname, '.pdf', sep = '')
    pdf(fn6, width = 5, height = 5)
    plot.roc(rocPlot, col="red", lwd = 4, print.thres = TRUE, print.thres.cex = 1.2, legacy.axes = TRUE,
             print.auc =TRUE, print.auc.col = "blue", print.auc.cex = 1.2, auc.polygon = TRUE, auc.polygon.col = "pink",
             print.auc.x=0.7, print.auc.y=0.1,
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
             print.auc.x=0.7, print.auc.y=0.1,
             main = colnames(roc.data)[j + 1],
             cex.axis=1.4, cex.lab=2, cex.main=2)
    dev.off()
    fn6 <- paste('Results/Figure4/HumanUrine_Raw.Counts_cpm.quantiles_', proteinname, '.pdf', sep = '')
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
auc_pvalue <- arrange(auc_pvalue, desc(auc))
write.csv(auc_pvalue, file = 'Results/Figure4/Files/HumanUrine_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.csv')
auc_pvalue <- filter(auc_pvalue, pvalue < 0.05)
write.csv(auc_pvalue, file = 'Results/Figure4/Files/HumanUrine_Raw.Counts_CPM.quantiles_ROC_AUC.pvalue.Pos.csv')
# Correlation of ROC-positive Protein with MMSE and MoCA-B
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
names(dt) <- names(dt) %>% str_replace_all('/', '.') # An error will be reported if the file name has '/' in it
for (i in 1:length(myPro)) {
  p <- ggplot(dt, aes_(x = as.name(colnames(dt)[i+4]))) + theme_bw() +
    geom_point(aes(y = MMSE, color = Group), size = 3, alpha = 0.8) + 
    scale_color_manual(values = c('#CC0000FF','#6699FFFF')) +
    stat_smooth(aes(y = MMSE), method='lm', formula = y~x, colour='red', linetype = 2, size = 1, alpha = 0.2) +
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
# Graphing the expression of ROC-positive Protein
for (i in 5:ncol(dt)) {
  myTitle <- gsub("/", "_", colnames(dt)[i])
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
# Joint ROC for 2-4 Pro
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
plot.roc(roc.pro5, col = "#DC0000CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#DC0000CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.2,
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()
pdf('Results/Figure4/HumanUrine_Raw.Counts_CPM.quantiles_ROC_2-5.pdf', height = 8, width = 8)
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
plot.roc(roc.pro5, col = "#DC0000CC", lwd = 3, add = T, print.thres = F, print.thres.cex = 1.2, legacy.axes = TRUE,
         print.auc =TRUE, print.auc.col = "#DC0000CC", print.auc.cex = 1.2,
         print.auc.x=0.3, print.auc.y=0.2,
         cex.axis=1.4, cex.lab=2, cex.main=2)
dev.off()

# -- ------------------------------------------------------------------------------------------------------------------
# 5.PCA analysis
# -- ------------------------------------------------------------------------------------------------------------------
library(FactoMineR)
library(factoextra)
# -- ------------------------------------------------------------------------------------------------------------------
# HumanUrine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data.Urine <- read.csv('Results/HumanUrine_Raw.Counts_cpm.quantiles.csv', row.names = 1)
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
                  mean.point = F, 
                  # label = "ind", # Display the label for each sample
                  geom.ind = "point",
                  repel = T, # Avoid overlap
                  pointsize = 3,      
                  pointshape = 21,   
                  fill.ind = data.pca$Group, 
                  alpha.ind = 0.8,
                  palette = c('#ED0000CC','#00468BCC'), # Set n colors for n groups
                  title = 'HumanUrine PCA') + theme_bw() + theme(legend.position = 'none')
ggsave(p, filename = 'Results/Figure4/HumanUrine_ROC.Pos_PCA.png', width = 4, height = 4.2)
ggsave(p, filename = 'Results/Figure4/HumanUrine_ROC.Pos_PCA.pdf', width = 4, height = 4.2)
# -- ------------------------------------------------------------------------------------------------------------------
# 4 body fluids combined
# -- ------------------------------------------------------------------------------------------------------------------
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
                  mean.point = F,
                  geom.ind = "point", 
                  repel = T, 
                  pointsize = 3,  
                  pointshape = 21,
                  fill.ind = data.pca$color,
                  alpha.ind = 0.8,
                  palette = c('#AE1F63CC','#D595A7CC','#925E9FFF','#ADB6B6CC','#CC9900CC','#F0E685CC','#6BD76BCC','#749B58CC'), # Set n colors for n groups
                  title = 'Human PCA') + theme_bw() + theme(legend.title=element_blank())
ggsave(p, filename = 'Results/Figure4/Human_PCA.png', width = 5, height = 4)
ggsave(p, filename = 'Results/Figure4/Human_PCA.pdf', width = 5, height = 4)


# The analysis for other groups is similar to the urine group
