# Set the maximum memory
if(.Platform$OS.type == "windows") withAutoprint({
  memory.size()
  memory.size(TRUE)
  memory.limit()
})

library(tidyverse) 
library(Seurat) 
library(harmony)

# -- ------------------------------------------------------------------------------------------------------------------
# 1.Because the number of EVs per sample in the original data is very large and sparse, 
# n EVs are synthesized into 1 super EV to perform Seurat analysis: Per10
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
sampleID <- read.csv('Mouse_SampleID_SampleName EV counts.csv')
sampleID <- sampleID %>% filter(Source == 'Urine')
pseudo.size = 10 # Every 10 EVs are combined into 1 EV
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('RawData/Mouse/',sampleID[i,1],'.total_ev_protein.csv')
  data <- read.csv(fn1, header = T)
  wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
  rownames(wide_data) <- wide_data$ev
  wide_data <- wide_data[, -c(1:2)]
  superEv10 <- floor(nrow(wide_data)/pseudo.size) 
  if (superEv10 > 0) {
    set.seed(2022)
    random.row <- sample(1:nrow(wide_data), superEv10*pseudo.size, replace = F)
    raw.data <- wide_data[random.row,] %>% as.data.frame()
    raw.data$Group <- rep(c(1:superEv10), each = pseudo.size)
    raw.data <- relocate(raw.data, Group)
    pseudo.data <- aggregate(list(raw.data[2:ncol(raw.data)]), by = list(raw.data$Group), FUN = sum)
    rownames(pseudo.data) <- str_c(sampleID[i,2], 'sEV', 1:superEv10, sep = '_')
    fn2 <- str_c('pseudoData/pseudoPer10/Mouse/',sampleID[i,1],'.sEv10.csv')
    write.csv(pseudo.data[-1], fn2)
    rm(data, wide_data,raw.data,pseudo.data)
  }
}

# -- ------------------------------------------------------------------------------------------------------------------
# 2.Create Seurat Objects
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
# Urine
sampleID <- read.csv('Mouse_SampleID_SampleName EV counts.csv')
sampleID <- sampleID %>% filter(Source == 'Urine')
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('pseudoData/pseudoPer10/Mouse/',sampleID[i,1],'.sEv10.csv')
  data <- read.csv(fn1, header = T, row.names = 1)
  data <- data %>% t() %>% as.data.frame()
  if (ncol(data) >= 100) {
    seuOb <- CreateSeuratObject(counts = data, project = sampleID[i,2], min.features = 2, min.cells = 30)
    # 给SeuratObjects添加metadata
    seuOb$Species <- sampleID[i,3]
    seuOb$Source <- sampleID[i,4]
    seuOb$Diseases <- sampleID[i,5]
    seuOb$Gender <- sampleID[i,6]
    seuOb$Age <- sampleID[i,7]
    seuOb$Condition <- sampleID[i,8]
    seuOb$Batch <- sampleID[i,9]
    fn2 <- str_c('SeuratObjects/pseudoPer10/Mouse/',sampleID[i,1],'.pseudoPer10.rds')
    write_rds(seuOb, fn2)
    rm(data, seuOb)
  }
}

# -- ------------------------------------------------------------------------------------------------------------------
# 3.SeuratObjects merge
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
# Urine
sampleID <- read.csv('Mouse_SampleID_SampleName EV counts.csv')
sampleID <- sampleID %>% filter(Source == 'Urine')
filenames <- str_c('SeuratObjects/pseudoPer10/Mouse/', sampleID$SampleID,'.pseudoPer10.rds')
seuList <- lapply(filenames, readRDS)
seuOb <- merge(seuList[[1]], seuList[-1], project = 'pseudoPer10')
head(seuOb@meta.data)
table(seuOb@meta.data$orig.ident)
sampleID <- str_split(rownames(seuOb@meta.data), '_s') %>% as.data.frame() %>% t() %>% as.data.frame()
seuOb@meta.data$orig.ident <- sampleID$V1
write_rds(seuOb, 'SeuratObjects/MouseUrine.pseudoPer10.Merge.rds')

# -- ------------------------------------------------------------------------------------------------------------------
# 4.SCTransform
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
seuOb <- readRDS('SeuratObjects/MouseUrine.pseudoPer10.Merge.rds')
seuOb <- seuOb %>% SCTransform(return.only.var.genes = F, verbose = F)
write_rds(seuOb, 'SeuratObjects/MouseUrine.pseudoPer10.Merge_SCT.rds')

# -- ------------------------------------------------------------------------------------------------------------------
# 5.Harmony
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
dir.create('Results/MouseUrine/Seurat')
seuOb <- readRDS('SeuratObjects/MouseUrine.pseudoPer10.Merge_SCT.rds')
seuOb <- RunPCA(seuOb, assay = "SCT", verbose = F) %>% 
  RunHarmony(assay.use = "SCT", reduction = "pca", group.by.vars = 'orig.ident', plot_convergence = T)
ggsave('Results/MouseUrine/Seurat/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Convergence.png', width = 6, height = 4)
write_rds(seuOb, 'SeuratObjects/MouseUrine.pseudoPer10.Merge_SCT_Harmony.rds')
table(seuOb@meta.data$orig.ident)

# -- ------------------------------------------------------------------------------------------------------------------
# 6.UMAP
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
seuOb <- readRDS('SeuratObjects/MouseUrine.pseudoPer10.Merge_SCT_Harmony.rds')
seuOb <- RunUMAP(seuOb, assay = "SCT", reduction = "harmony", dims = 1:20, verbose = F) %>% 
  RunTSNE(assay = "SCT", reduction = "harmony", dims = 1:20, check_duplicates = F) %>% 
  FindNeighbors(assay = "SCT", reduction = "harmony", dims = 1:20, verbose = F)
for (a in c(0.1, 0.3)){
  seuOb <- FindClusters(seuOb, resolution = a, verbose = F)
}
write_rds(seuOb, 'SeuratObjects/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_Res0.3.rds')
library(qs)
qsave(seuOb, file = 'SeuratObjects/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_Res0.3.qs')

# -- ------------------------------------------------------------------------------------------------------------------
# 7.only a portion of the EV is extracted for the Dimensionality reduction
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
seuOb <- qread('SeuratObjects/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_Res0.3.qs')
Idents(seuOb) <- 'orig.ident'
set.seed(2023)
sub.seuOb <- subset(seuOb, downsample = 1000)
Idents(sub.seuOb) <- 'SCT_snn_res.0.1'
p1 <- DimPlot(sub.seuOb, reduction = "umap", label = TRUE, raster = F, repel = T) + NoLegend()
p2 <- DimPlot(sub.seuOb, reduction = "tsne", label = TRUE, raster = F, repel = T)
p3 <- p1 + p2
ggsave(p3, filename = 'Results/MouseUrine/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.1.1000.png', 
       height = 5, width = 11)
ggsave(p3, filename = 'Results/MouseUrine/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.1.1000.pdf', 
       height = 5, width = 11)
Idents(seuOb) <- 'SCT_snn_res.0.1'
seuOb.markers <- FindAllMarkers(seuOb, assay = 'RNA', min.pct = 0.2, verbose = F)
top5 <- seuOb.markers %>% group_by(cluster) %>% top_n(n = 5, wt = -p_val_adj)
p4 <- DoHeatmap(sub.seuOb, features = top5$gene) + NoLegend()
ggsave(p4, file = 'Results/MouseUrine/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_0.1.1000.clusters.markers.Top5.DoHeatmap.png',
       height = 7, width = 20)
p5 <- DotPlot(sub.seuOb, features = unique(top5$gene), cols = c("lightgrey", "red")) + RotatedAxis()
ggsave(p5, file = 'Results/MouseUrine/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_0.1.1000.clusters.markers.Top5.DotPlot.png',
       height = 8, width = 8)
ggsave(p5, file = 'Results/MouseUrine/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_0.1.1000.clusters.markers.Top5.DotPlot.pdf',
       height = 8, width = 8)

# -- ------------------------------------------------------------------------------------------------------------------
# 8.average expression of each subcluster
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
dir.create('Results/MouseUrine/Seurat/Expression')
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
seuOb <- qread('SeuratObjects/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_Res0.3.qs')
seuOb@meta.data$orig.ident %>% table()
DefaultAssay(seuOb) 
Idents(seuOb) <- 'SCT_snn_res.0.1' 
table(seuOb@active.ident)
for (i in 0:20) {
  gc()
  sub_seuOb <- subset(seuOb, seurat_clusters == i)
  avrExp <- AverageExpression(sub_seuOb,group.by = 'orig.ident')
  fn <- str_c('Results/MouseUrine/Seurat/Expression/pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster',i,'.Expression.SCT.csv')
  write.csv(avrExp[['SCT']], fn)
}

# -- ------------------------------------------------------------------------------------------------------------------
# 9.Perform ROC analysis on Train's data and do AUC trend plots for the joint ROC of 1-5 proteins
# -- ------------------------------------------------------------------------------------------------------------------
library(tidyverse) 
library(glmnet)
library(caret)
library(gbm)
library(pROC)
dir.create('Results/MouseUrine/Logistic')
dir.create('Results/MouseUrine/caret')
dir.create('Results/MouseUrine/caret/RDS')
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
# a.Filtering variables with Logidtic regression
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
for( j in 0:20) {
  fn1 <- str_c('Results/MouseUrine/Seurat/Expression/pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster',j,'.Expression.SCT.csv')
  data <- read.csv(fn1, row.names = 1)
  data <- data %>% t() %>% as.data.frame()
  data$ID <- rownames(data) %>% factor
  data$Group <- ifelse(str_detect(rownames(data), 'D'),'AD','Con') %>% factor(levels = c('AD','Con'))
  data <- relocate(data, ID, Group)
  # a.Dividing the training set and the test set:
  data.ID <- levels(data$ID) 
  data.ID <- data.frame(ID = data.ID)
  data.ID$Group <- ifelse(str_detect(data.ID$ID, 'D'), 'AD', 'Con')
  data.ID$Class <- ifelse(str_detect(data.ID$ID, '15'), 15,
                          ifelse(str_detect(data.ID$ID, '10'), 10, 3))
  data.ID$Class <- factor(data.ID$Class, levels = c(3,10,15), labels = c('3M', '10M', '15M'))
  data.ID$Split <- str_c(data.ID$Class, data.ID$Group, sep = '.')
  data.ID$Split <- factor(data.ID$Split, levels = c('3M.AD', '3M.Con', '10M.AD', '10M.Con', '15M.AD', '15M.Con'))
  set.seed(2020)
  train.ID <- createDataPartition(data.ID$Split, p = 0.5, list = FALSE)
  train.ID <- data.ID$ID[train.ID]
  data.train <- data %>% filter(ID %in% train.ID) %>% dplyr::select(!ID)
  data.test <- data %>% filter(!ID %in% train.ID) %>% dplyr::select(!ID)
  data.train <- data.train %>% select_if(is.numeric) %>% select_if(~sum(.) > 0)
  data.train$Group <- ifelse(str_detect(rownames(data.train), 'D'),'AD','Con') %>% factor(levels = c('AD','Con'))
  data.train <- relocate(data.train, Group)
  my.var <- names(data.train)
  my.var <- my.var[-1]
  res <- c()
  # b.Filtering variables with Logidtic regression
  for (i in 1:length(my.var)) {
    fit <- glm(substitute(Group ~ x, list(x = as.name(my.var[i]))), family = binomial(), data = data.train)
    coef.data <- summary(fit)$coef
    or.data <- exp(cbind('OR' = coef(fit), confint(fit)))
    res1 <- cbind(coef.data, or.data) %>% as.data.frame()
    rownames(res1)[1] <- str_c('Intercept', my.var[i], sep = '.')
    res <- rbind(res, res1[-1,]) %>% as.data.frame()
  }
  res <- res %>% arrange(`Pr(>|z|)`)
  fn <- str_c('Results/MouseUrine/Logistic/pseudoPer10.Res0.1.Cluster',j,'_Train.2020_Logistic.csv')
  write.csv(res, fn)
  myPro <- rownames(res)[1:3]
  data.train <- data.train %>% select(Group, all_of(myPro))
  data.test <- data.test %>% select(Group, all_of(myPro))
  # c.Normalize the data and fill in the missing values, this can be done with the preProcess function
  preProcValues <- preProcess(data.train, method = c('center','scale')) # 用于进行其他验证数据集的标准化
  train.transformed <- predict(preProcValues, data.train)
  test.transformed <- predict(preProcValues, data.test)
  # Save data
  fn2 <- str_c('Results/MouseUrine/caret/RDS/pseudoPer10.Res.0.1.Cluster',j,'_Train.2020_LogisticTop3_caret.RData')
  save(data.train, data.test, preProcValues, train.transformed, test.transformed,
       file = fn2)
}
# -- ------------------------------------------------------------------------------------------------------------------
# d.Model Training and Tuning
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
myMethod <- c('glmnet','kknn','C5.0','rf','AdaBag','nnet','svmLinear','svmPoly','svmRadial','nb')
fitControl = trainControl(method = 'cv', number = 5,
                          classProbs = T, summaryFunction = twoClassSummary,
                          search = 'random')
for( j in 0:20) {
  fn2 <- str_c('Results/MouseUrine/caret/RDS/pseudoPer10.Res.0.1.Cluster',j,'_Train.2020_LogisticTop3_caret.RData')
  load(file = fn2)
  for (i in 1:length(myMethod)) {
    set.seed(2022)
    Fit <- train(Group ~ ., data = train.transformed, method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
                 trControl = fitControl, 
                 tuneLength = 10, 
                 metric = 'ROC', verbose = F) 
    fn3 <- str_c('Results/MouseUrine/caret/RDS/pseudoPer20.Res0.1.Cluster',j,'_Train.2020_', myMethod[i], '.Fit.rds')
    saveRDS(Fit, fn3)
  }
}
# -- ------------------------------------------------------------------------------------------------------------------
# e.Model prediction and evaluation
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
for ( j in 0:20) {
  fn2 <- str_c('Results/MouseUrine/caret/RDS/pseudoPer10.Res.0.1.Cluster',j,'_Train.2020_LogisticTop3_caret.RData')
  load(file = fn2)
  fn4 <- str_c('Results/MouseUrine/caret/RDS/pseudoPer20.Res0.1.Cluster',j,'_Train.2020_*.Fit.rds')
  filenames <- Sys.glob(fn4)
  fn5 <- str_c('Results/MouseUrine/caret/RDS/pseudoPer20.Res0.1.Cluster',j,'_Train.2020_')
  myMethod <- filenames %>% str_remove(fn5) %>% str_remove('.Fit.rds')
  res <- data.frame()
  for ( a in 1:10) {
    Fit <- readRDS(filenames[a])
    predict.train <- predict(Fit, newdata = train.transformed)
    predict.test <- predict(Fit, newdata = test.transformed)
    # The ROC analysis should distinguish the direction, otherwise it is easy to get the opposite AUC
    roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
    roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
    Fit.roc <- data.frame(Train = roc.train, Test = roc.test) %>% mutate(method = myMethod[a])
    rownames(Fit.roc) <- myMethod[a]
    res <- rbind(Fit.roc, res)
    fn7 <- str_c('Results/MouseUrine/caret/RDS/pseudoPer10.Res0.1.Cluster',j,'_Train.2020_caret.ROC.csv')
    write.csv(res, fn7)
  }
}
filenames <- Sys.glob('Results/MouseUrine/caret/RDS/pseudoPer10.Res0.1.Cluster*_Train.2020_caret.ROC.csv')
clusters <- filenames %>% str_remove('Results/MouseUrine/caret/RDS/pseudoPer10.Res0.1.Cluster') %>%
  str_remove('_Train.2020_caret.ROC.csv')
data <- read.csv(filenames[1])
data$Cluster <- clusters[1]
for (i in 2:length(clusters)) {
  temp <- read.csv(filenames[i])
  temp$Cluster <- clusters[i]
  data <- rbind(data, temp)
}
write_csv(data, file = 'Results/MouseUrine/caret/pseudoPer10.Res0.1.subClusters_Train.2020_caret.ROC.csv')
data <- data %>% filter(Test >= 0.65)
write_csv(data, file = 'Results/MouseUrine/caret/pseudoPer10.Res0.1.subClusters_Train.2020_caret.ROC_0.65.csv')

# -- ------------------------------------------------------------------------------------------------------------------
# 10.Roc visualization for each subcluster
# -- ------------------------------------------------------------------------------------------------------------------
library(pheatmap)
rm(list=ls())
gc()
data1 <- read.csv('Results/MouseUrine/caret/pseudoPer10.Res0.1.subClusters_Train.1_caret.ROC.csv')
data1$ID <- str_c(data1$Cluster, data1$method, sep = '_')
data2 <- read.csv('Results/MouseUrine/caret/pseudoPer10.Res0.1.subClusters_Train.2_caret.ROC.csv')
data2$ID <- str_c(data2$Cluster, data2$method, sep = '_')
data3 <- read.csv('Results/MouseUrine/caret/pseudoPer10.Res0.1.subClusters_Train.3_caret.ROC.csv')
data3$ID <- str_c(data3$Cluster, data3$method, sep = '_')
data <- full_join(data1, data2, by = 'ID')
data <- full_join(data, data3, by = 'ID')
data <- data %>% dplyr::select(ID, Cluster, method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data$Train.Mean <- rowMeans(data[4:6])
data$Train.sd <- apply(data[4:6], 1, sd)
data$Test.Mean <- rowMeans(data[7:9])
data$Test.sd <- apply(data[7:9], 1, sd)
data <- data %>% arrange(desc(Test.Mean), Cluster, method)
data.Urine <- filter(data, Test.Mean >= 0.75)
rownames(data.Urine) <- data.Urine$ID
data.Urine.heatmap <- data.Urine %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)
# Plot
data.heatmap <- data.Urine.heatmap
data.heatmap <- data.heatmap %>% round(., 3)
data.heatmap$Train.Show <- str_c(data.heatmap$Train.Mean, '±', data.heatmap$Train.sd, sep = ' ')
data.heatmap$Test.Show <- str_c(data.heatmap$Test.Mean, '±', data.heatmap$Test.sd, sep = ' ')
data.heatmap[1:2] %>% range()
breakList <- seq(0.75, 1, 0.05)
pheatmap(data.heatmap[1:2], color = topo.colors(9, alpha = 0.6), 
         cellwidth = 40, cellheight = 10,
         # display_numbers = T, 
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6], 
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         # gaps_row = c(10, 29, 43),
         filename = 'Results/MouseUrine/MouseUrine.subClusters_caret.ROC.pdf') 

# -- ------------------------------------------------------------------------------------------------------------------
# 11.Protein expression of the best subcluster
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data1 <- read.csv('Results/MouseUrine/Logistic/pseudoPer10.Res0.1.Cluster20_Train.1_Logistic.csv')
data2 <- read.csv('Results/MouseUrine/Logistic/pseudoPer10.Res0.1.Cluster20_Train.2_Logistic.csv')
data3 <- read.csv('Results/MouseUrine/Logistic/pseudoPer10.Res0.1.Cluster20_Train.3_Logistic.csv')
my_features <- c(data1$X[1:3], data2$X[1:3], data3$X[1:3]) %>% unique()
dt <- read.csv('Results/MouseUrine/Seurat/Expression/pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster20.Expression.SCT.csv', row.names = 1)
dt <- dt %>% t() %>% as.data.frame()
dt.sub <- dt %>% select(all_of(my_features))
dt.sub$Group <- ifelse(str_detect(rownames(dt.sub), 'D'), 'AD', 'Con')
dt.sub <- relocate(dt.sub, Group)
for (i in 2:ncol(dt.sub)) {
  myTitle <- colnames(dt.sub)[i]
  p2 <- ggplot(dt.sub, aes(x = Group, fill = Group)) + theme_bw() +
    # geom_violin(aes_(y = as.name(colnames(dt.sub)[i])), scale = 'width', width = 0.7, alpha = 0.7) +
    geom_boxplot(aes_(y = as.name(colnames(dt.sub)[i])), width = 0.6, size = 0.6) +  
    geom_jitter(aes_(y = as.name(colnames(dt.sub)[i])), alpha = 0.6, size = 2) + 
    geom_signif(aes_(y = as.name(colnames(dt.sub)[i])), comparisons = list(c("AD","Con")), 
                map_signif_level = F, step_increase = 0.05, tip_length = 0.01, 
                test = "t.test", textsize = 2.5) +
    xlab('') + ylab('SCT Value') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/MouseUrine/MouseUrine.pseudoPer10.Dim1-20.0.1.Cluster20.',myTitle, '.Pro.png')
  ggsave(p2, filename = fn2, width = 2.5, height = 3.5)
  fn2 <- str_c('Results/MouseUrine/MouseUrine.pseudoPer10.Dim1-20.0.1.Cluster20.',myTitle, '.Pro.pdf')
  ggsave(p2, filename = fn2, width = 2.5, height = 3.5)
}

# -- ------------------------------------------------------------------------------------------------------------------
# Visualization of EV Counts
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
seuOb <- qread('SeuratObjects/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_Res0.3.qs')
table(seuOb$orig.ident)
sub.seuOb_20 <- subset(seuOb, SCT_snn_res.0.1 == 20)
table(sub.seuOb_20$orig.ident)
a <- table(sub.seuOb_20@meta.data$orig.ident) %>% as.data.frame()
dt <- a
names(dt) <- c('Samples', 'EV')
dt$Group <- ifelse(str_detect(dt$Samples, 'D'), 'AD', 'Con')
p <- ggplot(dt, aes(x = Group, y = EV, fill = Group)) + theme_bw() +
  geom_boxplot(width = 0.8, size = 0.6) +
  geom_jitter(width = 0.3, size = 2, alpha = 0.8) +
  geom_signif(comparisons = list(c("AD","Con")), 
              map_signif_level = T, step_increase = 0.05, tip_length = 0.01, 
              test = "wilcox.test", textsize = 2.5) +
  xlab('') + ylab('EV Counts') + ggtitle('Urine Cluster 20') +
  scale_fill_manual(values = c('#CC0000FF','#6699FFFF'))
ggsave(plot = p, width = 3.5, height = 4,
       filename = 'Results/MouseUrine/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_0.1_Cluster20.png')
ggsave(plot = p, width = 3.5, height = 4,
       filename = 'Results/MouseUrine/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_0.1_Cluster20.pdf')

# PLAUR+ EV Counts
sub.seuOb_PLAUR <- subset(seuOb, PLAUR > 0)
a <- table(sub.seuOb_PLAUR@meta.data$orig.ident, sub.seuOb_PLAUR@meta.data$SCT_snn_res.0.1)
write.csv(a, 'Results/MouseUrine/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_0.1_PLAUR_subCluster.csv')
dt <- a %>% as.data.frame()
names(dt) <- c('Samples', 'Cluster', 'EV')
dt <- dt %>% filter(Cluster %in% 0:20)
dt$Group <- ifelse(str_detect(dt$Samples, 'D'), 'AD', 'Con')
dt$Cluster <- factor(dt$Cluster, levels = 0:20)
dt <- dt %>% filter(Cluster %in% c(0,1,3,19,20))
p <- ggplot(dt, aes(x = Group, y = EV, fill = Group)) + theme_bw() +
  geom_boxplot(width = 0.8, size = 0.6) +
  geom_jitter(width = 0.3, size = 0.5, alpha = 0.8) +
  geom_signif(comparisons = list(c("AD","Con")), 
              map_signif_level = T, step_increase = 0.05, tip_length = 0.01, 
              test = "wilcox.test", textsize = 2.5) +
  xlab('') + ylab('EV Counts') + ggtitle('PLAUR') +
  scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
  theme(legend.position = "none") + facet_grid(~Cluster) + 
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.border = element_blank())
ggsave(plot = p, width = 4, height = 4,
       filename = 'Results/MouseUrine/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_0.1_PLAUR_Clusters_less.png')
ggsave(plot = p, width = 4, height = 4,
       filename = 'Results/MouseUrine/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_0.1_PLAUR_Clusters_less.pdf')

# NECTIN1+ EV Counts
sub.seuOb_NECTIN1 <- subset(seuOb, NECTIN1 > 0)
a <- table(sub.seuOb_NECTIN1@meta.data$orig.ident, sub.seuOb_NECTIN1@meta.data$SCT_snn_res.0.1)
write.csv(a, 'Results/MouseUrine/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_0.1_NECTIN1_subCluster.csv')
dt <- a %>% as.data.frame()
names(dt) <- c('Samples', 'Cluster', 'EV')
dt <- dt %>% filter(Cluster %in% 0:20)
dt$Group <- ifelse(str_detect(dt$Samples, 'D'), 'AD', 'Con')
dt$Cluster <- factor(dt$Cluster, levels = 0:20)
dt <- dt %>% filter(Cluster %in% c(0,1,3,12,20))
p <- ggplot(dt, aes(x = Group, y = EV, fill = Group)) + theme_bw() +
  geom_boxplot(width = 0.8, size = 0.6) +
  geom_jitter(width = 0.3, size = 0.5, alpha = 0.8) +
  geom_signif(comparisons = list(c("AD","Con")), 
              map_signif_level = T, step_increase = 0.05, tip_length = 0.01, 
              test = "wilcox.test", textsize = 2.5) +
  xlab('') + ylab('EV Counts') + ggtitle('NECTIN1') +
  scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
  theme(legend.position = "none") + facet_grid(~Cluster) + 
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.border = element_blank())
ggsave(plot = p, width = 4, height = 4,
       filename = 'Results/MouseUrine/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_0.1_NECTIN1_Clusters_less.png')
ggsave(plot = p, width = 4, height = 4,
       filename = 'Results/MouseUrine/MouseUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_0.1_NECTIN1_Clusters_less.pdf')


# Similar codes have been omitted




