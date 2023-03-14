# -- ------------------------------------------------------------------------------------------------------------------
getwd()
# 设置最大的内存
if(.Platform$OS.type == "windows") withAutoprint({
  memory.size()
  memory.size(TRUE)
  memory.limit()
})
# options(future.globals.maxSize = 30*1024^3) #设置为30G
setwd("D:/2.AD.PBA.CAM.Article.V3.Proj")

library(tidyverse) # 处理数据神包
library(Seurat) # 单细胞数据分析包
library(harmony) # 整合数据

# -- ------------------------------------------------------------------------------------------------------------------
# 1.因为原始数据每个样本的EV数非常多，且非常稀疏，所以将n个EV合成为1个super EV来进行Seurat分析。Per5/10/20/30
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
# Urine和Tear用Per10
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID <- sampleID %>% filter(Source != 'Plasma') %>% filter(Source != 'Saliva')
pseudo.size = 10 # 每10个EV合成一个EV
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('RawData/',sampleID[i,1],'.total_ev_protein.csv')
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
    fn2 <- str_c('pseudoData/pseudoPer10/',sampleID[i,1],'.sEv10.csv')
    write.csv(pseudo.data[-1], fn2)
    rm(data, wide_data,raw.data,pseudo.data)
  }
}
# Saliva Per20
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID <- sampleID %>% filter(Source == 'Saliva')
pseudo.size = 20 # 每20个EV合成一个EV
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('RawData/',sampleID[i,1],'.total_ev_protein.csv')
  data <- read.csv(fn1, header = T)
  wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
  rownames(wide_data) <- wide_data$ev
  wide_data <- wide_data[, -c(1:2)]
  superEv20 <- floor(nrow(wide_data)/pseudo.size) 
  if (superEv20 > 0) {
    set.seed(2022)
    random.row <- sample(1:nrow(wide_data), superEv20*pseudo.size, replace = F)
    raw.data <- wide_data[random.row,] %>% as.data.frame()
    raw.data$Group <- rep(c(1:superEv20), each = pseudo.size)
    raw.data <- relocate(raw.data, Group)
    pseudo.data <- aggregate(list(raw.data[2:ncol(raw.data)]), by = list(raw.data$Group), FUN = sum)
    rownames(pseudo.data) <- str_c(sampleID[i,2], 'sEV', 1:superEv20, sep = '_')
    fn2 <- str_c('pseudoData/pseudoPer20/',sampleID[i,1],'.sEv20.csv')
    write.csv(pseudo.data[-1], fn2)
    rm(data, wide_data,raw.data,pseudo.data)
  }
}
# Plasma Per30
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID <- sampleID %>% filter(Source == 'Plasma')
pseudo.size = 30 # 每30个EV合成一个EV
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('RawData/',sampleID[i,1],'.total_ev_protein.csv')
  data <- read.csv(fn1, header = T)
  wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
  rownames(wide_data) <- wide_data$ev
  wide_data <- wide_data[, -c(1:2)]
  superEv30 <- floor(nrow(wide_data)/pseudo.size) 
  if (superEv30 > 0) {
    set.seed(2022)
    random.row <- sample(1:nrow(wide_data), superEv30*pseudo.size, replace = F)
    raw.data <- wide_data[random.row,] %>% as.data.frame()
    raw.data$Group <- rep(c(1:superEv30), each = pseudo.size)
    raw.data <- relocate(raw.data, Group)
    pseudo.data <- aggregate(list(raw.data[2:ncol(raw.data)]), by = list(raw.data$Group), FUN = sum)
    rownames(pseudo.data) <- str_c(sampleID[i,2], 'sEV', 1:superEv30, sep = '_')
    fn2 <- str_c('pseudoData/pseudoPer30/',sampleID[i,1],'.sEv30.csv')
    write.csv(pseudo.data[-1], fn2)
    rm(data, wide_data,raw.data,pseudo.data)
  }
}

# 由于之前已经跑过这部分的内容了，所以直接把之前的结果COPY过来
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID <- sampleID %>% filter(Source != 'Plasma') %>% filter(Source != 'Saliva')
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('D:/2.AD.PBA.CAM.Article.V1.Proj/pseudoData/pseudoPer10/',sampleID[i,1],'.sEv10.csv')
  file.copy(from = fn1, to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/pseudoData/pseudoPer10/')
}
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID <- sampleID %>% filter(Source == 'Saliva')
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('D:/2.AD.PBA.CAM.Article.V1.Proj/pseudoData/pseudoPer20/',sampleID[i,1],'.sEv20.csv')
  file.copy(from = fn1, to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/pseudoData/pseudoPer20/')
}
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID <- sampleID %>% filter(Source == 'Plasma')
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('D:/2.AD.PBA.CAM.Article.V1.Proj/pseudoData/pseudoPer30/',sampleID[i,1],'.sEv30.csv')
  file.copy(from = fn1, to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/pseudoData/pseudoPer30/')
}

# -- ------------------------------------------------------------------------------------------------------------------
# 2.创建Seurat Objects
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
# Urine Tear
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID <- sampleID %>% filter(Source != 'Plasma') %>% filter(Source != 'Saliva')
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('pseudoData/pseudoPer10/',sampleID[i,1],'.sEv10.csv')
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
    fn2 <- str_c('SeuratObjects/pseudoPer10/',sampleID[i,1],'.pseudoPer10.rds')
    write_rds(seuOb, fn2)
    rm(data, seuOb)
  }
}
# Saliva
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID <- sampleID %>% filter(Source == 'Saliva')
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('pseudoData/pseudoPer20/',sampleID[i,1],'.sEv20.csv')
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
    fn2 <- str_c('SeuratObjects/pseudoPer20/',sampleID[i,1],'.pseudoPer20.rds')
    write_rds(seuOb, fn2)
    rm(data, seuOb)
  }
}
# Plasma
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID <- sampleID %>% filter(Source == 'Plasma')
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('pseudoData/pseudoPer30/',sampleID[i,1],'.sEv30.csv')
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
    fn2 <- str_c('SeuratObjects/pseudoPer30/',sampleID[i,1],'.pseudoPer30.rds')
    write_rds(seuOb, fn2)
    rm(data, seuOb)
  }
}
# 由于之前已经跑过这部分的内容了，所以直接把之前的结果COPY过来
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID <- sampleID %>% filter(Source != 'Plasma') %>% filter(Source != 'Saliva')
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/pseudoPer10/',sampleID[i,1],'.pseudoPer10.rds')
  file.copy(from = fn1, to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/pseudoPer10/')
}
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID <- sampleID %>% filter(Source == 'Saliva')
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/pseudoPer20/',sampleID[i,1],'.pseudoPer20.rds')
  file.copy(from = fn1, to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/pseudoPer20/')
}
sampleID <- read.csv('Human_SampleID_SampleName.csv')
sampleID <- sampleID %>% filter(Source == 'Plasma')
for (i in 1:nrow(sampleID)) {
  gc()
  print(sampleID[i,2])
  fn1 <- str_c('D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/pseudoPer30/',sampleID[i,1],'.pseudoPer30.rds')
  file.copy(from = fn1, to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/pseudoPer30/')
}

# -- ------------------------------------------------------------------------------------------------------------------
# 3.将SeuratObjects merge在一起
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
sampleID <- read_csv('Human_SampleID_SampleName.csv', show_col_types = F)
# HumanPlasma
sampleID.sub <- sampleID %>% filter(Source == 'Plasma')
filenames <- str_c('SeuratObjects/pseudoPer30/', sampleID.sub$SampleID,'.pseudoPer30.rds')
seuList <- lapply(filenames, readRDS)
seuOb <- merge(seuList[[1]], seuList[-1], project = 'pseudoPer30')
head(seuOb@meta.data)
table(seuOb@meta.data$orig.ident)
sampleID <- str_split(rownames(seuOb@meta.data), '_s') %>% as.data.frame() %>% t() %>% as.data.frame()
seuOb@meta.data$orig.ident <- sampleID$V1
write_rds(seuOb, 'SeuratObjects/HumanPlasma.AllSamples.pseudoPer30.Merge.rds')
# HumanSaliva
sampleID.sub <- sampleID %>% filter(Source == 'Saliva')
filenames <- str_c('SeuratObjects/pseudoPer20/', sampleID.sub$SampleID,'.pseudoPer20.rds')
seuList <- lapply(filenames, readRDS)
seuOb <- merge(seuList[[1]], seuList[-1], project = 'pseudoPer20')
head(seuOb@meta.data)
table(seuOb@meta.data$orig.ident)
sampleID <- str_split(rownames(seuOb@meta.data), '_s') %>% as.data.frame() %>% t() %>% as.data.frame()
seuOb@meta.data$orig.ident <- sampleID$V1
write_rds(seuOb, 'SeuratObjects/HumanSaliva.AllSamples.pseudoPer20.Merge.rds')
# HumanUrine
sampleID.sub <- sampleID %>% filter(Source == 'Urine')
filenames <- str_c('SeuratObjects/pseudoPer10/', sampleID.sub$SampleID,'.pseudoPer10.rds')
seuList <- lapply(filenames, readRDS)
seuOb <- merge(seuList[[1]], seuList[-1], project = 'pseudoPer10')
head(seuOb@meta.data)
table(seuOb@meta.data$orig.ident)
sampleID <- str_split(rownames(seuOb@meta.data), '_s') %>% as.data.frame() %>% t() %>% as.data.frame()
seuOb@meta.data$orig.ident <- sampleID$V1
write_rds(seuOb, 'SeuratObjects/HumanUrine.AllSamples.pseudoPer10.Merge.rds')
# HumanTear
sampleID.sub <- sampleID %>% filter(Source == 'Tear')
filenames <- str_c('SeuratObjects/pseudoPer10/', sampleID.sub$SampleID,'.pseudoPer10.rds')
seuList <- lapply(filenames, readRDS)
seuOb <- merge(seuList[[1]], seuList[-1], project = 'pseudoPer10')
head(seuOb@meta.data)
table(seuOb@meta.data$orig.ident)
sampleID <- str_split(rownames(seuOb@meta.data), '_s') %>% as.data.frame() %>% t() %>% as.data.frame()
seuOb@meta.data$orig.ident <- sampleID$V1
write_rds(seuOb, 'SeuratObjects/HumanTear.AllSamples.pseudoPer10.Merge.rds')

file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/HumanPlasma.AllSamples.pseudoPer30.Merge.rds',
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/HumanSaliva.AllSamples.pseudoPer20.Merge.rds', 
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/HumanTear.AllSamples.pseudoPer10.Merge.rds', 
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/HumanUrine.AllSamples.pseudoPer10.Merge.rds', 
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/')

# -- ------------------------------------------------------------------------------------------------------------------
# 4.SCTransform标准化
# -- ------------------------------------------------------------------------------------------------------------------
# return.only.var.genes = TRUE ，只对差异基因进行scale的结果。设置为FALSE相当于三步法的ScaleData(features=all.genes)。
# min_cells = 5是SCTransform默认的过滤基因的参数
rm(list=ls())
gc()
filenames <- Sys.glob('SeuratObjects/*.Merge.rds')
for (i in filenames) {
  seuOb <- readRDS(i)
  seuOb <- seuOb %>% SCTransform(return.only.var.genes = F, verbose = F)
  fn <- str_c(str_remove(i, '.rds'), '_SCT.rds')
  write_rds(seuOb, fn)
}
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/HumanPlasma.AllSamples.pseudoPer30.Merge_SCT.rds', 
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/HumanSaliva.AllSamples.pseudoPer20.Merge_SCT.rds', 
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/HumanTear.AllSamples.pseudoPer10.Merge_SCT.rds', 
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/HumanUrine.AllSamples.pseudoPer10.Merge_SCT.rds', 
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/')

# -- ------------------------------------------------------------------------------------------------------------------
# 5.Harmony整合
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
filenames <- Sys.glob('SeuratObjects/*.Merge_SCT.rds')
for (i in filenames) {
  seuOb <- readRDS(i)
  seuOb <- RunPCA(seuOb, assay = "SCT", verbose = F) %>% 
    RunHarmony(assay.use = "SCT", reduction = "pca", group.by.vars = 'orig.ident', plot_convergence = T)
  fn1 <- str_c(str_remove(i, '.rds'), '_Harmony.Convergence.png')
  ggsave(fn1, width = 6, height = 4)
  fn2 <- str_c(str_remove(i, '.rds'), '_Harmony.rds')
  write_rds(seuOb, fn2)
}
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/HumanPlasma.AllSamples.pseudoPer30.Merge_SCT_Harmony.rds', 
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/HumanSaliva.AllSamples.pseudoPer20.Merge_SCT_Harmony.rds', 
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/HumanTear.AllSamples.pseudoPer10.Merge_SCT_Harmony.rds', 
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/')
file.copy(from = 'D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/HumanUrine.AllSamples.pseudoPer10.Merge_SCT_Harmony.rds', 
          to = 'D:/2.AD.PBA.CAM.Article.V3.Proj/SeuratObjects/')
seuOb <- readRDS('D:/2.AD.PBA.CAM.Article.V1.Proj/SeuratObjects/HumanSaliva.AllSamples.pseudoPer20.Merge_SCT_Harmony.rds')
table(seuOb@meta.data$orig.ident)
# -- ------------------------------------------------------------------------------------------------------------------
# 6.降维UMAP/TSNE，这一步还是用所有的EV 比较好，因为最后需要获得所有的avrExp，并且用部分EV的话得到的降维坐标会不一样
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
filenames <- Sys.glob('SeuratObjects/*.Merge_SCT_Harmony.rds')
for (i in filenames) {
  seuOb <- readRDS(i)
  seuOb <- RunUMAP(seuOb, assay = "SCT", reduction = "harmony", dims = 1:20, verbose = F) %>% 
    RunTSNE(assay = "SCT", reduction = "harmony", dims = 1:20, check_duplicates = F) %>% 
    FindNeighbors(assay = "SCT", reduction = "harmony", dims = 1:20, verbose = F)
  for (a in c(0.1, 0.3)){
    seuOb <- FindClusters(seuOb, resolution = a, verbose = F)
  }
  fn <- str_c(str_remove(i, '.rds'), '.Dim1-20_Res0.3.rds')
  write_rds(seuOb, fn)
}
# -- ------------------------------------------------------------------------------------------------------------------
# 7.可视化的时候只抽取一部分EV来作图
# -- ------------------------------------------------------------------------------------------------------------------
# HumanUrine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
seuOb <- readRDS('SeuratObjects/HumanUrine.AllSamples.pseudoPer10.Merge_SCT_Harmony.Dim1-20_Res0.3.rds')
Idents(seuOb) <- 'orig.ident'
set.seed(2022)
sub.seuOb <- subset(seuOb, downsample = 1000)
Idents(sub.seuOb) <- 'SCT_snn_res.0.1'
p1 <- DimPlot(sub.seuOb, reduction = "umap", label = TRUE, raster = F, repel = T) + NoLegend()
p2 <- DimPlot(sub.seuOb, reduction = "tsne", label = TRUE, raster = F, repel = T)
p3 <- p1 + p2
ggsave(p3, filename = 'Results/Figure6/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.1.1000.png', 
       height = 5, width = 11)
ggsave(p3, filename = 'Results/Figure6/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.1.1000.pdf', 
       height = 5, width = 11)
Idents(seuOb) <- 'SCT_snn_res.0.1'
seuOb.markers <- FindAllMarkers(seuOb, assay = 'RNA', min.pct = 0.5, verbose = F)
top5 <- seuOb.markers %>% group_by(cluster) %>% top_n(n = 5, wt = -p_val_adj)
p4 <- DoHeatmap(sub.seuOb, features = top5$gene) + NoLegend()
ggsave(p4, file = 'Results/Figure6/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.1.1000.clusters.markers.Top5.DoHeatmap.png',
       height = 10, width = 10)
p5 <- DotPlot(sub.seuOb, features = unique(top5$gene), cols = c("lightgrey", "red")) + RotatedAxis()
ggsave(p5, file = 'Results/Figure6/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.1.1000.clusters.markers.Top5.DotPlot.png',
       height = 8, width = 8)
ggsave(p5, file = 'Results/Figure6/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.1.1000.clusters.markers.Top5.DotPlot.pdf',
       height = 8, width = 10)

# -- ------------------------------------------------------------------------------------------------------------------
# HumanPlasma
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
seuOb <- readRDS('SeuratObjects/HumanPlasma.AllSamples.pseudoPer30.Merge_SCT_Harmony.Dim1-20_Res0.3.rds')
Idents(seuOb) <- 'orig.ident'
set.seed(2022)
sub.seuOb <- subset(seuOb, downsample = 1000)
Idents(sub.seuOb) <- 'SCT_snn_res.0.3'
p1 <- DimPlot(sub.seuOb, reduction = "umap", label = TRUE, raster = F, repel = T) + NoLegend()
p2 <- DimPlot(sub.seuOb, reduction = "tsne", label = TRUE, raster = F, repel = T)
p3 <- p1 + p2
ggsave(p3, filename = 'Results/Figure6/HumanPlasma.pseudoPer30.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.3.1000.png', 
       height = 5, width = 11)
ggsave(p3, filename = 'Results/Figure6/HumanPlasma.pseudoPer30.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.3.1000.pdf', 
       height = 5, width = 11)
Idents(seuOb) <- 'SCT_snn_res.0.3'
seuOb.markers <- FindAllMarkers(seuOb, assay = 'RNA', min.pct = 0.5, verbose = F)
top5 <- seuOb.markers %>% group_by(cluster) %>% top_n(n = 5, wt = -p_val_adj)
top5$cluster <- top5$cluster %>% factor(levels = 0:23)
top5 <- top5 %>% arrange(cluster)
sub.seuOb@meta.data %>% head()
sub.seuOb@active.ident %>% head()
sub.seuOb@active.ident <- sub.seuOb@active.ident %>% factor(levels = c(0:23))
p5 <- DotPlot(sub.seuOb, features = unique(top5$gene), cols = c("lightgrey", "red")) + RotatedAxis()
ggsave(p5, file = 'Results/Figure6/HumanPlasma.pseudoPer30.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.3.1000.clusters.markers.Top5.DotPlot.png',
       height = 8, width = 8)
ggsave(p5, file = 'Results/Figure6/HumanPlasma.pseudoPer30.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.3.1000.clusters.markers.Top5.DotPlot.pdf',
       height = 8, width = 10)

# -- ------------------------------------------------------------------------------------------------------------------
# HumanTear
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
seuOb <- readRDS('SeuratObjects/HumanTear.AllSamples.pseudoPer10.Merge_SCT_Harmony.Dim1-20_Res0.3.rds')
Idents(seuOb) <- 'orig.ident'
set.seed(2022)
sub.seuOb <- subset(seuOb, downsample = 1000)
Idents(sub.seuOb) <- 'SCT_snn_res.0.3'
sub.seuOb@meta.data %>% head()
table(sub.seuOb$SCT_snn_res.0.3)
p1 <- DimPlot(sub.seuOb, reduction = "umap", label = TRUE, raster = F, repel = T) + NoLegend()
p2 <- DimPlot(sub.seuOb, reduction = "tsne", label = TRUE, raster = F, repel = T)
p3 <- p1 + p2
ggsave(p3, filename = 'Results/Figure6/HumanTear.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.3.1000.png',
       height = 5, width = 11)
ggsave(p3, filename = 'Results/Figure6/HumanTear.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.3.1000.pdf',
       height = 5, width = 11)
Idents(seuOb) <- 'SCT_snn_res.0.3'
seuOb.markers <- FindAllMarkers(seuOb, assay = 'RNA', min.pct = 0.5, verbose = F)
write.csv(seuOb.markers, 'Results/Figure6/Files/HumanTear.pseudoPer10.Merge_SCT_Harmony.Dim1-20_Res0.3.clusters.markers.csv')
top5 <- seuOb.markers %>% group_by(cluster) %>% top_n(n = 5, wt = -p_val_adj)
top5$cluster <- top5$cluster %>% factor(levels = 0:18)
top5 <- top5 %>% arrange(cluster)
sub.seuOb@meta.data %>% head()
sub.seuOb@active.ident %>% head()
sub.seuOb@active.ident <- sub.seuOb@active.ident %>% factor(levels = c(0:19))
p5 <- DotPlot(sub.seuOb, features = unique(top5$gene), cols = c("lightgrey", "red")) + RotatedAxis()
ggsave(p5, file = 'Results/Figure6/HumanTear.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.3.1000.clusters.markers.Top5.DotPlot.png',
       height = 8, width = 8)
ggsave(p5, file = 'Results/Figure6/HumanTear.pseudoPer10.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.3.1000.clusters.markers.Top5.DotPlot.pdf',
       height = 8, width = 10)

# -- ------------------------------------------------------------------------------------------------------------------
# HumanSaliva
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
seuOb <- readRDS('SeuratObjects/HumanSaliva.AllSamples.pseudoPer20.Merge_SCT_Harmony.Dim1-20_Res0.3.rds')
table(seuOb$orig.ident)
Idents(seuOb) <- 'orig.ident'
set.seed(2022)
sub.seuOb <- subset(seuOb, downsample = 1000)
Idents(sub.seuOb) <- 'SCT_snn_res.0.1'
p1 <- DimPlot(sub.seuOb, reduction = "umap", label = TRUE, raster = F, repel = T) + NoLegend()
p2 <- DimPlot(sub.seuOb, reduction = "tsne", label = TRUE, raster = F, repel = T)
p3 <- p1 + p2
ggsave(p3, filename = 'Results/Figure6/HumanSaliva.pseudoPer20.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.1.1000.png', 
       height = 5, width = 11)
ggsave(p3, filename = 'Results/Figure6/HumanSaliva.pseudoPer20.Merge_SCT_Harmony.Dim1-20_UMAP_tSNE_0.1.1000.pdf', 
       height = 5, width = 11)
Idents(seuOb) <- 'SCT_snn_res.0.1'
seuOb.markers <- FindAllMarkers(seuOb, assay = 'RNA', min.pct = 0.5, verbose = F)
write.csv(seuOb.markers, 'Results/Figure6/Files/HumanSaliva.pseudoPer20.Merge_SCT_Harmony.Dim1-20_Res0.1.clusters.markers.csv')
top5 <- seuOb.markers %>% group_by(cluster) %>% top_n(n = 5, wt = -p_val_adj)
top5$cluster <- top5$cluster %>% factor(levels = 0:23)
top5 <- top5 %>% arrange(cluster)
sub.seuOb@meta.data %>% head()
sub.seuOb@active.ident %>% head()
sub.seuOb@active.ident <- sub.seuOb@active.ident %>% factor(levels = c(0:23))
p5 <- DotPlot(sub.seuOb, features = unique(top5$gene), cols = c("lightgrey", "red")) + RotatedAxis()
ggsave(p5, file = 'Results/Figure6/HumanSaliva.pseudoPer20.Merge_SCT_Harmony.Dim1-20_Res0.1.1000.clusters.markers.Top5.DotPlot.png',
       height = 8, width = 8)
ggsave(p5, file = 'Results/Figure6/HumanSaliva.pseudoPer20.Merge_SCT_Harmony.Dim1-20_Res0.1.1000.clusters.markers.Top5.DotPlot.pdf',
       height = 8, width = 10)

# -- ------------------------------------------------------------------------------------------------------------------
# 8.对每个亚群提取平均表达值
# -- ------------------------------------------------------------------------------------------------------------------
# HumanUrine
dir.create('Results/Figure6/HumanUrine')
dir.create('Results/Figure6/HumanUrine/Expression')
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
seuOb <- readRDS('SeuratObjects/HumanUrine.AllSamples.pseudoPer10.Merge_SCT_Harmony.Dim1-20_Res0.3.rds')
# 按照亚群来分析差异表达
seuOb@meta.data$orig.ident %>% table()
DefaultAssay(seuOb) #查看默认数据集
Idents(seuOb) <- 'SCT_snn_res.0.1' 
table(seuOb@active.ident)
for (i in 0:20) {
  gc()
  sub_seuOb <- subset(seuOb, seurat_clusters == i)
  avrExp <- AverageExpression(sub_seuOb,group.by = 'orig.ident')
  fn1 <- str_c('Results/Figure6/HumanUrine/Expression/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster',i,'.Expression.SCT.csv')
  write.csv(avrExp[['SCT']], fn1)
  fn2 <- str_c('Results/Figure6/HumanUrine/Expression/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster',i,'.Expression.RNA.csv')
  write.csv(avrExp[['RNA']], fn2)
}
# -- ------------------------------------------------------------------------------------------------------------------
# HumanSaliva
dir.create('Results/Figure6/HumanSaliva')
dir.create('Results/Figure6/HumanSaliva/Expression')
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
seuOb <- readRDS('SeuratObjects/HumanSaliva.AllSamples.pseudoPer20.Merge_SCT_Harmony.Dim1-20_Res0.3.rds')
# 按照亚群来分析差异表达
seuOb@meta.data$orig.ident %>% table()
DefaultAssay(seuOb) #查看默认数据集
Idents(seuOb) <- 'SCT_snn_res.0.1' 
table(seuOb@active.ident)
for (i in 0:23) {
  gc()
  sub_seuOb <- subset(seuOb, seurat_clusters == i)
  avrExp <- AverageExpression(sub_seuOb,group.by = 'orig.ident')
  fn1 <- str_c('Results/Figure6/HumanSaliva/Expression/HumanSaliva.pseudoPer20.Merge_SCT_Harmony.Dim1-20.0.1.Cluster',i,'.Expression.SCT.csv')
  write.csv(avrExp[['SCT']], fn1)
  fn2 <- str_c('Results/Figure6/HumanSaliva/Expression/HumanSaliva.pseudoPer20.Merge_SCT_Harmony.Dim1-20.0.1.Cluster',i,'.Expression.RNA.csv')
  write.csv(avrExp[['RNA']], fn2)
}
# -- ------------------------------------------------------------------------------------------------------------------
# HumanTear
dir.create('Results/Figure6/HumanTear')
dir.create('Results/Figure6/HumanTear/Expression')
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
seuOb <- readRDS('SeuratObjects/HumanTear.AllSamples.pseudoPer10.Merge_SCT_Harmony.Dim1-20_Res0.3.rds')
# 按照亚群来分析差异表达
seuOb@meta.data$orig.ident %>% table()
DefaultAssay(seuOb) #查看默认数据集
Idents(seuOb) <- 'SCT_snn_res.0.3' 
table(seuOb@active.ident)
for (i in 0:19) {
  gc()
  sub_seuOb <- subset(seuOb, seurat_clusters == i)
  avrExp <- AverageExpression(sub_seuOb,group.by = 'orig.ident')
  fn1 <- str_c('Results/Figure6/HumanTear/Expression/HumanTear.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.3.Cluster',i,'.Expression.SCT.csv')
  write.csv(avrExp[['SCT']], fn1)
  fn2 <- str_c('Results/Figure6/HumanTear/Expression/HumanTear.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.3.Cluster',i,'.Expression.RNA.csv')
  write.csv(avrExp[['RNA']], fn2)
}
# -- ------------------------------------------------------------------------------------------------------------------
# HumanPlasma
dir.create('Results/Figure6/HumanPlasma')
dir.create('Results/Figure6/HumanPlasma/Expression')
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
seuOb <- readRDS('SeuratObjects/HumanPlasma.AllSamples.pseudoPer30.Merge_SCT_Harmony.Dim1-20_Res0.3.rds')
# 按照亚群来分析差异表达
seuOb@meta.data$orig.ident %>% table()
DefaultAssay(seuOb) #查看默认数据集
Idents(seuOb) <- 'SCT_snn_res.0.3' 
table(seuOb@active.ident)
for (i in 0:24) {
  gc()
  sub_seuOb <- subset(seuOb, seurat_clusters == i)
  avrExp <- AverageExpression(sub_seuOb,group.by = 'orig.ident')
  fn1 <- str_c('Results/Figure6/HumanPlasma/Expression/HumanPlasma.pseudoPer30.Merge_SCT_Harmony.Dim1-20.0.3.Cluster',i,'.Expression.SCT.csv')
  write.csv(avrExp[['SCT']], fn1)
  fn2 <- str_c('Results/Figure6/HumanPlasma/Expression/HumanPlasma.pseudoPer30.Merge_SCT_Harmony.Dim1-20.0.3.Cluster',i,'.Expression.RNA.csv')
  write.csv(avrExp[['RNA']], fn2)
}

# -- ------------------------------------------------------------------------------------------------------------------
# 9.对每个亚群进行caret机器学习
library(tidyverse) # 处理数据神包
library(glmnet)
library(caret)
library(gbm)
library(pROC)
# -- ------------------------------------------------------------------------------------------------------------------
# HumanSaliva #其他组见Step14
dir.create('Results/Figure6/HumanSaliva/Logistic')
dir.create('Results/Figure6/HumanSaliva/caret')
dir.create('Results/Figure6/HumanSaliva/caret/RDS')
# -- ------------------------------------------------------------------------------------------------------------------
# a.用Logidtic回归筛选变量
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
for( j in 0:23) {
  fn1 <- str_c('Results/Figure6/HumanSaliva/Expression/HumanSaliva.pseudoPer20.Merge_SCT_Harmony.Dim1-20.0.1.Cluster',j,'.Expression.SCT.csv')
  data <- read.csv(fn1, row.names = 1)
  data <- data %>% t() %>% as.data.frame()
  data$ID <- rownames(data) %>% factor
  data$Group <- ifelse(str_detect(rownames(data), 'D'),'AD','Con') %>% factor(levels = c('AD','Con'))
  data <- relocate(data, ID, Group)
  # 划分训练集与测试集:
  data.ID <- levels(data$ID) # 只有factor才能用levels这个函数
  data.ID <- data.frame(ID = data.ID)
  data.ID$Group <- ifelse(str_detect(data.ID$ID, 'D'), 'AD', 'Con')
  # set.seed(1234)
  # set.seed(5678)
  # set.seed(2020)
  # set.seed(2021)
  set.seed(2022)
  train.ID <- createDataPartition(data.ID$Group, p = 0.6, list = FALSE)
  train.ID <- data.ID$ID[train.ID]
  data.train <- data %>% filter(ID %in% train.ID) %>% dplyr::select(!ID)
  data.test <- data %>% filter(!ID %in% train.ID) %>% dplyr::select(!ID)
  # 去掉全为0的列：
  data.train <- data.train %>% select_if(is.numeric) %>% select_if(~sum(.) > 0)
  data.train$Group <- ifelse(str_detect(rownames(data.train), 'D'),'AD','Con') %>% factor(levels = c('AD','Con'))
  data.train <- relocate(data.train, Group)
  my.var <- names(data.train)
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
  res <- res %>% arrange(`Pr(>|z|)`)
  fn <- str_c('Results/Figure6/HumanSaliva/Logistic/AllSamples.pseudoPer20.Res0.1.Cluster',j,'_Train.2022_Logistic.csv')
  write.csv(res, fn)
  myPro <- rownames(res)[1:3]
  data.train <- data.train %>% select(Group, all_of(myPro))
  data.test <- data.test %>% select(Group, all_of(myPro))
  # 将数据进行标准化并补足缺失值，这时可以用preProcess命令，缺省参数是标准化数据。?preProcess
  preProcValues <- preProcess(data.train, method = c('center','scale')) # 用于进行其他验证数据集的标准化
  train.transformed <- predict(preProcValues, data.train)
  test.transformed <- predict(preProcValues, data.test)
  # 保存数据
  fn2 <- str_c('Results/Figure6/HumanSaliva/caret/RDS/AllSamples.pseudoPer20.Res.0.1.Cluster',j,'_Train.2022_LogisticTop3_caret.RData')
  save(data.train, data.test, preProcValues, train.transformed, test.transformed,
       file = fn2)
}
# -- ------------------------------------------------------------------------------------------------------------------
# b.模型训练及调参
# -- ------------------------------------------------------------------------------------------------------------------
# 常用的分类算法包括：NBC（Naive Bayesian Classifier，朴素贝叶斯分类）算法、
# LR（Logistic Regress，逻辑回归）算法、C5.0 决策树算法、SVM（Support Vector Machine，支持向量机）算法、
# KNN(K-Nearest Neighbor，K 最近邻近)算法、ANN（Artificial Neural Network，人工神经网络）算法
# Stochastic Gradient Boosting(GBM)
rm(list=ls())
gc()
myMethod <- c('glmnet','kknn','C5.0','rf','AdaBag','nnet','svmLinear','svmPoly','svmRadial','nb') #gbm方法报错了
fitControl = trainControl(method = 'cv', number = 5,
                          classProbs = T, summaryFunction = twoClassSummary,
                          search = 'random') # 随机调参
for( j in 0:23) {
  fn2 <- str_c('Results/Figure6/HumanSaliva/caret/RDS/AllSamples.pseudoPer20.Res.0.1.Cluster',j,'_Train.2022_LogisticTop3_caret.RData')
  load(file = fn2)
  for (i in 1:length(myMethod)) {
    set.seed(2022)
    Fit <- train(Group ~ ., data = train.transformed, method = myMethod[i], # 不调参的话可以使用默认参数，只用这一行代码就行
                 trControl = fitControl, 
                 tuneLength = 10, # 随机调参的随机次数
                 metric = 'ROC', verbose = F) 
    fn3 <- str_c('Results/Figure6/HumanSaliva/caret/RDS/pseudoPer20.Res0.1.Cluster',j,'_Train.2022.', myMethod[i], '.Fit.rds')
    saveRDS(Fit, fn3)
  }
}

# -- ------------------------------------------------------------------------------------------------------------------
# c.模型预测及评价
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
for ( j in 0:23) {
  fn2 <- str_c('Results/Figure6/HumanSaliva/caret/RDS/AllSamples.pseudoPer20.Res.0.1.Cluster',j,'_Train.2022_LogisticTop3_caret.RData')
  load(file = fn2)
  fn4 <- str_c('Results/Figure6/HumanSaliva/caret/RDS/pseudoPer20.Res0.1.Cluster',j,'_Train.2022.*.Fit.rds')
  filenames <- Sys.glob(fn4)
  fn5 <- str_c('Results/Figure6/HumanSaliva/caret/RDS/pseudoPer20.Res0.1.Cluster',j,'_Train.2022.')
  myMethod <- filenames %>% str_remove(fn5) %>% str_remove('.Fit.rds')
  res <- data.frame()
  for ( a in 1:10) {
    Fit <- readRDS(filenames[a])
    predict.train <- predict(Fit, newdata = train.transformed)
    predict.test <- predict(Fit, newdata = test.transformed)
    # ROC分析的时候要区分方向，否则容易得到相反的AUC
    roc.train <- roc(train.transformed$Group, as.numeric(predict.train), levels = c('Con','AD'), direction = '>')$auc
    roc.test <- roc(test.transformed$Group, as.numeric(predict.test), levels = c('Con','AD'), direction = '>')$auc
    Fit.roc <- data.frame(Train = roc.train, Test = roc.test) %>% mutate(method = myMethod[a])
    rownames(Fit.roc) <- myMethod[a]
    res <- rbind(Fit.roc, res)
    fn7 <- str_c('Results/Figure6/HumanSaliva/caret/RDS/pseudoPer20.Res0.1.Cluster',j,'_Train.2022.caret.ROC.csv')
    write.csv(res, fn7)
  }
}
# -- ------------------------------------------------------------------------------------------------------------------
# d.将各个亚群的合并成一个总的caret.ROC
# -- ------------------------------------------------------------------------------------------------------------------
filenames <- Sys.glob('Results/Figure6/HumanSaliva/caret/RDS/pseudoPer20.Res0.1.Cluster*_Train.2022.caret.ROC.csv')
clusters <- filenames %>% str_remove('Results/Figure6/HumanSaliva/caret/RDS/pseudoPer20.Res0.1.Cluster') %>%
  str_remove('_Train.2022.caret.ROC.csv')
data <- read.csv(filenames[1])
data$Cluster <- clusters[1]
for (i in 2:length(clusters)) {
  temp <- read.csv(filenames[i])
  temp$Cluster <- clusters[i]
  data <- rbind(data, temp)
}
write_csv(data, file = 'Results/Figure6/Files/HumanSaliva.AllSamples.pseudoPer20.Res0.1.subClusters_Train.2022.caret.ROC.csv')
data <- data %>% filter(Test >= 0.65)
write_csv(data, file = 'Results/Figure6/Files/HumanSaliva.AllSamples.pseudoPer20.Res0.1.subClusters_Train.2022.caret.ROC_0.65.csv')


# -- ------------------------------------------------------------------------------------------------------------------
# 10.根据Tang的要求对Train的数据进行ROC分析，并做1-5个蛋白的联合ROC的AUC趋势图
# -- ------------------------------------------------------------------------------------------------------------------
dir.create('Results/Figure6')
dir.create('Results/Figure6/Files')
library(tidyverse) # 处理数据神包
library(glmnet) # 进行正则技术所需
library(pROC) # ROC分析与绘图
library(ggsci) # 期刊配色，单个离散配色模块中最多只有51种颜色
library(scales) # show_col
library(patchwork) # 拼接图
library(ggbreak) # 坐标轴截断
library(ggsignif) # 计算差异显著性并添加显著性线段
library(ggthemes) # 主题设置
library(pheatmap) # 热图
# -- ------------------------------------------------------------------------------------------------------------------
# HumanUrine #其他组没做
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
for( j in 0:20) {
  fn1 <- str_c('Results/HumanUrine/subCluster/Expression/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster',j,'.Expression.SCT.csv')
  data <- read.csv(fn1, row.names = 1)
  data <- data %>% t() %>% as.data.frame()
  data$ID <- rownames(data) %>% factor
  data$Group <- ifelse(str_detect(rownames(data), 'D'),'AD','Con') %>% factor(levels = c('AD','Con'))
  data <- relocate(data, ID, Group)
  # 划分训练集与测试集:
  data.ID <- levels(data$ID) # 只有factor才能用levels这个函数
  data.ID <- data.frame(ID = data.ID)
  data.ID$Group <- ifelse(str_detect(data.ID$ID, 'D'), 'AD', 'Con')
  set.seed(2022)
  train.ID <- createDataPartition(data.ID$Group, p = 0.6, list = FALSE)
  train.ID <- data.ID$ID[train.ID]
  data.train <- data %>% filter(ID %in% train.ID) %>% dplyr::select(!ID)
  data.test <- data %>% filter(!ID %in% train.ID) %>% dplyr::select(!ID)
  # 去掉全为0的列：
  data.train <- data.train %>% select_if(is.numeric) %>% select_if(~sum(.) > 0)
  data.train$Group <- ifelse(str_detect(rownames(data.train), 'D'),'AD','Con') %>% factor(levels = c('AD','Con'))
  data.train <- relocate(data.train, Group)
  fn2 <- str_c('Results/HumanUrine/subCluster/caret/RDS/pseudoPer10.Res0.1.Cluster',j,'_Train.2022_Logistic.csv')
  dt <- read.csv(fn2, row.names = 1)
  myPro <- rownames(dt)[1:5]
  roc.data <- data.train %>% select(Group, all_of(myPro))
  roc.data$Group <- roc.data$Group %>% factor(levels = c('Con','AD'))
  roc.data <- relocate(roc.data, Group)
  f2 <- glm(Group ~ roc.data[,2] + roc.data[,3], data = roc.data, family = binomial())
  roc.data$Pro2 <- predict(f2, newdata = roc.data, type = "response")
  f3 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4], data = roc.data, family = binomial())
  roc.data$Pro3 <- predict(f3, newdata = roc.data, type = "response")
  f4 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4] + roc.data[,5], data = roc.data, family = binomial())
  roc.data$Pro4 <- predict(f4, newdata = roc.data, type = "response")
  f5 <- glm(Group ~ roc.data[,2] + roc.data[,3] + roc.data[,4] + roc.data[,5] + roc.data[,6], data = roc.data, family = binomial())
  roc.data$Pro5 <- predict(f5, newdata = roc.data, type = "response")
  auc <- c()
  z <- c()
  pvalue <- c()
  rocPlot1 <- roc(roc.data[,1] ~ roc.data[,2], levels = c('Con','AD'), direction = '<', data = roc.data)
  auc[1] <- rocPlot1$auc
  se <- sqrt(pROC::var(rocPlot1))
  b <- auc[1] - 0.5
  z[1] <- (b / se)
  pvalue[1] <- 2 * pt(-abs(z[1]), df=Inf)
  for(a in 2:5){
    rocPlot <- roc(Group ~ roc.data[, a + 5], levels = c('Con','AD'), direction = '<', data = roc.data)
    auc[a] <- rocPlot$auc
    se <- sqrt(pROC::var(rocPlot))
    b <- auc[a] - 0.5
    z[a] <- (b / se)
    pvalue[a] <- 2 * pt(-abs(z[a]), df=Inf)
  }
  auc_pvalue <- data.frame(auc, z, pvalue)
  rownames(auc_pvalue) <- str_c('Pro', 1:5)
  fn3 <- str_c('Results/Figure6/Files/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster',j,'.Train_2022.ROC.csv')
  write.csv(auc_pvalue, file = fn3)
}
# -- ------------------------------------------------------------------------------------------------------------------
# 将各个亚群的Train ROC合并成一个总的csv
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
filenames <- Sys.glob('Results/Figure6/Files/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster*.Train_2022.ROC.csv')
clusters <- filenames %>% str_remove('Results/Figure6/Files/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster') %>%
  str_remove('.Train_2022.ROC.csv')
data <- read.csv(filenames[1])
data$Cluster <- clusters[1]
for (i in 2:length(clusters)) {
  temp <- read.csv(filenames[i])
  if (nrow(temp) > 0) {
    temp$Cluster <- clusters[i]
    data <- rbind(data, temp)
  }
}
write_csv(data, file = 'Results/Figure6/HumanUrine.pseudoPer10_SCT.Harmony.Dim1-20.0.1.subClusters_Train_2022.ROC.csv')
# -- ------------------------------------------------------------------------------------------------------------------
# 取3次random的结果对combine区分效果进行图示
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data1 <- read.csv('Results/Figure6/HumanUrine.pseudoPer10_SCT.Harmony.Dim1-20.0.1.subClusters_Train_5678.ROC.csv')
data1$ID <- str_c(data1$Cluster, data1$X, sep = '_')
data2 <- read.csv('Results/Figure6/HumanUrine.pseudoPer10_SCT.Harmony.Dim1-20.0.1.subClusters_Train_1234.ROC.csv')
data2$ID <- str_c(data2$Cluster, data2$X, sep = '_')
data3 <- read.csv('Results/Figure6/HumanUrine.pseudoPer10_SCT.Harmony.Dim1-20.0.1.subClusters_Train_2020.ROC.csv')
data3$ID <- str_c(data3$Cluster, data3$X, sep = '_')
data <- full_join(data1, data2, by = 'ID')
data <- full_join(data, data3, by = 'ID')
data <- data %>% dplyr::select(ID, Cluster, X, auc, auc.x, auc.y)
data$Mean <- rowMeans(data[4:6])
data$Mean2 <- ifelse(data$Mean > 0.5, data$Mean, 1-data$Mean)
range(data$Mean2)
data <- data %>% arrange(Cluster, X)
# 作图
show_col(pal_igv(alpha = 0.8)(56))
# ggplot2绘图点的形状不够用怎么办？
shape_level <- nlevels(factor(data$Cluster))
if (shape_level < 15){
  shapes = (0:shape_level) %% 15
} else{
  shapes = c(0:14, c((15:shape_level) %% 110 + 18))
} 
p <- ggplot(data[1:50,], mapping = aes(x = X, y = Mean2)) + theme_bw() + 
  geom_point(aes(colour = factor(Cluster), shape = factor(Cluster)), size = 2) + 
  scale_shape_manual(values = shapes) +
  geom_line(aes(colour = factor(Cluster), group = factor(Cluster)), size = 0.5) + 
  scale_color_igv() + 
  xlab('Protein Number') + ylab('Mean AUC of 3 Random') + 
  theme(legend.position="none")
ggsave(p, filename = 'Results/Figure6/HumanUrine.pseudoPer10.20-0.1.subClusters.10_Train.ROC.AUC.png', 
       width = 3.6, height = 3.6)
ggsave(p, filename = 'Results/Figure6/HumanUrine.pseudoPer10.20-0.1.subClusters.10_Train.ROC.AUC.pdf', 
       width = 3.6, height = 3.6)
p <- ggplot(data, mapping = aes(x = X, y = Mean2)) + theme_bw() + 
  geom_point(aes(colour = factor(Cluster), shape = factor(Cluster)), size = 2) + 
  scale_shape_manual(values = shapes) +
  geom_line(aes(colour = factor(Cluster), group = factor(Cluster)), size = 0.5) + 
  scale_color_igv() + 
  xlab('Protein Number') + ylab('Mean AUC of 3 Random') + 
  theme(legend.position="none")
ggsave(p, filename = 'Results/Figure6/HumanUrine.pseudoPer10.20-0.1.subClusters_Train.ROC.AUC.png', 
       width = 3.6, height = 3.6)
ggsave(p, filename = 'Results/Figure6/HumanUrine.pseudoPer10.20-0.1.subClusters_Train.ROC.AUC.pdf', 
       width = 3.6, height = 3.6)


# -- ------------------------------------------------------------------------------------------------------------------
# 11.将各个亚群的caret ROC合并成一个总的csv，求Mean
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data1 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.5678.caret.ROC.csv')
data1$ID <- str_c(data1$Cluster, data1$method, sep = '_')
data2 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.2020.caret.ROC.csv')
data2$ID <- str_c(data2$Cluster, data2$method, sep = '_')
data3 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.1234.caret.ROC.csv')
data3$ID <- str_c(data3$Cluster, data3$method, sep = '_')
# data4 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.2021.caret.ROC.csv')
# data4$ID <- str_c(data4$Cluster, data4$method, sep = '_')
# data5 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.2022.caret.ROC.csv')
# data5$ID <- str_c(data5$Cluster, data5$method, sep = '_')
data <- full_join(data1, data2, by = 'ID')
data <- full_join(data, data3, by = 'ID')
# data <- full_join(data, data4, by = 'ID')
# data <- full_join(data, data5, by = 'ID')
# data <- data %>% dplyr::select(ID, Cluster, method, Train, Train.x, Train.y, Train.x.x, Train.y.y, Test, Test.x, Test.y, Test.x.x, Test.y.y)
# -- ------------------------------------------------------------------------------------------------------------------
# 取几个数中最大的3个
# df <- apply(data[9:13], 1, function(x){sort(x, decreasing=T)[1:3]}) %>% t() %>% as.data.frame()
# names(df) <- str_c('Test', 1:3, sep = '.')
# df$Test.Mean <- rowMeans(df)
# dt <- cbind(data, df)
# dt <- dt %>% arrange(desc(Test.Mean))
# -- ------------------------------------------------------------------------------------------------------------------
data <- data %>% dplyr::select(ID, Cluster, method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data$Train.Mean <- rowMeans(data[4:6])
data$Train.sd <- apply(data[4:6], 1, sd)
data$Test.Mean <- rowMeans(data[7:9])
data$Test.sd <- apply(data[7:9], 1, sd)
data <- data %>% arrange(desc(Test.Mean), Cluster, method)
data.Urine <- filter(data, Test.Mean > 0.75)
rownames(data.Urine) <- data.Urine$ID
data.Urine.heatmap <- data.Urine %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)



data1 <- read.csv('Results/Figure6/Files/HumanSaliva.AllSamples.pseudoPer20.Res0.1.subClusters_Train.2021.caret.ROC.csv')
data1$ID <- str_c(data1$Cluster, data1$method, sep = '_')
data2 <- read.csv('Results/Figure6/Files/HumanSaliva.AllSamples.pseudoPer20.Res0.1.subClusters_Train.2022.caret.ROC.csv')
data2$ID <- str_c(data2$Cluster, data2$method, sep = '_')
data3 <- read.csv('Results/Figure6/Files/HumanSaliva.AllSamples.pseudoPer20.Res0.1.subClusters_Train.5678.caret.ROC.csv')
data3$ID <- str_c(data3$Cluster, data3$method, sep = '_')
data <- full_join(data1, data2, by = 'ID')
data <- full_join(data, data3, by = 'ID')
data <- data %>% dplyr::select(ID, Cluster, method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data$Train.Mean <- rowMeans(data[4:6])
data$Train.sd <- apply(data[4:6], 1, sd)
data$Test.Mean <- rowMeans(data[7:9])
data$Test.sd <- apply(data[7:9], 1, sd)
data <- data %>% arrange(desc(Test.Mean), Cluster, method)
data.Saliva <- filter(data, Test.Mean > 0.6)
rownames(data.Saliva) <- data.Saliva$ID
data.Saliva.heatmap <- data.Saliva %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)


data1 <- read.csv('Results/HumanPlasma/subCluster/caret/AllSamples.pseudoPer30.Res0.3.subClusters_Train.2020.caret.ROC.csv')
data1$ID <- str_c(data1$Cluster, data1$method, sep = '_')
data2 <- read.csv('Results/HumanPlasma/subCluster/caret/AllSamples.pseudoPer30.Res0.3.subClusters_Train.2022.caret.ROC.csv')
data2$ID <- str_c(data2$Cluster, data2$method, sep = '_')
data3 <- read.csv('Results/HumanPlasma/subCluster/caret/AllSamples.pseudoPer30.Res0.3.subClusters_Train.5678.caret.ROC.csv')
data3$ID <- str_c(data3$Cluster, data3$method, sep = '_')
data <- full_join(data1, data2, by = 'ID')
data <- full_join(data, data3, by = 'ID')
data <- data %>% dplyr::select(ID, Cluster, method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data$Train.Mean <- rowMeans(data[4:6])
data$Train.sd <- apply(data[4:6], 1, sd)
data$Test.Mean <- rowMeans(data[7:9])
data$Test.sd <- apply(data[7:9], 1, sd)
data <- data %>% arrange(desc(Test.Mean), Cluster, method)
data.Plasma <- filter(data, Test.Mean > 0.65)
rownames(data.Plasma) <- data.Plasma$ID
data.Plasma.heatmap <- data.Plasma %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)



data1 <- read.csv('Results/HumanTear/subCluster/caret/AllSamples.pseudoPer10.Res0.3.subClusters_Train.12345.caret.ROC.csv')
data1$ID <- str_c(data1$Cluster, data1$method, sep = '_')
data2 <- read.csv('Results/HumanTear/subCluster/caret/AllSamples.pseudoPer10.Res0.3.subClusters_Train.2022.caret.ROC.csv')
data2$ID <- str_c(data2$Cluster, data2$method, sep = '_')
data3 <- read.csv('Results/HumanTear/subCluster/caret/AllSamples.pseudoPer10.Res0.3.subClusters_Train.5678.caret.ROC.csv')
data3$ID <- str_c(data3$Cluster, data3$method, sep = '_')
data <- full_join(data1, data2, by = 'ID')
data <- full_join(data, data3, by = 'ID')
data <- data %>% dplyr::select(ID, Cluster, method, Train, Train.x, Train.y, Test, Test.x, Test.y)
data$Train.Mean <- rowMeans(data[4:6])
data$Train.sd <- apply(data[4:6], 1, sd)
data$Test.Mean <- rowMeans(data[7:9])
data$Test.sd <- apply(data[7:9], 1, sd)
data <- data %>% arrange(desc(Test.Mean), Cluster, method)
data.Tear <- filter(data, Test.Mean > 0.65)
rownames(data.Tear) <- data.Tear$ID
data.Tear.heatmap <- data.Tear %>% dplyr::select(Train.Mean, Test.Mean, Train.sd, Test.sd)


# 作图
data.heatmap <- rbind(data.Urine.heatmap, data.Saliva.heatmap, data.Tear.heatmap, data.Plasma.heatmap)
data.heatmap <- data.heatmap %>% round(., 3)
data.heatmap$Train.Show <- str_c(data.heatmap$Train.Mean, '±', data.heatmap$Train.sd, sep = ' ')
data.heatmap$Test.Show <- str_c(data.heatmap$Test.Mean, '±', data.heatmap$Test.sd, sep = ' ')
data.heatmap[1:2] %>% range()
breakList <- seq(0.6, 1, 0.05)
pheatmap(data.heatmap[1:2], color = topo.colors(9, alpha = 0.6), 
         cellwidth = 40, cellheight = 10,
         # display_numbers = T, 
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(10, 29, 43),
         filename = 'Results/Figure6/Human.subClusters_caret.ROC.png') # 保存，自动调整纸张大小
pheatmap(data.heatmap[1:2], color = topo.colors(9, alpha = 0.6), 
         cellwidth = 40, cellheight = 10,
         # display_numbers = T, 
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(10, 29, 43),
         filename = 'Results/Figure6/Human.subClusters_caret.ROC.pdf') # 保存，自动调整纸张大小
pheatmap(data.heatmap[1:2], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         # display_numbers = T, 
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(10, 29, 43),
         filename = 'Results/Figure6/Human.subClusters_caret.ROC.png') # 保存，自动调整纸张大小
pheatmap(data.heatmap[1:2], color = pal_material(palette = 'red', n = length(breakList))(length(breakList)), 
         cellwidth = 40, cellheight = 10,
         # display_numbers = T, 
         number_color = 'black', number_format = '%.2f', show_rownames = T,
         display_numbers = data.heatmap[5:6], #设置热图区分标记
         fontsize_number = 6, fontsize_row = 8, fontsize_col = 10,
         cluster_rows = F, cluster_cols = F, treeheight_row = 0,
         legend_breaks = breakList,
         gaps_row = c(10, 24, 44),
         filename = 'Results/Figure6/Human.subClusters_caret.ROC.pdf') # 保存，自动调整纸张大小


# -- ------------------------------------------------------------------------------------------------------------------
# 展示最大AUC亚群的AD vs Con的差异
# -- ------------------------------------------------------------------------------------------------------------------
# Urine
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
data1 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.1234_Logistic.Pos.csv')
data2 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.5678_Logistic.Pos.csv')
data3 <- read.csv('Results/HumanUrine/subCluster/caret/AllSamples.pseudoPer10.Res0.1.subClusters_Train.2020_Logistic.Pos.csv')
# Cluster == 3
data1 <- data1 %>% filter(Cluster == 3)
data2 <- data2 %>% filter(Cluster == 3)
data3 <- data3 %>% filter(Cluster == 3)
my_features <- c(data1$X[1:3], data2$X[1:3], data3$X[1:3]) %>% unique()
dt <- read.csv('Results/HumanUrine/subCluster/Expression/HumanUrine.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.1.Cluster3.Expression.SCT.csv', row.names = 1)
dt <- dt %>% t() %>% as.data.frame()
dt.sub <- dt %>% select(all_of(my_features))
dt.sub$Group <- ifelse(str_detect(rownames(dt.sub), 'D'), 'AD', 'Con')
dt.sub <- relocate(dt.sub, Group)
for (i in 2:ncol(dt.sub)) {
  myTitle <- colnames(dt.sub)[i]
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(dt.sub, aes(x = Group, fill = Group)) + theme_bw() +
    # geom_violin(aes_(y = as.name(colnames(dt.sub)[i])), scale = 'width', width = 0.7, alpha = 0.7) +
    geom_boxplot(aes_(y = as.name(colnames(dt.sub)[i])), width = 0.6, size = 0.6) +  
    geom_jitter(aes_(y = as.name(colnames(dt.sub)[i])), alpha = 0.6, size = 2) + 
    geom_signif(aes_(y = as.name(colnames(dt.sub)[i])), comparisons = list(c("AD","Con")), 
                map_signif_level = F, step_increase = 0.05, tip_length = 0.01, 
                test = "t.test", textsize = 2.5) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure6/HumanUrine.pseudoPer10.Dim1-20.0.1.Cluster3.',myTitle, '.Pro.png')
  ggsave(p2, filename = fn2, width = 2.5, height = 3.5)
  fn2 <- str_c('Results/Figure6/HumanUrine.pseudoPer10.Dim1-20.0.1.Cluster3.',myTitle, '.Pro..pdf')
  ggsave(p2, filename = fn2, width = 2.5, height = 3.5)
}


# -- ------------------------------------------------------------------------------------------------------------------
# Saliva
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
# Cluster == 14
data1 <- read.csv('Results/Figure6/HumanSaliva/Logistic/AllSamples.pseudoPer20.Res0.1.Cluster14_Train.2021_Logistic.csv')
data2 <- read.csv('Results/Figure6/HumanSaliva/Logistic/AllSamples.pseudoPer20.Res0.1.Cluster14_Train.2022_Logistic.csv')
data3 <- read.csv('Results/Figure6/HumanSaliva/Logistic/AllSamples.pseudoPer20.Res0.1.Cluster14_Train.5678_Logistic.csv')
my_features <- c(data1$X[1:3], data2$X[1:3], data3$X[1:3]) %>% unique()
dt <- read.csv('Results/Figure6/HumanSaliva/Expression/HumanSaliva.pseudoPer20.Merge_SCT_Harmony.Dim1-20.0.1.Cluster14.Expression.SCT.csv', row.names = 1)
dt <- dt %>% t() %>% as.data.frame()
dt.sub <- dt %>% select(all_of(my_features))
dt.sub$Group <- ifelse(str_detect(rownames(dt.sub), 'D'), 'AD', 'Con')
dt.sub <- relocate(dt.sub, Group)
for (i in 2:ncol(dt.sub)) {
  myTitle <- colnames(dt.sub)[i]
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(dt.sub, aes(x = Group, fill = Group)) + theme_bw() +
    geom_boxplot(aes_(y = as.name(colnames(dt.sub)[i])), width = 0.5, size = 0.6) +  
    geom_jitter(aes_(y = as.name(colnames(dt.sub)[i])), width = 0.2, alpha = 0.6, size = 2) + 
    # geom_violin(aes_(y = as.name(colnames(dt.sub)[i])), scale = 'width', width = 0.7, alpha = 0.7) +
    geom_signif(aes_(y = as.name(colnames(dt.sub)[i])), comparisons = list(c("AD","Con")), 
                map_signif_level = F, step_increase = 0.05, tip_length = 0.01, 
                test = "t.test", textsize = 2.5) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure6/HumanSaliva.pseudoPer20.Dim1-20.0.1.Cluster14.',myTitle, '.Pro.png')
  ggsave(p2, filename = fn2, width = 2.5, height = 3.5)
  fn2 <- str_c('Results/Figure6/HumanSaliva.pseudoPer20.Dim1-20.0.1.Cluster14.',myTitle, '.Pro..pdf')
  ggsave(p2, filename = fn2, width = 2.5, height = 3.5)
}


# -- ------------------------------------------------------------------------------------------------------------------
# Plasma
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
# Cluster == 10
data1 <- read.csv('Results/HumanPlasma/subCluster/caret/RDS/pseudoPer30.Res0.3.Cluster10_Train.1234_Logistic.csv')
data2 <- read.csv('Results/HumanPlasma/subCluster/caret/RDS/pseudoPer30.Res0.3.Cluster10_Train.5678_Logistic.csv')
data3 <- read.csv('Results/HumanPlasma/subCluster/caret/RDS/pseudoPer30.Res0.3.Cluster10_Train.2022_Logistic.csv')
my_features <- c(data1$X[1:3], data2$X[1:3], data3$X[1:3]) %>% unique()
dt <- read.csv('Results/HumanPlasma/subCluster/Expression/HumanPlasma.AllSamples.pseudoPer30.Merge_SCT_Harmony.Dim1-20.0.3.Cluster10.Expression.SCT.csv', row.names = 1)
dt <- dt %>% t() %>% as.data.frame()
dt.sub <- dt %>% select(all_of(my_features))
dt.sub$Group <- ifelse(str_detect(rownames(dt.sub), 'D'), 'AD', 'Con')
dt.sub <- relocate(dt.sub, Group)
for (i in 2:ncol(dt.sub)) {
  myTitle <- colnames(dt.sub)[i]
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(dt.sub, aes(x = Group, fill = Group)) + theme_bw() +
    geom_boxplot(aes_(y = as.name(colnames(dt.sub)[i])), width = 0.5, size = 0.6) +  
    geom_jitter(aes_(y = as.name(colnames(dt.sub)[i])), width = 0.2, alpha = 0.6, size = 2) + 
    # geom_violin(aes_(y = as.name(colnames(dt.sub)[i])), scale = 'width', width = 0.7, alpha = 0.7) +
    geom_signif(aes_(y = as.name(colnames(dt.sub)[i])), comparisons = list(c("AD","Con")), 
                map_signif_level = F, step_increase = 0.05, tip_length = 0.01, 
                test = "t.test", textsize = 2.5) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure6/HumanPlasma.pseudoPer30.Dim1-20.0.3.Cluster10.',myTitle, '.Pro.png')
  ggsave(p2, filename = fn2, width = 2.5, height = 3.5)
  fn2 <- str_c('Results/Figure6/HumanPlasma.pseudoPer30.Dim1-20.0.3.Cluster10.',myTitle, '.Pro..pdf')
  ggsave(p2, filename = fn2, width = 2.5, height = 3.5)
}


# -- ------------------------------------------------------------------------------------------------------------------
# Tear
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
# Cluster == 7
data1 <- read.csv('Results/HumanTear/subCluster/caret/RDS/pseudoPer10.Res0.3.Cluster7_Train.1234_Logistic.csv')
data2 <- read.csv('Results/HumanTear/subCluster/caret/RDS/pseudoPer10.Res0.3.Cluster7_Train.5678_Logistic.csv')
data3 <- read.csv('Results/HumanTear/subCluster/caret/RDS/pseudoPer10.Res0.3.Cluster10_Train.2022_Logistic.csv')
my_features <- c(data1$X[1:3], data2$X[1:3], data3$X[1:3]) %>% unique()
dt <- read.csv('Results/HumanTear/subCluster/Expression/HumanTear.pseudoPer10.Merge_SCT_Harmony.Dim1-20.0.3.Cluster7.Expression.SCT.csv', row.names = 1)
dt <- dt %>% t() %>% as.data.frame()
dt.sub <- dt %>% select(all_of(my_features))
dt.sub$Group <- ifelse(str_detect(rownames(dt.sub), 'D'), 'AD', 'Con')
dt.sub <- relocate(dt.sub, Group)
for (i in 2:ncol(dt.sub)) {
  myTitle <- colnames(dt.sub)[i]
  # 在aes()函数中无法编程，需要用aes_()函数传递表达式的结果
  p2 <- ggplot(dt.sub, aes(x = Group, fill = Group)) + theme_bw() +
    geom_boxplot(aes_(y = as.name(colnames(dt.sub)[i])), width = 0.5, size = 0.6) +  
    geom_jitter(aes_(y = as.name(colnames(dt.sub)[i])), width = 0.2, alpha = 0.6, size = 2) + 
    # geom_violin(aes_(y = as.name(colnames(dt.sub)[i])), scale = 'width', width = 0.7, alpha = 0.7) +
    geom_signif(aes_(y = as.name(colnames(dt.sub)[i])), comparisons = list(c("AD","Con")), 
                map_signif_level = F, step_increase = 0.05, tip_length = 0.01, 
                test = "t.test", textsize = 2.5) +
    xlab('') + ylab('') + ggtitle(myTitle) +
    scale_fill_manual(values = c('#CC0000FF','#6699FFFF')) +
    theme(legend.position = "none")
  fn2 <- str_c('Results/Figure6/HumanTear.pseudoPer10.Dim1-20.0.3.Cluster7.',myTitle, '.Pro.png')
  ggsave(p2, filename = fn2, width = 2.5, height = 3.5)
  fn2 <- str_c('Results/Figure6/HumanTear.pseudoPer10.Dim1-20.0.3.Cluster7.',myTitle, '.Pro..pdf')
  ggsave(p2, filename = fn2, width = 2.5, height = 3.5)
}


