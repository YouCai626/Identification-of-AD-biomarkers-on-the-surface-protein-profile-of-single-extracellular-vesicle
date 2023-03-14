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
# 可视化常用包
library(tidyverse) # 包含了ggplot2等包
library(ggsci) # 期刊配色，单个离散配色模块中最多只有51种颜色
library(scales) # show_col
library(patchwork) # 拼接图
library(ggbreak) # 坐标轴截断
library(ggsignif) # 计算差异显著性并添加显著性线段
library(ggthemes) # 主题设置
library(stringr) # 文本处理包，包含在tidyverse中

# -- ------------------------------------------------------------------------------------------------------------------
# FigS5.数据的频率分布图，按照百分比展示(每个体液种类的展示几个样本)
# -- ------------------------------------------------------------------------------------------------------------------
# Mouse的三种体液
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
sampleID <- read.csv('Mouse_SampleID_SampleName EV counts.csv')
head(sampleID)
i = 41
fn1 <- str_c('RawData/Mouse/',sampleID[i,1],'.total_ev_protein.csv')
data <- read.csv(fn1, header = T)
wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
wide_data[1:10,1:10]
wide.data <- wide_data[,-c(1:2)]
wide.data <- ifelse(wide.data > 0, 1, 0)
wide.data[1:10,1:10]
sum.data <- rowSums(wide.data)
proKind.Serum <- table(sum.data) %>% as.data.frame()
names(proKind.Serum) <- c('ProKind', 'Frequency')
str(proKind.Serum)
proKind.Serum$Serum <- proKind.Serum$Frequency/sum(proKind.Serum$Frequency)*100

i = 65
fn1 <- str_c('RawData/Mouse/',sampleID[i,1],'.total_ev_protein.csv')
data <- read.csv(fn1, header = T)
wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
wide_data[1:10,1:10]
wide.data <- wide_data[,-c(1:2)]
wide.data <- ifelse(wide.data > 0, 1, 0)
wide.data[1:10,1:10]
sum.data <- rowSums(wide.data)
proKind.Urine <- table(sum.data) %>% as.data.frame()
names(proKind.Urine) <- c('ProKind', 'Frequency')
str(proKind.Urine)
proKind.Urine$Urine <- proKind.Urine$Frequency/sum(proKind.Urine$Frequency)*100

i = 17
fn1 <- str_c('RawData/Mouse/',sampleID[i,1],'.total_ev_protein.csv')
data <- read.csv(fn1, header = T)
wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
wide_data[1:10,1:10]
wide.data <- wide_data[,-c(1:2)]
wide.data <- ifelse(wide.data > 0, 1, 0)
wide.data[1:10,1:10]
sum.data <- rowSums(wide.data)
proKind.NIS <- table(sum.data) %>% as.data.frame()
names(proKind.NIS) <- c('ProKind', 'Frequency')
str(proKind.NIS)
proKind.NIS$NIS <- proKind.NIS$Frequency/sum(proKind.NIS$Frequency)*100

proKind <- inner_join(proKind.Serum, proKind.Urine, by = 'ProKind')
proKind <- inner_join(proKind, proKind.NIS, by = 'ProKind')
proKind <- proKind %>% dplyr::select(ProKind, Serum, Urine, NIS) %>% filter(ProKind %in% c(1:10))

data.plot <- pivot_longer(proKind, cols = Serum:NIS, names_to = 'Type', values_to = 'Percentage')
data.plot$Type <- data.plot$Type %>% factor(levels = c('Serum', 'Urine', 'NIS'))

p <- ggplot(data.plot, aes(x = ProKind, y = Percentage, color = Type, shape = Type)) + 
  geom_point(size = 3) + theme_bw() +
  xlab('Protein Kind') + ylab('Percentage') + ggtitle('Percentage of Protein Kind') +
  theme(axis.text.x = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.title = element_text(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 10),
        legend.title = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_y_break(c(16, 77), scales = 0.3)
ggsave(p, filename = 'Results/FigureS5/Mouse_Protein Kind Percentage.png', width = 6, height = 8)
ggsave(p, filename = 'Results/FigureS5/Mouse_Protein Kind Percentage.pdf', width = 6, height = 8)
p <- p + stat_summary(fun = mean, colour = "red", geom = "point", shape = '*', size = 6, aes(group = 1)) + 
  stat_summary(fun = mean, colour = "black", geom = "line", size = 0.8, linetype = 'dashed', aes(group = 1))
ggsave(p, filename = 'Results/FigureS5/Mouse_Protein Kind Percentage.2.png', width = 6, height = 8)
ggsave(p, filename = 'Results/FigureS5/Mouse_Protein Kind Percentage.2.pdf', width = 6, height = 8)

p <- ggplot(data.plot, aes(x = ProKind, y = Percentage, color = Type, shape = Type)) + 
  geom_point(size = 3) + theme_bw() +
  xlab('Protein Kind') + ylab('Percentage') + ggtitle('Percentage of Protein Kind') +
  theme(axis.text.x = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.title = element_text(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 10),
        legend.title = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16))
p <- p + stat_summary(fun = mean, colour = "red", geom = "point", shape = '*', size = 8, aes(group = 1)) + 
  stat_summary(fun = mean, colour = "black", geom = "line", size = 0.8, linetype = 'dashed', aes(group = 1))
ggsave(p, filename = 'Results/FigureS5/Mouse_Protein Kind Percentage.3.png', width = 6, height = 8)
ggsave(p, filename = 'Results/FigureS5/Mouse_Protein Kind Percentage.3.pdf', width = 6, height = 8)



# -- ------------------------------------------------------------------------------------------------------------------
# Human的4种体液
# -- ------------------------------------------------------------------------------------------------------------------
rm(list=ls())
gc()
sampleID <- read.csv('Human_SampleID_SampleName.csv')
head(sampleID)
i = 9
fn1 <- str_c('RawData/',sampleID[i,1],'.total_ev_protein.csv')
data <- read.csv(fn1, header = T)
wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
wide_data[1:10,1:10]
wide.data <- wide_data[,-c(1:2)]
wide.data <- ifelse(wide.data > 0, 1, 0)
wide.data[1:10,1:10]
sum.data <- rowSums(wide.data)
proKind.Plasma <- table(sum.data) %>% as.data.frame()
names(proKind.Plasma) <- c('ProKind', 'Frequency')
str(proKind.Plasma)
proKind.Plasma$Plasma <- proKind.Plasma$Frequency/sum(proKind.Plasma$Frequency)*100

i = 120
fn1 <- str_c('RawData/',sampleID[i,1],'.total_ev_protein.csv')
data <- read.csv(fn1, header = T)
wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
wide_data[1:10,1:10]
wide.data <- wide_data[,-c(1:2)]
wide.data <- ifelse(wide.data > 0, 1, 0)
wide.data[1:10,1:10]
sum.data <- rowSums(wide.data)
proKind.Urine <- table(sum.data) %>% as.data.frame()
names(proKind.Urine) <- c('ProKind', 'Frequency')
str(proKind.Urine)
proKind.Urine$Urine <- proKind.Urine$Frequency/sum(proKind.Urine$Frequency)*100

i = 47
fn1 <- str_c('RawData/',sampleID[i,1],'.total_ev_protein.csv')
data <- read.csv(fn1, header = T)
wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
wide_data[1:10,1:10]
wide.data <- wide_data[,-c(1:2)]
wide.data <- ifelse(wide.data > 0, 1, 0)
wide.data[1:10,1:10]
sum.data <- rowSums(wide.data)
proKind.Saliva <- table(sum.data) %>% as.data.frame()
names(proKind.Saliva) <- c('ProKind', 'Frequency')
str(proKind.Saliva)
proKind.Saliva$Saliva <- proKind.Saliva$Frequency/sum(proKind.Saliva$Frequency)*100

i = 84
fn1 <- str_c('RawData/',sampleID[i,1],'.total_ev_protein.csv')
data <- read.csv(fn1, header = T)
wide_data <- spread(data, key = 'variable', value = 'value', fill = 0, convert = T)
wide_data[1:10,1:10]
wide.data <- wide_data[,-c(1:2)]
wide.data <- ifelse(wide.data > 0, 1, 0)
wide.data[1:10,1:10]
sum.data <- rowSums(wide.data)
proKind.Tear <- table(sum.data) %>% as.data.frame()
names(proKind.Tear) <- c('ProKind', 'Frequency')
str(proKind.Tear)
proKind.Tear$Tear <- proKind.Tear$Frequency/sum(proKind.Tear$Frequency)*100

proKind <- inner_join(proKind.Plasma, proKind.Urine, by = 'ProKind')
proKind <- inner_join(proKind, proKind.Saliva, by = 'ProKind')
proKind <- inner_join(proKind, proKind.Tear, by = 'ProKind')
proKind <- proKind %>% dplyr::select(ProKind, Plasma, Urine, Saliva, Tear) %>% filter(ProKind %in% c(1:10))

data.plot <- pivot_longer(proKind, cols = Plasma:Tear, names_to = 'Type', values_to = 'Percentage')
data.plot$Type <- data.plot$Type %>% factor(levels = c('Plasma', 'Urine', 'Saliva', 'Tear'))

p <- ggplot(data.plot, aes(x = ProKind, y = Percentage, color = Type, shape = Type)) + 
  geom_point(size = 3) + theme_bw() +
  xlab('Protein Kind') + ylab('Percentage') + ggtitle('Percentage of Protein Kind') +
  theme(axis.text.x = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.title = element_text(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 10),
        legend.title = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16)) +
  scale_y_break(c(22, 60), scales = 0.3)
ggsave(p, filename = 'Results/FigureS5/Human_Protein Kind Percentage.png', width = 6, height = 8)
ggsave(p, filename = 'Results/FigureS5/Human_Protein Kind Percentage.pdf', width = 6, height = 8)
p <- p + stat_summary(fun = mean, colour = "red", geom = "point", shape = '*', size = 10, aes(group = 1)) + 
  stat_summary(fun = mean, colour = "black", geom = "line", size = 0.8, linetype = 'dashed', aes(group = 1))
ggsave(p, filename = 'Results/FigureS5/Human_Protein Kind Percentage.2.png', width = 6, height = 8)
ggsave(p, filename = 'Results/FigureS5/Human_Protein Kind Percentage.2.pdf', width = 6, height = 8)

p <- ggplot(data.plot, aes(x = ProKind, y = Percentage, color = Type, shape = Type)) + 
  geom_point(size = 3) + theme_bw() +
  xlab('Protein Kind') + ylab('Percentage') + ggtitle('Percentage of Protein Kind') +
  theme(axis.text.x = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.title = element_text(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 10),
        legend.title = element_text(colour = 'black', size = 12),
        title = element_text(colour = 'black', size = 16))
p <- p + stat_summary(fun = mean, colour = "red", geom = "point", shape = '*', size = 10, aes(group = 1)) + 
  stat_summary(fun = mean, colour = "black", geom = "line", size = 0.8, linetype = 'dashed', aes(group = 1))
ggsave(p, filename = 'Results/FigureS5/Human_Protein Kind Percentage.3.png', width = 6, height = 8)
ggsave(p, filename = 'Results/FigureS5/Human_Protein Kind Percentage.3.pdf', width = 6, height = 8)



