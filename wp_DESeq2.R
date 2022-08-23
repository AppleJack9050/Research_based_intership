#setwd("D:/Intership")

## Loading
library("readxl")
library("magrittr")
library(dplyr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(stats)

## read data
my_data <- read_excel("1-s2.0-S1074761317300729-mmc2.xlsx", skip=2) %>% data.frame() %>% select("Gene.Names", "WT", "...9","...14","...15")
# change colnames
colnames(my_data) <- c("Gene", "0ha", "0hb", "16ha", "16hb") 
# remove line one and get gene id
my_data <- my_data[-1,]
gene_id <- my_data[, 1]
my_data <- my_data[, -1]

## Differential Expression
coldata <- data.frame(condition = factor(rep(c('0h', '16h'), each = 2), levels = c('0h', '16h')))
new_data <- data.frame(apply(my_data, 2, function(x) as.numeric(as.character(x))))
# get dds matrix
dds <- DESeqDataSetFromMatrix(countData = round(new_data), colData = coldata, design= ~condition)
# DESeq
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
# get result
res <- results(dds1, contrast = c('condition', '0h', '16h'))


## Output result
# write to local
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res1, 'control_treat.DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)

# gene filter
res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

res1[which(res1$log2FoldChange >= 0.5 & res1$padj < 0.05),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -0.5 & res1$padj < 0.05),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 0.5 | res1$padj >= 0.05),'sig'] <- 'none'

# overall differential expression gene
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = 'control_treat.DESeq2.select.txt', sep = '\t', col.names = NA, quote = FALSE)

# output up and down
res1_up <- subset(res1, sig == 'up')
res1_down <- subset(res1, sig == 'down')

write.table(res1_up, file = 'control_treat.DESeq2.up.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(res1_down, file = 'control_treat.DESeq2.down.txt', sep = '\t', col.names = NA, quote = FALSE)


## plot
p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 1) +  # scatter plot
  scale_color_manual(values = c('red', 'gray', 'green'), limits = c('up', 'none', 'down')) +  #colour
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'control vs treat', color = '') +  # axis
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #background colour, grid line, title
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  # setting threshold
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(-12, 12) + ylim(0, 35)  # define boundary
p


## Compare with original data set 
df <- read_excel("1-s2.0-S1074761317300729-mmc2.xlsx", skip=2) %>% 
  data.frame() %>% 
  select("Gene.Names", "ANOVA.P.values", "BH.FDR") %>% 
  drop_na()

# jaccard index within selected protein and anova_id
selected_id <- seq(1, 1704, 1)

deseq2_id <- rownames(res1_select)

intersect1 <- intersect(deseq2_id, selected_id)
union1 <- union(deseq2_id, selected_id)

iou1 <- length(intersect1) / length(union1)
iou1

# jaccard index within selected protein and anova_id
anova_id <- rownames(subset(df, BH.FDR < 0.05))
union2 <- union(deseq2_id, anova_id)

intersect2 <- intersect(deseq2_id, anova_id)
iou2 <- length(intersect2) / length(union2)
iou2

## Q-value
q_value <- p.adjust(p = df$ANOVA.P.values, method = 'fdr')

