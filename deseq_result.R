## load package
library("readxl")
library("magrittr")
library(dplyr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(stats)
library(WGCNA)
library(RColorBrewer)
library(pheatmap)
library(genefilter) 

## read and process data
df <- read_excel("1-s2.0-S1074761317300729-mmc2.xlsx", skip=2) %>% data.frame() %>% 
  select("Gene.Names", "WPC", "WT", "...9", "...10", "...11", "...12", "...13", "...14","...15") %>%
  distinct(Gene.Names, .keep_all = TRUE)  # remove replicated row

# change col and row names
col_names <- c("Gene_id", "WPC", t(df[1, 3:ncol(df)]))
colnames(df) <- col_names
df <- filter(df, !is.na(df$Gene_id))

# data preparation
gene_id <- df$Gene_id # gene id
expr_df <- df[, 3:ncol(df)] %>% mutate_all(function(x) as.numeric(as.character(x)))

## Deseq2
coldata <- data.frame(condition = factor(col_names[-c(1, 2)]), levels = factor(rep(c('0h',  '2h', '8h', '16h'), each=2)))
# get dds matrix
dds <- DESeqDataSetFromMatrix(countData = round(expr_df), colData = coldata, design= ~levels)
# DESeq
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
vsd <- varianceStabilizingTransformation(dds1)

wpn_vsd <- getVarianceStabilizedData(dds1)

# get result
res <- results(dds1, contrast = c('levels', '16h', '0h'))
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

# output result
res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.05),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.05),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.05),'sig'] <- 'none'

# overall differential expression gene
res1_select <- subset(res1, sig %in% c('up', 'down'))
up <- res1_select[res1_select$sig == 'up',]
down <- res1_select[res1_select$sig == 'down',]

# volcano plot
p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 1) +  # scatter plot
  scale_color_manual(values = c('red', 'gray', 'green'), limits = c('up', 'none', 'down')) +  #colour
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'control vs treat', color = '') +  # axis
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #background colour, grid line, title
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.3, 0.3), lty = 3, color = 'black') +  # setting threshold
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(-12, 12) + ylim(0, 35)  # define boundary
p

# expression matrix 
# up
up_id <- as.numeric(rownames(up))
up_df <- wpn_vsd[up_id, ]
pheatmap(up_df, scale="row", cluster_rows = FALSE, cluster_cols = FALSE, labels_row=NA)
# down
down_id <- as.numeric(rownames(down))
down_df <- wpn_vsd[down_id, ]
pheatmap(down_df, scale="row", cluster_rows = FALSE, cluster_cols = FALSE, labels_row=NA)
# all
de_id <- sort(c(up_id, down_id))
de_df <- wpn_vsd[de_id, ]
pheatmap(de_df, scale="row", cluster_rows = FALSE, cluster_cols = FALSE, labels_row=NA)

## WGCNA
input_mat = t(de_df)

allowWGCNAThreads() # allow muti thread

# Choose a set of soft-thresholding powers
powers = c(1:10, seq(11, 30, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

# visual scale free model fit
par(mfrow = c(1,2))
k = softConnectivity(datExpr = input_mat, power=1)
hist(k)
scaleFreePlot(k, main="Check scale free topology")

# from the plot we choose 70 as ß here
picked_power = 1
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

mergedColors = labels2colors(netwk$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

wpc <- factor(netwk$colors)

## pheatmap
wpc_df <- data.frame(WPC = as.numeric(wpc), row.names=de_id)
pheatmap(de_df, scale="row", cluster_rows = FALSE, cluster_cols = FALSE, labels_row=NA, annotation_row = wpc_df)

## Step3: PPI
de_gene <- data.frame(gene_id = gene_id[de_id])

wpc1_id <- data.frame(gene_id = de_gene[wpc_df$WPC==1, ])
wpc2_id <- data.frame(gene_id = de_gene[wpc_df$WPC==2, ])
wpc3_id <- data.frame(gene_id = de_gene[wpc_df$WPC==3, ])

# wpc1
string_db <- STRINGdb$new( version="11.5", 
                           species=10090,   # mice 10090 
                           input_directory="") 

dat_map1 <- string_db$map(my_data_frame=wpc1_id, 
                          my_data_frame_id_col_names="gene_id", 
                          removeUnmappedRows = TRUE )
hits1 <- dat_map1$STRING_id 
png("deseq2_PPI_wpc1.png",units="in",width = 10, height = 10, res=400)
string_db$plot_network(hits1)
dev.off()
# wpc2
dat_map2 <- string_db$map(my_data_frame=wpc2_id, 
                          my_data_frame_id_col_names="gene_id", 
                          removeUnmappedRows = TRUE )
hits2 <- dat_map2$STRING_id 

png("deseq2_PPI_wpc2.png",units="in",width = 10, height = 10, res=400)
string_db$plot_network(hits2)

# wpc3
dat_map3 <- string_db$map(my_data_frame=wpc3_id, 
                          my_data_frame_id_col_names="gene_id", 
                          removeUnmappedRows = TRUE )
hits3 <- dat_map3$STRING_id 
png("deseq2_PPI_wpc3.png",units="in",width = 10, height = 10, res=400)
string_db$plot_network(hits3)
dev.off()

