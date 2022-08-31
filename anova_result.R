setwd("/Users/scz/Documents/Intership/Anova sample")
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
library(STRINGdb) 
library(igraph)

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

## differential expression
extract_p <- function(df) {   # function
  data1 <- data.frame("Time" = as.factor(rep(c(0, 2, 8, 16), each=2)), "exp" = c(t(df)))
  
  one.way <- aov(exp ~ Time, data = data1)
  pvalue <- summary(one.way)[[1]][["Pr(>F)"]][1]
  
  return(pvalue)
}

## Step1: Differential expression
# pvalue and qvalue
pvalue <- apply(expr_df, 1, extract_p) 
names(pvalue) <- 1:length(pvalue)
qvalue <- p.adjust(p = pvalue, method = 'fdr')

df_list <- qvalue < 0.005
de_df <- expr_df[df_list, ] # de matrix
de_id<- rownames(de_df)        # de gene id 

# jaccard index between paper result and my anova result
selected_id <- seq(1, 1608, 1)
intersect <- intersect(de_id, selected_id)
union <- union(de_id, selected_id)

iou <- length(intersect) / length(union)
iou


## Step2: WGCNA
input_mat = t(de_df)

allowWGCNAThreads() # allow muti thread

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

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

# from the plot we choose 6 as ß here
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
wpc_df <- data.frame(WPC = as.numeric(wpc))
pheatmap(de_df, scale="row", cluster_rows = FALSE, cluster_cols = FALSE, labels_row=NA, annotation_row = wpc_df)

## Step3: PPI
de_gene <- data.frame(gene_id = gene_id[df_list])

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
png("anova_PPI_wpc1.png",units="in",width = 10, height = 10, res=400)
string_db$plot_network(hits1)

# wpc2
dat_map2 <- string_db$map(my_data_frame=wpc2_id, 
                          my_data_frame_id_col_names="gene_id", 
                          removeUnmappedRows = TRUE )
hits2 <- dat_map2$STRING_id 
png("anova_PPI_wpc2.png",units="in",width = 10, height = 10, res=400)
string_db$plot_network(hits2)

# wpc3
dat_map3 <- string_db$map(my_data_frame=wpc3_id, 
                          my_data_frame_id_col_names="gene_id", 
                          removeUnmappedRows = TRUE )
hits3 <- dat_map3$STRING_id 
png("anova_PPI_wpc3.png",units="in",width = 10, height = 10, res=400)
string_db$plot_network(hits3)
dev.off()
