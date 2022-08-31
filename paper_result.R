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
df <- df[2:nrow(df), ] 

# data preparation
gene_id <- df$Gene_id # gene id
expr_df <- df[, 3:ncol(df)] %>% mutate_all(function(x) as.numeric(as.character(x)))

## differential expression
de_df <- expr_df[1:1608, ]
wpc <- df$WPC[1:1608]

## pheatmap
wpc_df <- as.data.frame(as.numeric(as.factor(wpc)))
colnames(wpc_df) <- "WPC"

pheatmap(de_df, scale="row", cluster_rows = FALSE, cluster_cols = FALSE, labels_row=NA, annotation_row = wpc_df)

## Step3: PPI 
de_gene <- data.frame(gene_id = gene_id[1:1608])

wpc1_id <- data.frame(gene_id = de_gene[wpc_df$WPC==1, ])
wpc2_id <- data.frame(gene_id = de_gene[wpc_df$WPC==2, ])
wpc3_id <- data.frame(gene_id = de_gene[wpc_df$WPC==3, ])
wpc4_id <- data.frame(gene_id = de_gene[wpc_df$WPC==4, ])
wpc5_id <- data.frame(gene_id = de_gene[wpc_df$WPC==5, ])
wpc6_id <- data.frame(gene_id = de_gene[wpc_df$WPC==6, ])

# wpc1
string_db <- STRINGdb$new( version="11.5", 
                           species=10090,   # mice 10090 
                           input_directory="") 

dat_map1 <- string_db$map(my_data_frame=wpc1_id, 
                          my_data_frame_id_col_names="gene_id", 
                          removeUnmappedRows = TRUE )
hits1 <- dat_map1$STRING_id 
png("PPI_wpc1.png",units="in",width = 10, height = 10, res=400)
string_db$plot_network(hits1)
dev.off()

# wpc2
dat_map2 <- string_db$map(my_data_frame=wpc2_id, 
                          my_data_frame_id_col_names="gene_id", 
                          removeUnmappedRows = TRUE )
hits2 <- dat_map2$STRING_id 
png("PPI_wpc2.png",units="in",width = 10, height = 10, res=400)
string_db$plot_network(hits2)
dev.off()

# wpc3
dat_map3 <- string_db$map(my_data_frame=wpc3_id, 
                          my_data_frame_id_col_names="gene_id", 
                          removeUnmappedRows = TRUE )
hits3 <- dat_map3$STRING_id 
png("PPI_wpc3.png",units="in",width = 10, height = 10, res=400)
string_db$plot_network(hits3)
dev.off()

# wpc4
dat_map4 <- string_db$map(my_data_frame=wpc4_id, 
                          my_data_frame_id_col_names="gene_id", 
                          removeUnmappedRows = TRUE )
hits4 <- dat_map4$STRING_id 
png("PPI_wpc4.png",units="in",width = 10, height = 10, res=400)
string_db$plot_network(hits4)
dev.off()

# wpc5
dat_map5 <- string_db$map(my_data_frame=wpc5_id, 
                          my_data_frame_id_col_names="gene_id", 
                          removeUnmappedRows = TRUE )
hits5 <- dat_map5$STRING_id 
png("PPI_wpc5.png",units="in",width = 10, height = 10, res=400)
string_db$plot_network(hits5)
dev.off()

# wpc6
dat_map6 <- string_db$map(my_data_frame=wpc6_id, 
                          my_data_frame_id_col_names="gene_id", 
                          removeUnmappedRows = TRUE )
hits6 <- dat_map5$STRING_id 
png("PPI_wpc6.png",units="in",width = 10, height = 10, res=400)
string_db$plot_network(hits6)
dev.off()





