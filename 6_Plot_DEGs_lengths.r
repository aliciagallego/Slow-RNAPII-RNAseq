#!/usr/bin/Rscript

# This script includes codes for the following functions:
# 1. Retrieve gene length from featureCounts output
# 2. Retrieve the list of DEGs computed by DESeq2
# 3. Assign gene lengths to DEGs
# 4. Plot gene length distributions in Up, Down and n.s. DEGs

# ----------------
# Load libraries |
# ----------------
library("ggplot2")

# -------------------------------
# featureCounts output to RPKMs |
# -------------------------------
# Open featureCounts output
COUNT_MAT=read.table('/featureCounts_path/featureCounts_file.txt',row.names=1,header=TRUE,sep='\t')
# keep first 5 fetaureCounts columns: Chr, Start, End, Strand, Length
gene_length <- COUNT_MAT[,1:5]

# ---------------------
# Select list of DEGs |
# ---------------------
# Open the DESeq2 output for the selected contrast
comparison <- "mutd0_vs_WTd0" # example: mutd0 vs WTd0 contrast
df <- read.table(paste0("/DESeq2_path/",comparison,"/",comparison,"_DESq2.txt"), header = TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(df) <- df$Gene_name
df$DESeq <- comparison

# -----------------------------
# Assign gene lengths to DEGs |
# -----------------------------
df_up <- gene_length[which(rownames(gene_length) %in% df$Gene_name[df$label == "Up"]),]
df_down <- gene_length[which(rownames(gene_length) %in% df$Gene_name[df$label == "Down"]),]
df_ns <- gene_length[which(rownames(gene_length) %in% df$Gene_name[df$label == "NS"]),]

df_up$DEG <- paste0(comparison, "_up")
df_down$DEG <- paste0(comparison, "_down")
df_ns$DEG <- paste0(comparison, "_ns")

degs_gene_length <- rbind(df_up,df_down,df_ns)

# ------
# Plot |
# ------
# Violin + boxplot
ggplot(degs_gene_length, aes(x = DEG, y = log2(Length), fill = DEG)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA, na.rm = TRUE, adjust = 0.8) +
  geom_boxplot(width = 0.4, alpha = 0.4, position = position_dodge(0.5), outlier.shape = NA, na.rm = TRUE, color = "black") +
  labs(title = paste0("Gene sizes of DEGs"),
       x = "",
       y = expression(Log[2] ~ "(base pairs)")) +
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.background = element_blank(),  
    axis.text.y = element_text(size = 14, color = "black"),  
    axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1), 
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm")
  ) 