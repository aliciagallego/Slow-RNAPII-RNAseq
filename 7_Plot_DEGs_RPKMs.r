#!/usr/bin/Rscript

# This script includes codes for the following functions:
# 1. Compute normalized gene expression levels (RPKMs) in the featureCounts output
# 2. Retrieve the list of DEGs computed by DESeq2
# 3. Assign normalized expression levels (RPKMs) to DEGs
# 4. Plot RPKM distributions in Up, Down and n.s. DEGs

# ----------------
# Load libraries |
# ----------------
library("ggplot2")

# -------------------------------
# featureCounts output to RPKMs |
# -------------------------------

# Open featureCounts output
# (first 5 columns must include: Chr, Start, End, Strand, Length)
COUNT_MAT=read.table('/featureCounts_path/featureCounts_file.txt',row.names=1,header=TRUE,sep='\t')

# Keep gene length and gene counts for all replicates
COUNT_MAT=COUNT_MAT[,(5:ncol(COUNT_MAT))]

# Take gene distance
L=COUNT_MAT[,1]

# Take gene counts
COUNT_MAT=COUNT_MAT[,2:ncol(COUNT_MAT)]

# Define RPKM function:
.rpkm=function(x, L){
  y= x *1000 * 1000000 / L / sum(x); 
  return(y)
}

# Apply RPK function to featureCounts output
NORM_COUNT_MAT=apply(COUNT_MAT, 2,.rpkm, L)
head(NORM_COUNT_MAT)

# ---------------------
# Select list of DEGs |
# ---------------------
# Open the DESeq2 output for the selected contrast
comparison <- "mutd0_vs_WTd0" # example: mutd0 vs WTd0 contrast
df <- read.table(paste0("/DESeq2_path/",comparison,"/",comparison,"_DESq2.txt"), header = TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(df) <- df$Gene_name
df$DESeq <- comparison

# ----------------------
# Assign RPKMs to DEGs |
# ----------------------
# Get RPKMs for Up, Down, and n.s. DEGs
df_up <- as.data.frame(NORM_COUNT_MAT[which(rownames(NORM_COUNT_MAT) %in% df$Gene_name[df$label == "Up"]),])
df_down <- as.data.frame(NORM_COUNT_MAT[which(rownames(NORM_COUNT_MAT) %in% df$Gene_name[df$label == "Down"]),])
df_ns <- as.data.frame(NORM_COUNT_MAT[which(rownames(NORM_COUNT_MAT) %in% df$Gene_name[df$label == "NS"]),])
df_up$DEG <- paste0(comparison, "_up")
df_down$DEG <- paste0(comparison, "_down")
df_ns$DEG <- paste0(comparison, "_ns")
  
# Generate one data frame per sample
for (k in 1:ncol(NORM_COUNT_MAT)) {
  table_up <- df_up[, c(k, ncol(df_up))]
  table_down <- df_down[, c(k, ncol(df_down))]
  table_ns <- df_ns[, c(k, ncol(df_ns))]
  
  new_label <- colnames(df_up[k])
  table_up$sample <- new_label
  table_down$sample <- new_label
  table_ns$sample <- new_label
    
  colnames(table_up)[1] <- "RPKM"
  colnames(table_down)[1] <- "RPKM"
  colnames(table_ns)[1] <- "RPKM"
  
  table_updownns <- rbind(table_up,table_down,table_ns)
  assign(paste0(comparison,"_",new_label), table_updownns)
}
# Expected data frame names for this example: wtd0, mutd0

# ------
# Plot |
# ------
# Combine samples of interest to plot
data <- rbind(wtd0,mutd0) # example: mutd0 vs WTd0 contrast
data$DEG_sample <- interaction(data$DEG, data$sample)

# Violin + boxplot
ggplot(data, aes(x = DEG, y = log2(RPKM), fill = sample)) +
  geom_violin(trim = FALSE, alpha = 0.4, colour = NA, na.rm = TRUE, adjust = 0.7) +  
  geom_boxplot(width = 0.4, alpha = 0.5, position = position_dodge(0.9), outlier.shape = NA, na.rm = TRUE, color = "black") + 
  labs(title = paste0("Normalized expression of DEGs"),
       x = "", y = "RPKM") +  
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