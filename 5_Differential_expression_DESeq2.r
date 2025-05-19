#!/usr/bin/Rscript

# This script includes codes for the following analyses:
# 1. DESeq2 analysis from featureCounts output
# 2. DESeq2 data visualization: PCA, volcano plot
# 3. GO-term analysis with clusterProfiler

# ----------------
# Load libraries |
# ----------------
library("dplyr")
library("tidyr")
library("tidyverse")
library("DESeq2")
library("vsn")
library("ggplot2")
library("clusterProfiler")
library("org.Mm.eg.db")

# ------------------------------------
# Open data from featureCouts output |
# ------------------------------------
cts <- read.csv("/featureCounts_path/featureCounts_file.csv", sep="\t",stringsAsFactors = F)

# --------------------
# Select column data |
# --------------------
cts2 <- cts[,2:ncol(cts)]

# --------------------------
# Create experiment labels |
# --------------------------
sample <- c("WTd0_rep1", "WTd0_rep2", "WTd1_rep1", "WTd1_rep2","WTd2_rep1", "WTd2_rep2", "WTd7_rep1", "WTd7_rep2",
            "mutd0_rep1", "mutd0_rep2", "mutd1_rep1", "mutd1_rep2","mutd2_rep1", "mutd2_rep2", "mutd7_rep1", "mutd7_rep2")
condition <- c("WTd0", "WTd0", "WTd1", "WTd1","WTd2", "WTd2", "WTd7", "WTd7",
               "mutd0", "mutd0", "mutd1", "mutd1","mutd2", "mutd2", "mutd7", "mutd7")
cell <- c("WT", "WT", "WT", "WT","WT", "WT", "WT", "WT",
          "mut", "mut", "mut", "mut","mut", "mut", "mut", "mut")
time <- c("d0", "d0", "d1", "d1", "d2", "d2","d7","d7",
          "d0", "d0", "d1", "d1", "d2", "d2","d7","d7")

# --------------------------------------
# Add rownames and convert into matrix |
# --------------------------------------
rownames(cts2) <- cts$Geneid
cts3 <-as.matrix(cts2)

# ------------------------------- 
# Construction of coldata table |
# -------------------------------
coldata <- data.frame(condition = condition, cell=cell, time=time)
rownames(coldata) <- sample
coldata$cell <- factor(coldata$cell)
all(rownames(coldata) %in% colnames(cts3))

# ---------------------------------------
# Create dds object for DESeq2 analysis |
# ---------------------------------------
dds <-DESeqDataSetFromMatrix(countData=cts3, 
                             colData=coldata, 
                             design = ~ condition)
dds <- DESeq(dds)

# -----------------
# Define contrast |
# -----------------
res2 <- results(dds, alpha=0.05, contrast=c("condition","mutd0","WTd0")) # example: mutd0 vs WTd0
res2 <- res2[order(res2$padj),]
summary(res2)

# ------------------------------
# Convert results to dataframe | 
# ------------------------------
res2_table <- as.data.frame(res2)
res2_table$Gene_name <- rownames(res2_table)

# Add significance label
for (i in 1:nrow(res2_table)) {
  if (!is.na(res2_table[i, "padj"]) & !is.na(res2_table[i, "log2FoldChange"])) {
    if (res2_table[i, "padj"] <= 0.05 & res2_table[i, "log2FoldChange"] > 0.5) {
      res2_table[i, "label"] <- "Up"} 
    else if (res2_table[i, "padj"] <= 0.05 & res2_table[i, "log2FoldChange"] < -0.5) {
      res2_table[i, "label"] <- "Down"} 
    else {
      res2_table[i, "label"] <- "NS"} 
    }
  if (is.na(res2_table[i, "padj"])) {
      res2_table[i, "label"] <- "remove"}
}

tosave <- res2_table[which(res2_table$label!="remove"),]

# Save DESeq info of all expressed genes 
comparison <- "mutd0_vs_WTd0" # example: mutd0 vs WTd0
setwd(paste0("/DESeq2_path/",comparison,"/"))
write.table(as.data.frame(tosave), file = paste0(comparison,"_DESq2.txt"),quote = F, sep="\t", col.names = T, row.names = F)

# ---------------------
# Data transformation |
# ---------------------
vsd <- vst(dds, blind=FALSE)

# --------
# PCA 2D |
# --------
vsd_matrix <- assay(vsd)
pca <- prcomp(t(vsd_matrix))
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
pcaData <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3])
pcaData$time <- rownames(pcaData) 

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=5) +
  geom_text(label=time,size=4, nudge_x =1.5, nudge_y =-.5, angle=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5, vjust = 0.5), 
    panel.background = element_rect(fill = "white", color = "black"))+
theme_classic()
dev.off()

# --------
# PCA 3D |
# --------
library("plotly")
library("reticulate")

system("pip install kaleido")
plot_ly(pcaData, x=~PC1, y=~PC2, z=~PC3, 
  type="scatter3d", mode="markers+text", marker=list(size = 8)
) %>%
  layout(
    scene = list(
      xaxis = list(title = paste0("PC1: ", percentVar[1], "%"),showgrid = FALSE),
      yaxis = list(title = paste0("PC2: ", percentVar[2], "%"),showgrid = FALSE),
      zaxis = list(title = paste0("PC3: ", percentVar[3], "%"),showgrid = FALSE)),
    title = list(text = "PCA plot", x = 0.5, y = 0.95, xanchor = 'center', yanchor = 'top', font = list(size = 15)),
    margin = list(l = 50, r = 50, b = 50, t = 50)
)

# --------------
# Volcano plot |
# --------------
# Gene name labels of top deregulated genes
topdegs <- head(tosave[order(tosave$padj), 'Gene_name'], 20)
tosave$Gene_labels <- ifelse(tosave$Gene_name %in% topdegs, tosave$Gene_name, NA)

ggplot(data=tosave, aes(x=log2FoldChange, y=-log10(padj), col=label, label=Gene_labels)) + 
  geom_point(size=3,alpha = 0.5) +
  labs(color = 'Condition',x= expression("log"[2]*" (fold change)"), y=expression("-log"[10]*" (p-adj)")) +
  geom_text(aes(label=Gene_labels), size=5, nudge_x =.9,nudge_y = 2.6, check_overlap = T) +  # Add gene labels
  theme_minimal() + 
  theme(
    axis.title.x = element_text(size = 14,colour = "black"),
    axis.title.y = element_text(size = 14,colour = "black"),
    axis.text.x = element_text(angle = 0, size = 14,colour = "black"),
    axis.text.y = element_text(hjust = 1, size = 14,colour = "black"),
    plot.title = element_text(size = 16,hjust = 0.5, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "right",
    plot.margin = margin(20, 20, 1, 10))

# ----------
# GO-terms |
# ----------
# Define list of Up or Down regulated genes
genesup_to_test <- res2_table[which(res2_table$label=="Up"),]$"Gene_name"
genesdown_to_test <- res2_table[which(res2_table$label=="Down"),]$"Gene_name"

# Select data to analyze
treatment <- "mutd0_vs_WTd0" # example: mutd0 vs WTd0

DEG <- "Up"
genes_to_test <- genesup_to_test

DEG <- "Down"
genes_to_test <- genesdown_to_test

# do GO analysis for Biological Processes
GO_results_BiologicalProcesses <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
as.data.frame(GO_results_BiologicalProcesses)
pdf(paste0("/path/GO_",DEG,"_",treatment,"_BP.pdf"), width = 10, height=10)
dotplot(GO_results_BiologicalProcesses, showCategory = 15, title=paste0("Biological Processes - ",DEG," in ",treatment), label_format = 30, font.size = 10)
dev.off()