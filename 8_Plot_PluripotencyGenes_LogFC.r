#!/usr/bin/Rscript

# This script includes codes for the following functions:
# 1. Define a list of naive pluripotency genes described in:
# Lynch CJ, et al. Global hyperactivation of enhancers stabilizes human and mouse naive pluripotency through inhibition of CDK8/19 Mediator kinases. Nat Cell Biol. 2020 Oct;22(10):1223-1238. 
# doi: 10.1038/s41556-020-0573-1
# 2. Retrieve the list of DEGs computed by DESeq2
# 3. Assign Log2(FC) values computed by DESeq2 to the pluripotency gene list
# 4. Plot heatmap of Log2(FC) values in pluripotency genes


# ----------------
# Load libraries |
# ----------------
library("pheatmap")

# --------------------
# Pluripotency genes |
# --------------------
pluripotent_genes <- c("Lin28a","Nr0b1","Myc","L1td1","Eras","Zfp42","Podxl","Pou5f1","Sox2","Tcl1",
                       "Pecam1","Dppa5a","Alpl","Tbx3","Klf4","Esrrb","Tdh","Nanog","Tfcp2l1")

# ---------------------
# Select list of DEGs |
# ---------------------
# Open the DESeq2 output for the selected contrasts
# examples for: 'd1 vs d0', 'd2 vs d0', and 'd7 vs d0' contrasts
list <- c("wtd1_wtd0", "wtd2_wtd0", "wtd7_wtd0","mutd1_mutd0", "mutd2_mutd0", "mutd7_mutd0")

for (i in list){
  df <- read.table(paste0("/DESeq2_path/",comparison,"/",comparison,"_DESq2.txt"), header = TRUE, sep="\t", stringsAsFactors=FALSE)
  rownames(df) <- df$Gene_name
  df_pluri <- df[pluripotent_genes,]
  assign(i,df_pluri)
}

# ------------------------------------------------------
# Assign Log2(FC) values to the pluripotency gene list |
# ------------------------------------------------------
# Combine all contrast of interest and get log2FoldChange values
all_samples_pluri <- cbind(wtd1_wtd0[,c(7,2)],wtd2_wtd0[,2],wtd7_wtd0[,2],
                           mutd1_mutd0[,2],mutd2_mutd0[,2],mutd7_mutd0[,2])
all_samples_pluri <- all_samples_pluri[,c(2:7)]

# Rename samples
colnames(all_samples_pluri) <- c("WT d1","WT d2","WT d7","mut d1","mut d2","mut d7")

# ------
# Plot |
# ------
# heatmap
pheatmap(t(all_samples_pluri),
         scale="none",
         show_rownames=T,
         show_colnames=T,
         cluster_rows=F,
         cluster_cols=T,
         gaps_row = c(3),
         cellwidth=20,
         cellheight=15,
         annotation_names_row=TRUE,
         legend=T,
         breaks = seq(-2.5, 2.5, length.out = 100),
         legend_labels = "Log2[FoldChange]",
         display_numbers=F,
         fontsize_row=12,
         fontsize_col=12,
         angle_col = "45",
         main = "Pluripotency genes")