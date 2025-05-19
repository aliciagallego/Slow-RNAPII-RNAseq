#!/usr/bin/Rscript

# This script includes codes for the following functions:
# 1. Define a function to remove Ns values from vast-tools INCLUSION_LEVELS.tab output table (quality filter)
# 2. Retrieve the list of PSI values of each event in all samples computed by vast-tools combine (own and other datasets)
# 3. Unify different types of exon events in one single category
# 4. Filter to keep only AS events (PSI values between 10%-90%)
# 5. Plot total numbers of AS events per sample

# -----------
# Libraries |
# -----------
library("tidyr")
library("dplyr")
library("ggplot2")
library("ggbreak")

# -----------------
# Filter function |
# -----------------
# Remove Ns from vast-tools output
removeNs = function(x){
  scores <- grep(".Q",names(x),fixed = T)
  score2 <- scores - 1
  for(i in score2){x[,i] <- ifelse(grepl("^N",x[,i+1]),NA,x[,i])}
  return(x)
}

# --------------------------------
# Open vast-tools combine output |
# --------------------------------
# Own dataset
setwd("/path/Splicing/vasttools/vast_out/")
df = read.delim("INCLUSION_LEVELS_FULL-mm10-8.tab",stringsAsFactors = F)
df2 = removeNs(df)
df2 = df2[,c(1:6,seq(7,ncol(df2),2))]
df3 = df2[complete.cases(df2),]

# Maslon et al. dataset
setwd("/path/Splicing/vasttools_Maslon/vast_out/")
df_maslon = read.delim("INCLUSION_LEVELS_FULL-mm10-2.tab",stringsAsFactors = F)
df_maslon2 = removeNs(df_maslon)
df_maslon2 = df_maslon2[,c(1:6,seq(7,ncol(df_maslon2),2))]
df_maslon3 = df_maslon2[complete.cases(df_maslon2),]

# --------------------------------------
# Total number of events - Own dataset |
# --------------------------------------
# Reunify all exon events ("ANN", "S", "C3", "C2", "C1", "MIC") 
df4 <- df3
df4 <- df4 %>%
  mutate(COMPLEX = ifelse(COMPLEX %in% c("ANN", "S", "C3", "C2", "C1", "MIC"), "Exon", 
                          ifelse(COMPLEX == "IR", "Intron", COMPLEX)))

# Quantify AS events (PSI 10-90%) in one table 
cols_to_check <- c("mutd1", "mutd2", "mutd7", "mutd0", "WTd1", "WTd2", "WTd7", "WTd0")

results_df <- do.call(rbind, lapply(cols_to_check, function(col) {
  df4 %>%
    filter(.data[[col]] > 10 & .data[[col]] < 90) %>%
    group_by(COMPLEX) %>%
    summarize(count = n(), .groups = 'drop') %>%
    mutate(Column = gsub("_", " ", col))
}))

# Filter to keep samples of interest: example for wtd0 and mutd0
filtered_results_df <- results_df %>%
  filter(Column %in% c("mutd0",  "WTd0"))

# Select data for plotting
data <- filtered_results_df

# Bar Plot 
ggplot(data, aes(x = COMPLEX, y = count, fill = Column)) +
  geom_bar(stat = "identity", position = position_dodge(width = .8), width = 0.7) + 
  scale_x_discrete(expand = expansion(mult = 0.3)) + 
  labs(title = "Total AS events (PSI 10-90%) - own data d0",
       x = "",
       y = "Number of AS events",
       fill = "Column") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.text.y = element_text(size = 14, color = "black"),  
    axis.text.x = element_text(size = 14, color = "black", angle = 0, hjust = .5), 
    axis.title.x = element_text(size = 14, color = "black"),  
    axis.title.y = element_text(size = 14, color = "black"),  
    axis.line = element_line(color = "black"),  
    axis.line.x = element_line(colour = "black", linewidth = .5),
    axis.line.y = element_line(colour = "black", linewidth = .5),
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.15, "cm")  
  ) 

# ------------------------------------------------
# Total number of events - Maslon et al. dataset |
# ------------------------------------------------
# Reunify all exon events ("ANN", "S", "C3", "C2", "C1", "MIC") 
df_maslon4 <- df_maslon3
df_maslon4 <- df_maslon4 %>%
  mutate(COMPLEX = ifelse(COMPLEX %in% c("ANN", "S", "C3", "C2", "C1", "MIC"), "Exon", 
                          ifelse(COMPLEX == "IR", "Intron", COMPLEX)))

# Quantify AS events (PSI 10-90%) in one table 
cols_to_check <- c("mut", "WT")

results_df <- do.call(rbind, lapply(cols_to_check, function(col) {
  df_maslon4 %>%
    filter(.data[[col]] > 10 & .data[[col]] < 90) %>%
    group_by(COMPLEX) %>%
    summarize(count = n(), .groups = 'drop') %>%
    mutate(Column = gsub("_", " ", col))
}))

# Select data for plotting
data <- results_df

# Bar Plot (ampty color bars)
ggplot(data, aes(x = COMPLEX, y = count, fill = "white", color = Column)) +
  geom_bar(stat = "identity", position = position_dodge(width = .8), width = 0.7, size = 1) + 
  scale_fill_manual(values = c("white")) + 
  scale_x_discrete(expand = expansion(mult = 0.3)) +
  labs(title = "Alternative spliced events (PSI 10-90%) - Maslon data",
       x = "",
       y = "Number of events",
       color = "Column") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.background = element_blank(),  
    axis.text.y = element_text(size = 14, color = "black"),  
    axis.text.x = element_text(size = 14, color = "black", angle = 0, hjust = .5),  
    axis.title.x = element_text(size = 14, color = "black"),  
    axis.title.y = element_text(size = 14, color = "black"),  
    axis.line = element_line(color = "black"),  
    axis.line.x = element_line(colour = "black", linewidth = .5),
    axis.line.y = element_line(colour = "black", linewidth = .5),
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.15, "cm")  
  ) +
  scale_y_break(c(100, 3500))  # Break between 3000 and 3100