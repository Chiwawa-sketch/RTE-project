# scripts/04_correlation.R
# Load required libraries
library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr) 
# ==============================================================================
# 1. Load Data
# ==============================================================================
if (!file.exists("results/tcell_with_rte.rds")) {
  stop("Input file 'results/tcell_with_rte.rds' not found. Please run script 03 first.")
}
tcell <- readRDS("results/tcell_with_rte.rds")
# ==============================================================================
# 2. Define Age Metadata (Young & Old ONLY)
# ==============================================================================
sample_metadata <- data.frame(
  sample_id = c(
    # Young (Ages: 27, 42, 23)
    "GSM4750303", "GSM4750304", "GSM4750305",
    # Old (Ages: 78, 72, 80, 88, 97, 100)
    "GSM4750309", "GSM4750310", "GSM4750311", 
    "GSM5684306", "GSM5684307", "GSM5684308"
  ),
  True_Age = c(
    27, 42, 23,           
    78, 72, 80, 88, 97, 100 
  )
)

# Merge metadata
meta_data <- tcell@meta.data
meta_data$barcode <- rownames(meta_data)
meta_data <- meta_data %>% left_join(sample_metadata, by = "sample_id")
rownames(meta_data) <- meta_data$barcode
tcell@meta.data <- meta_data
# ==============================================================================
# 3. Calculate Proportions (All Samples)
# ==============================================================================
naive_meta <- tcell@meta.data %>% filter(celltype == "Naive T")
target_genes <- c("SOX4", "IKZF2", "TOX", "TOX2", "NREP", 
                  "AUTS2", "PECAM1", "CD38", "ZMIZ1", "PTK7")
prop_data_list <- list()
for (g in target_genes) {
  pos_string <- paste0(g, "+ Naive T")
  df <- naive_meta %>%
    # Calculate for ALL groups (including Cord/Frail which have NA Age)
    group_by(sample_id, True_Age, Group) %>%
    summarise(
      Total_Naive = n(),
      Count_Pos = sum(.data[[g]] == pos_string, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Percent_Positive = (Count_Pos / Total_Naive) * 100,
      Gene = g
    )
  prop_data_list[[g]] <- df
}
all_samples_df <- bind_rows(prop_data_list)
# ==============================================================================
# 4. Export Data for Figure 1G
# ==============================================================================
donut_data <- all_samples_df %>%
  group_by(Group, Gene) %>%
  summarise(
    Mean_Percent = mean(Percent_Positive),
    SD_Percent = sd(Percent_Positive),
    .groups = "drop"
  )
write_csv(donut_data, "results/SourceData_FigureG_Donuts.csv")
# ==============================================================================
# 5. Plot Correlations (supFigure 1f) - YOUNG AND OLD ONLY
# ==============================================================================
message("Generating supFigure 1f (Correlation)...")
cor_df <- all_samples_df %>% filter(!is.na(True_Age))
write_csv(cor_df, "results/SourceData_supFig1f_Correlation.csv")
pdf("results/supFig1f_Correlations.pdf", width = 15, height = 8)
p <- ggplot(cor_df, aes(x = Percent_Positive, y = True_Age)) +
  geom_point(aes(color = Group), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", fill = "grey", alpha = 0.2, linewidth = 0.8) +
  stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "top", size = 3.5) + 
  facet_wrap(~Gene, scales = "free_x", ncol = 5) +
  theme_classic(base_size = 14) +
  labs(title = "Correlation: Age vs RTE+ Frequency (Young & Old Only)",
       x = "RTE+ Cell Percentage (%)", y = "Donor Age (Years)") +
  scale_color_manual(values = c("Young"="#4DAF4A", "Old"="#E41A1C"))
print(p)
dev.off()
# ==============================================================================
# 6. Plot Figure 1G (Donuts/Radial Bars)
# ==============================================================================
message("Generating Figure 1G (Donut Charts)...")
# Create a sub-directory for these plots
if(!dir.exists("results/Figure1G_Donuts")) dir.create("results/Figure1G_Donuts")
# Set colors matching the original figure
# Cord=Blue, Young=Green, Old=Cyan/Blue-ish, Frail=Red
group_colors <- c("Cord_blood" = "#377EB8", 
                  "Young"      = "#00BA38", 
                  "Old"        = "#619CFF", 
                  "Frail"      = "#F8766D")

# Define plot order
donut_data$Group <- factor(donut_data$Group, levels = c("Cord_blood", "Young", "Old", "Frail"))
# Loop through each gene to create a separate Donut plot
for (g in target_genes) {
  # Filter data for current gene
  plot_data <- donut_data %>% filter(Gene == g)
  # We construct a 'rect' for the positive part and a background track if needed.
  # Logic: ymax = Percentage, ymin = 0.
  p_donut <- ggplot(plot_data, aes(x = 2, y = Mean_Percent, fill = Group)) +
    # The colored bar representing the percentage
    geom_bar(stat = "identity", color = "white") +
    # Limits for the y-axis (0 to 100%)
    ylim(0, 100) +
    # Transform to polar coordinates to make it a donut
    coord_polar(theta = "y", start = 0) +
    # Facet to separate the 4 groups nicely
    facet_wrap(~Group, ncol = 4) +
    # Manual Colors
    scale_fill_manual(values = group_colors) +
    # Limits for x-axis creates the "Hole" in the middle
    xlim(0.5, 2.5) +
    theme_void() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    labs(title = g) +
    # Add text label in the middle
    geom_text(aes(label = paste0(round(Mean_Percent, 2), "%"), y = 0), 
              x = 0.5, size = 4, fontface = "bold", color = "black")
  # Save Plot
  ggsave(filename = file.path("results/Figure1G_Donuts", paste0("Donut_", g, ".pdf")), 
         plot = p_donut, width = 8, height = 3)
}

message("Analysis Complete.")
message("1. supFigure 1f (Correlation) saved to: results/supFig1f_Correlations.pdf")
message("2. Figure 1G (Donuts) saved to: results/Figure1G_Donuts/")