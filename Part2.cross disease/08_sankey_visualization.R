library(ggplot2)
library(ggalluvial)
library(readr)
library(dplyr)

# ==============================================================================
# Configuration
# ==============================================================================
# Resolve paths relative to script location
# Assumes script is in: [Project_Root]/scripts/Part2_CrossDisease/
args <- commandArgs(trailingOnly=FALSE)
file_arg <- grep("--file=", args, value=TRUE)

if (length(file_arg) > 0) {
  script_dir <- dirname(sub("--file=", "", file_arg))
} else {
  script_dir <- getwd()
}

base_dir <- normalizePath(file.path(script_dir, "..", ".."))
data_dir <- file.path(base_dir, "results", "part2_cross_disease", "MEBOCOST")
plot_dir <- file.path(data_dir, "SankeyPlots")

if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ==============================================================================
# Data Loading
# ==============================================================================
f_sender   <- file.path(data_dir, "SourceData_Fig3H_Sender.csv")
f_receiver <- file.path(data_dir, "SourceData_Fig3G_Receiver.csv")

if(!file.exists(f_sender)) stop("Sender data file missing (SourceData_Fig3H_Sender.csv).")
if(!file.exists(f_receiver)) stop("Receiver data file missing (SourceData_Fig3G_Receiver.csv).")

df_sender   <- read_csv(f_sender, show_col_types = FALSE)
df_receiver <- read_csv(f_receiver, show_col_types = FALSE)

# ==============================================================================
# Visualization
# ==============================================================================
generate_sankey <- function(data, plot_title) {
  
  if(!"Count" %in% colnames(data)) data$Count <- 1
  
  cols <- c("Healthy"="#4DAF4A", "RA"="#E41A1C", "PSA"="#984EA3", "SPA"="#377EB8")
  
  p <- ggplot(data,
              aes(axis1 = Status, axis2 = Sender, axis3 = Sensor, 
                  axis4 = Metabolite, axis5 = Receiver, y = Count)) +
    geom_alluvium(aes(fill = Status), width = 1/12, alpha = 0.7, knot.pos = 0.3) +
    geom_stratum(width = 1/6, fill = "grey95", color = "black", alpha = 0.8) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
              size = 3, fontface = "bold", min.y = 3) +
    scale_x_discrete(limits = c("Status", "Sender", "Sensor", "Metabolite", "Receiver"), 
                     expand = c(.05, .05)) +
    scale_fill_manual(values = cols) +
    labs(title = plot_title, y = "Interaction Frequency") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          legend.position = "none")
  
  return(p)
}

# Generate and Save
p1 <- generate_sankey(df_sender, "RTE Metabolic Secretion (RTE as Sender)")
ggsave(file.path(plot_dir, "Fig3H_Sankey_Sender.pdf"), p1, width = 14, height = 8)

p2 <- generate_sankey(df_receiver, "RTE Metabolic Uptake (RTE as Receiver)")
ggsave(file.path(plot_dir, "Fig3G_Sankey_Receiver.pdf"), p2, width = 14, height = 8)

message("Plots generated in: ", plot_dir)
