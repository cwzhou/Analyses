home_dir = "~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper3_RE"
date_result_list <- "2025-03-05"
file_path <- file.path(home_dir, "output", date_result_list, "mff")

# Generate file indices, excluding 801-900
file_indices <- setdiff(1:1100, 801:900)

# Generate full file paths
file_names <- file.path(file_path, paste0("stacked.mff_sim", file_indices, ".rds"))

# Read all .rds files into a list
data_list <- lapply(file_names, readRDS)

# Optionally combine if they have the same structure
library(dplyr)
combined_data <- bind_rows(data_list)
View(combined_data)
