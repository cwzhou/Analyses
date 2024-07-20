survival_table = means_table$survival %>% as.data.frame()
endpoint_table = means_table$endpoint %>% as.data.frame()
# Combine tables side by side
combined_table <- dplyr::bind_cols(
  survival_table,
  endpoint_table
)

# Rename columns to have unique names
colnames(combined_table)[1:2] <- c("means_survival", "sds_survival")
colnames(combined_table)[3:4] <- c("means_endpoint", "sds_endpoint")
combined_table = combined_table %>% mutate(diff_means = means_survival-means_endpoint)
colnames(combined_table) <- NULL
# Convert to xtable format
xtable_format <- xtable::xtable(combined_table)

# Print the LaTeX code
print(xtable_format, include.rownames = FALSE)

print(xtable_format,
      only.contents=TRUE,
      include.rownames=TRUE,
      type="latex",
      # digits(tbl) <- c(4,4,4),
      file="Scratch/tblout.tex")


tbl <- xtable(means_table$survival,
              caption = sprintf("Mean and Sd of Survival Endpoint for all Simulations (n.sim = %s)", full$settings$n.sim))
