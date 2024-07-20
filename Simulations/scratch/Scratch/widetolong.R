library(reshape2)
# Reshape the data into long format
df_long <- melt(r00, id = "sim")

# Extract method, A, and B from variable column
df_long2 <- transform(df_long,
                      method = sub("^(.*?)_.*$", "\\1", variable),
                      data_type = sub("^.*?_(.*?)_.*$", "\\1", variable),
                      endpoint_type = sub("^.*_([^_]+)$", "\\1", variable),
                      proptrt = value)

# Create the final data frame
df_long3 <- df_long2[, c("method", "data_type", "endpoint_type", "proptrt")]
View(df_long3)
