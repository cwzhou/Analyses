fold_trt_df <- data.frame(
  obs = testing_dataset_surv$A,
  czmk = opt.rule.CZMK[[1]],  # Extract the first (or correct) column
  zom = opt.rule.ZOM[[1]]     # Extract the first (or correct) column
)

library(dplyr)
summary_table_fold <- data.frame(
  Trt = c("Trt 0", "Trt 1", "Trt 2"),
  CZMK = c(round(mean(fold_trt_df$czmk == 0) * 100, 2), round(mean(fold_trt_df$czmk == 1) * 100, 2), round(mean(fold_trt_df$czmk == 2) * 100, 2)),
  ZOM = c(round(mean(fold_trt_df$zom == 0) * 100, 2), round(mean(fold_trt_df$zom == 1) * 100, 2), round(mean(fold_trt_df$zom == 2) * 100, 2)),
  Observed = c(round(mean(fold_trt_df$obs == 0) * 100, 2), round(mean(fold_trt_df$obs == 1) * 100, 2), round(mean(fold_trt_df$obs == 2) * 100, 2))
)

# Print the summary table to check
print(summary_table_fold)

