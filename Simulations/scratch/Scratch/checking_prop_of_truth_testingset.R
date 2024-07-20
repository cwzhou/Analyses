# check how many of the dataset for cause 1 is equal to true best treatmnet
#cause1
check_cause1_prop = function(dataset, cause, comp_dataset){# Find row indices in dataset where var1 == 1
row_indices <- which(dataset$status == cause)
# Index dataset based on the row indices
filtered_dataset <- dataset[row_indices, ]
# Index comp_dataset based on the row indices
filtered_comp_dataset <- comp_dataset[row_indices, ]

# Find the mean
result_mean <- mean(filtered_dataset$action == filtered_comp_dataset$best.action)
return(result_mean)
}

for (cause in 1:2){
  message("\n\nCause",cause,": OBS vs TRUTH for TESTING: ", round(check_cause1_prop(rep_obs,cause,test_truth),2))
  message("Cause",cause,": CZMK vs TRUTH for TESTING: ", round(check_cause1_prop(rep_czmk,cause,test_truth),2))
  message("Cause",cause,": CSK vs TRUTH for TESTING: ", round(check_cause1_prop(rep_csk,cause,test_truth),2))
  message("Cause",cause,": ZOM vs TRUTH for TESTING: ", round(check_cause1_prop(rep_zom,cause,test_truth),2))
}

