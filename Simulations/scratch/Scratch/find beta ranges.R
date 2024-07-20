library(beepr)
# simple_exp setting with 2 covariates
# lambda11 = exp(z1*beta111+z2*beta121);
# lambda10 = exp(z1*beta110+z2*beta120). 
# lambda21 = exp(z1*beta211+z2*beta221);
# lambda20 = exp(z1*beta210+z2*beta220). 
# Lam1 = lambda11+lambda21. 
# Lam0 = lambda10+lambda11. 
# beta110 = -1, beta120 = -1, beta111 = -1, beta121 = 0.2. 
# t0 = 0.6. 
# z1 and z2 are generated from random standard uniform distribution. 
# WANT TO FIND ranges of beta211,221,210,220 such that the following equation is always true: 
# (lambda10/Lam0) - exp(-t0*lambda20)+(lambda20/Lam0)*exp(-t0*Lam0) < 
#   (lambda11/Lam1) - exp(-t0*lambda21)+(lambda21/Lam1)*exp(-t0*Lam1). 

# Given values
beta110 <- -0.3 # treatment 0
beta120 <- -0.5 # treatment 0
beta111 <- 0.8 # treatment 1
beta121 <- -0.2 # treatment 1
t0 <- 0.6

# Generate random values for z1 and z2
set.seed(123)  # Set seed for reproducibility
og_z1 <- runif(1)
og_z2 <- runif(1)

# Function to check the inequality
check_cif_inequality <- function(z1, z2, beta210, beta220, beta211, beta221) {
  lambda10 <- exp(z1 * beta110 + z2 * beta120)
  lambda20 <- exp(z1 * beta210 + z2 * beta220)
  lambda11 <- exp(z1 * beta111 + z2 * beta121)
  lambda21 <- exp(z1 * beta211 + z2 * beta221)
  
  Lam0 <- lambda10 + lambda11
  Lam1 <- lambda11 + lambda21
  
  left_side <- (lambda10 / Lam0) - exp(-t0 * lambda20) + (lambda20 / Lam0) * exp(-t0 * Lam0)
  right_side <- (lambda11 / Lam1) - exp(-t0 * lambda21) + (lambda21 / Lam1) * exp(-t0 * Lam1)
  
  return(left_side < right_side)
}

# Grid search for beta211, beta221, beta210, beta220
possible_values_210 <- seq(0.1, 1.5, by = 0.1)
possible_values_220 <- seq(-1, 1, by = 0.1)
possible_values_211 <- seq(-1, 1, by = 0.1)
possible_values_221 <- seq(-1.5, -0.1, by = 0.1)
valid_combinations0 <- expand.grid(beta210 = possible_values_210, 
                                  beta220 = possible_values_220,
                                  beta211 = possible_values_211, 
                                  beta221 = possible_values_221,
                                  beta1110 = beta110, 
                                  beta120 = beta120, 
                                  beta111 = beta111, 
                                  beta121 = beta121, 
                                  z1 = og_z1, 
                                  z2 = og_z2, 
                                  t0 = 0.6)
possible_values1 <- seq(0, 1, by = 0.1)[-1]
z_combinations <- expand.grid(z1 = possible_values1, 
                              z2 = possible_values1)

# Check the inequality for each combination
valid_combinations1 <- valid_combinations0[apply(valid_combinations0, 1, function(row) {
  check_cif_inequality(row[9], row[10], row[1], row[2], row[3], row[4])
}), ]
# Print the valid combinations
print(head(valid_combinations1)); View(valid_combinations1)


library(dplyr)
# Create a new variable 'abs_diff' representing the absolute difference
your_dataset <- valid_combinations1 %>%
  mutate(abs_diff_b2a0 = abs(beta211 - beta221),
         abs_diff_b2a1 = abs(beta210 - beta220))
max20 = max(your_dataset$abs_diff_b2a0) - 0.5
max21 = max(your_dataset$abs_diff_b2a1) - 0.5
# Find the row(s) where the absolute difference is the largest
max_diff_rows <- your_dataset %>% 
  filter(abs_diff_b2a0 >= max20 | abs_diff_b2a1 >= max21)

# Print or use the resulting dataset
print(max_diff_rows)

sample_num = 10
valid_combinations = max_diff_rows[sample(nrow(max_diff_rows), sample_num), ]
valid_combinations

new_beta2_poss = data.frame()
count = 0
for (ix in 1:nrow(valid_combinations)){
  tmp_beta2 = valid_combinations[ix,1:4]
  logic <- do.call(check_cif_inequality, c(list(z1 = og_z1, z2 = og_z2), as.list(tmp_beta2)))
  if (!logic){
    stop("error")
  } else{
    test = data.frame(i = 1:nrow(z_combinations), new = NA)
    for (i in 1:nrow(z_combinations)){
      z1 <- z_combinations[i,1]
      z2 <- z_combinations[i,2]
      logic_newz = do.call(check_cif_inequality, 
                           c(list(z1 = z1, z2 = z2), 
                             as.list(tmp_beta2)))
      test[i, "new"] = logic_newz
      # if (!logic_newz){
      #   message(i,": this beta2 setting doesn't work for iteration ", i)
      #   cat("z1:",z1)
      #   cat("\nz2:",z2,"\n")
      # } # end of !logic_newz
  
    } # end of for 1:nrow(z_combinations)

    if (all(test$new)) {
      count = count+1
      message("Iteration ", ix)
      print(paste("All Z combinations are valid: beta2:", toString(tmp_beta2)))
      new_beta2_poss[count,1] = ix
      new_beta2_poss[count,2:5] = tmp_beta2
    } #else {
      # print(paste("Some Z combinations are not valid for iteration", ix))
    # }
  } # end of else
} # end of for 1:nrow(valid combination)
new_beta2_poss$beta110 = beta110
new_beta2_poss$beta120 = beta120
new_beta2_poss$beta111 = beta111
new_beta2_poss$beta121 = beta121
new_beta2_poss$t0 = t0
colnames(new_beta2_poss)[1] = "Iteration"
View(new_beta2_poss)
beep(2)
