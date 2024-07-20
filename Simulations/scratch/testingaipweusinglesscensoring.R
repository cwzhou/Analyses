df = pad_df
# Calculate the total count of people with status == 0
total_count <- df %>%
  filter(status == 0) %>%
  summarize(total_count = n())

# Extract the count value
total_count <- total_count$total_count

# Calculate the count of people to retain (i.e., 10% of the total count)
retain_count <- round(0.1 * total_count)

# Filter the data frame to retain only 10% of the people with status == 0
filtered_df <- df %>%
  filter(!(status == 0 & row_number() > retain_count))

# Now 'filtered_df' contains 90% of the people with status == 0 removed
table(filtered_df$status)
table(filtered_df$Trt)

filtered_df1 = filtered_df %>% 
  dplyr::select(-c("ischemia", "woundClass", "maxRutherfordClass", "Race")) %>%
  mutate(Trt = ifelse(Trt == "Open", 0, 1))

train_aipwe = aipwe_data_format(filtered_df1)
View(filtered_df1$Trt)
View(train_aipwe$data$A)

# getting rid of some 0s at end dataset
aif = aipwe.fit(data_list = train_aipwe,
                pp.v = (ncol(train_aipwe$Z1)-1)/2, #minus trt; divided by 2 b/c of interactions
                tau1 = as.numeric(t0_aipwe),
                tune = c(0.001,0.01,0.5,1,
                         seq(0.1,400,
                             length.out=16))) #eta0-eta_{ncov}

# package example dataset
aipwe.fit(data_list = tempD,
          pp.v = 30, #minus trt; divided by 2 b/c of interactions
          tau1 = 1,
          tune = c(0.001,0.01,0.5,1,
                   seq(0.1,400,
                       length.out=16)))

# Error in solve.default(h, z[[2]]) : 
#   system is computationally singular: reciprocal condition number = 2.02936e-16
# In addition: Warning message:
#   In coxph.fit(X, Y, istrat, offset, init, control, weights = weights,  :
#                  Ran out of iterations and did not converge