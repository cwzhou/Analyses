# Define your list of IDs
# id_list <- c(17, 4, 58, 73, 73, 24, 45, 36, 36, 36, 36, 12, 19, 52, 52)
# Define your new list of IDs
id_list <- c(25, 65, 41, 4, 58, 1, 73, 73, 24, 36, 36, 36, 36, 12, 19, 52, 52)

# Extract rows where ID matches any value in id_list
filtered_a0 <- a0[a0$ID %in% id_list, ]

# View the filtered dataset
# View(filtered_a0)

tp_surv = c(0, sort(timePointsSurvival), tau)

# Initialize results
at_risk <- numeric(length(tp_surv))
death_events <- numeric(length(tp_surv))

# Loop through each time point
for (i in seq_along(tp_surv)) {
  time <- tp_surv[i]

  # Find unique IDs of people at risk at this time point
  at_risk_ids <- unique(filtered_a0$ID[filtered_a0$L_open <= time & filtered_a0$R_closed >= time])
  at_risk[i] <- length(at_risk_ids)

  # Count the number of death events at this time point
  death_event_ids <- unique(filtered_a0$ID[filtered_a0$R_closed == time & filtered_a0$indD == 1])
  death_events[i] <- length(death_event_ids)
}

# Combine results into a data frame for clarity
risk_death_summary <- data.frame(
  TimePoint = tp_surv,
  AtRisk = at_risk,
  DeathEvents = death_events
)

# Print the results
print(risk_death_summary)
