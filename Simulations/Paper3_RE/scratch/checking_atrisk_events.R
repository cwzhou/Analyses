id_list <- c(
  22, 22, 22, 22, 35, 258, 57, 57, 96, 23, 33, 33, 260, 260, 97, 97, 97, 97,
  100, 242, 210, 171, 41, 29, 29, 81, 81, 81, 233, 16, 2, 2, 2, 2, 47, 90, 71,
  71, 294, 271, 145, 280, 148, 159, 228, 235, 24, 300, 298, 26, 26, 30, 180,
  116, 230, 17, 8, 217, 217, 72, 72, 72, 186, 133, 192, 192, 192, 291, 291, 45,
  94, 94, 94, 94, 252, 252, 252, 252, 182, 182, 182, 281, 161, 161, 46, 46, 259,
  68, 143, 53, 53, 53, 64, 64, 64, 64, 102, 102, 102, 248, 248, 248, 263, 264,
  264, 264, 264, 99, 99, 99, 99, 282, 282, 104, 136, 136, 136, 136, 15, 7, 7,
  285, 220, 183, 183, 183, 37, 37, 37, 37, 83, 137, 137, 137, 266, 266, 201, 201,
  208, 61, 61, 61, 202, 288, 288, 105, 221, 251, 85, 146, 240
)



# id_list_full <- c(
#   5, 13, 13, 15, 15, 15, 15, 18, 24, 26, 27, 27, 27, 27, 29, 29, 29, 30, 34, 34, 34, 34, 35, 35, 36, 36,
#   39, 40, 41, 41, 44, 44, 44, 44, 47, 47, 49, 51, 56, 56, 56, 56, 58, 59, 60, 61, 61, 63, 63, 63, 63, 64,
#   65, 68, 68, 68, 68, 69, 69, 69, 69, 71, 71, 71, 74, 75, 76, 76, 76, 77, 78, 80, 84, 91, 96, 96, 100, 100,
#   101, 101, 101, 101, 104, 104, 106, 109, 111, 113, 113, 113, 113, 120, 120, 122, 124, 128, 129, 131, 131,
#   131, 133, 135, 135, 135, 135, 140, 140, 142, 142, 146, 154, 156, 160, 161, 162, 164, 165, 166, 166, 166,
#   172, 178, 178, 180, 180, 180, 181, 181, 182, 184, 184, 186, 187, 188, 188, 188, 189, 190, 190, 190, 192,
#   193, 193, 199, 199, 199, 199, 200, 200, 200, 200, 202, 202, 202, 203, 205, 206, 208, 209, 209, 216, 220,
#   220, 220, 223, 226, 230, 230, 230, 230, 234, 234, 234, 235, 237, 237, 237, 237, 238, 239, 240, 241, 241,
#   241, 241, 243, 245, 245, 245, 248, 249, 252, 253, 253, 254, 254, 254, 254, 255, 259, 261, 262, 268, 268,
#   269, 269, 269, 269, 271, 273, 274, 275, 279, 281, 282, 282, 282, 282, 283, 288, 290, 290, 290, 290, 291,
#   291, 291, 292, 292, 292, 292, 294, 295, 295, 297, 298, 298, 299, 300, 300, 300, 300
# )

# View(data_to_use %>%
#        filter(ID %in% id_list_full) %>%
#        dplyr::select(ID, Z1, IndD, IndR))

# data_to_use %>%
#   filter(ID %in% id_list_full) %>%
#   # group_by(ID) %>%
#   mutate(ID_new = ifelse(row_number() == 1, row_number(), 0)) %>%
#   mutate(ID_new = ifelse(lag(ID) == ID, lag(ID_new),
#                          ID_new+1)) %>%  # Ensure ID_new is sequential for each ID group
#   dplyr::select(ID_new, ID, Z1, IndD, IndR) %>%
#   View()

# data_to_use %>%
#   filter(ID %in% id_list_full) %>%  # Filter based on ID list
#   mutate(ID_new = 1 + cumsum(lag(ID, default = first(ID)) != ID)) %>%  # Start ID_new from 1
#   dplyr::select(ID_new, ID, Z1, IndD, IndR) %>%  # Select the desired columns
#   View()




# Extract rows where ID matches any value in id_list
filtered_a0 <- a0[a0$ID %in% id_list, ]
# View(filtered_a0)

tp_surv = c(0, sort(timePointsSurvival), tau0)
tp_end = c(0, sort(timePointsEndpoint), tau0)

# Define the function to summarize risk events
summarize_risk <- function(tp_surv, indicator_column, time_column, event_column, filtered_data) {
  # Initialize results
  at_risk <- numeric(length(tp_surv))
  event_count <- numeric(length(tp_surv))

  # Loop through each time point
  for (i in seq_along(tp_surv)) {
    time <- tp_surv[i]

    # Find unique IDs of people at risk at this time point
    at_risk_ids <- unique(filtered_data$ID[filtered_data[[time_column]] <= time & filtered_data[[event_column]] >= time])
    at_risk[i] <- length(at_risk_ids)

    # Count the number of events at this time point based on the indicator
    event_ids <- unique(filtered_data$ID[filtered_data[[event_column]] == time & filtered_data[[indicator_column]] == 1])
    event_count[i] <- length(event_ids)
  }

  # Combine results into a data frame for clarity
  summary_df <- data.frame(
    TimePoint = tp_surv,
    AtRisk = at_risk,
    EventCount = event_count
  )

  # Return the summary data frame
  return(summary_df)
}

# Example usage for 'risk_death_summary' with IndD:
risk_death_summary <- summarize_risk(tp_surv, 'IndD', 'L_open', 'R_closed', filtered_a0)
# Example usage for 'risk_RE_summary' with IndR:
risk_RE_summary <- summarize_risk(tp_end, 'IndR', 'L_open', 'R_closed', filtered_a0)
View(risk_RE_summary)
View(risk_death_summary)



# below is old, manual.
# # Initialize results
# at_risk <- numeric(length(tp_surv))
# death_events <- numeric(length(tp_surv))
#
# # Loop through each time point
# for (i in seq_along(tp_surv)) {
#   time <- tp_surv[i]
#
#   # Find unique IDs of people at risk at this time point
#   at_risk_ids <- unique(filtered_a0$ID[filtered_a0$L_open <= time & filtered_a0$R_closed >= time])
#   at_risk[i] <- length(at_risk_ids)
#
#   # Count the number of death events at this time point
#   death_event_ids <- unique(filtered_a0$ID[filtered_a0$R_closed == time & filtered_a0$IndD == 1])
#   death_events[i] <- length(death_event_ids)
# }
#
# # Combine results into a data frame for clarity
# risk_death_summary <- data.frame(
#   TimePoint = tp_surv,
#   AtRisk = at_risk,
#   DeathEvents = death_events
# )
#
# # Print the results
# # print(head(risk_death_summary,10))
# View(risk_death_summary)


# rightCases_loop with 207 records
rightCases_loop <- c(190, 96, 156, 118, 119, 120, 182, 183, 184, 185, 211, 228,
                     229, 230, 231, 112, 212, 86, 167, 168, 169, 170, 220, 65,
                     25, 26, 157, 110, 94, 45, 87, 195, 196, 197, 198, 201, 175,
                     176, 177, 178, 192, 27, 203, 204, 180, 15, 16, 17, 171, 172,
                     173, 9, 44, 144, 145, 146, 147, 88, 89, 90, 91, 152, 153,
                     154, 113, 54, 55, 56, 57, 28, 46, 47, 70, 95, 117, 121,
                     187, 188, 189, 79, 80, 81, 82, 48, 49, 50, 51, 101, 98, 99,
                     100, 37, 2, 3, 130, 131, 161, 138, 139, 140, 108, 109, 53,
                     134, 135, 136, 174, 29, 30, 77, 78, 67, 68, 69, 19, 20, 21,
                     22, 186, 179, 214, 85, 158, 72, 122, 123, 35, 36, 193, 194,
                     116, 102, 103, 104, 105, 124, 125, 126, 219, 221, 222, 223,
                     224, 155, 92, 93, 1, 141, 38, 148, 149, 150, 151, 165, 11,
                     12, 13, 14, 58, 59, 60, 61, 75, 76, 232, 18, 191, 10, 52,
                     199, 215, 216, 217, 218, 106, 107, 23, 24, 213, 166, 137,
                     4, 5, 6, 7, 200, 74, 39, 40, 41, 42, 225, 226, 227, 62, 63,
                     64, 162, 163, 164, 83, 84, 66, 202, 142, 143, 210)

# rightPeople_loop with size 107
rightPeople_loop <- c(18, 19, 20, 21, 21, 21, 22, 22, 22, 22, 23, 24, 24, 24,
                      24, 25, 26, 27, 28, 28, 28, 28, 29, 30, 31, 31, 32, 33, 34,
                      35, 36, 37, 37, 37, 37, 38, 39, 39, 39, 39, 40, 41, 42, 42,
                      43, 44, 44, 44, 45, 45, 45, 46, 47, 48, 48, 48, 48, 49, 49,
                      49, 49, 50, 50, 50, 51, 52, 52, 52, 52, 53, 54, 54, 55, 56,
                      57, 58, 59, 59, 59, 60, 60, 60, 60, 61, 61, 61, 61, 62, 63,
                      63, 63, 64, 65, 65, 66, 66, 67, 68, 68, 68, 69, 69, 70, 71,
                      71, 71, 72, 73, 73, 74, 74, 75, 75, 75, 76, 76, 76, 76, 77,
                      78, 79, 80, 81, 82, 83, 83, 84, 84, 85, 85, 86, 87, 87, 87,
                      87, 88, 88, 88, 89, 90, 90, 90, 90, 91, 92, 92, 93, 94, 95,
                      96, 96, 96, 96, 97, 98, 98, 98, 98, 99, 99, 99, 99, 100,
                      100, 101, 102, 103, 104, 105, 106, 107, 107, 107, 107, 108,
                      108, 109, 109, 110, 111, 112, 113, 113, 113, 113, 114, 115,
                      116, 116, 116, 116, 117, 117, 117, 118, 118, 118, 119, 119,
                      119, 120, 120, 121, 122, 123, 123, 124)

# rightPeople_loop_og (same as rightCases_loop)
rightPeople_loop_og <- rightCases_loop

