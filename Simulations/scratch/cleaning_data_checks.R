proc3 <- proc2 %>%
  mutate(open_date = case_when(
    procGroup == "Open" ~ PROCEDURE_DATE,
    TRUE ~ NA
  ),
  endovascular_date = case_when(
    procGroup == "Endovascular" ~ PROCEDURE_DATE,
    TRUE ~ NA
  ),
  major_amputation_date = case_when(
    procGroup == "Major amputation" ~ PROCEDURE_DATE,
    TRUE ~ NA
  )
  )

proc3 %>%
  group_by(CURRENT_MRN) %>%
  summarize_all(.funs = list(~first(na.omit(.)))) %>%
  View()


dups_death = duplicated(death2$CURRENT_MRN); sum(dups_death)
dups_demo = duplicated(demo2$CURRENT_MRN); sum(dups_demo)
dups_comor = duplicated(comor0$CURRENT_MRN); sum(dups_comor)

comor0[dups_comor,] %>% arrange(CURRENT_MRN) %>% View()
multiple_rows_ids <- unique(death2$CURRENT_MRN[duplicated(death2$CURRENT_MRN)])

# Subset the original proc2 frame based on IDs with multiple rows
subset_df <- df[df$ID %in% multiple_rows_ids, ]




t1 <- demo4 %>% arrange(CURRENT_MRN)
t2 <- death3 %>% arrange(CURRENT_MRN)
# Merge or join the data frames based on the ID variable
merged_data <- merge(t1, t2, by = "CURRENT_MRN", suffixes = c("_t1", "_t2")) %>%
  dplyr::select(-c(GENDER, ETHNICITY, PATIENT_RACE_1))
# Check if birth dates are the same for corresponding IDs
same_birth_dates <- merged_data$BIRTH_DATE_t1 == merged_data$BIRTH_DATE_t2
same_birth_dates
# Filter the merged dataset to include only individuals with different birth dates
different_birth_dates <- merged_data[!same_birth_dates, ]
different_birth_dates

# Perform anti-join to find rows in t2 that didn't match with t1 based on CURRENT_MRN
rows_not_in_merged_data <- anti_join(t2, merged_data, by = "CURRENT_MRN")
# Display the rows from t2 that didn't make it into merged_data
rows_not_in_merged_data


same_birth_dates <- (merged_data$EPIC_DEATH_DATE_t1 == merged_data$EPIC_DEATH_DATE_t2) | (is.na(merged_data$EPIC_DEATH_DATE_t1) & is.na(merged_data$EPIC_DEATH_DATE_t2))
same_birth_dates
# Filter the merged dataset to include only individuals with different birth dates
different_birth_dates <- merged_data[!same_birth_dates, ]
different_birth_dates

# Perform anti-join to find rows in t2 that didn't match with t1 based on CURRENT_MRN
rows_not_in_merged_data <- anti_join(t2, merged_data, by = "CURRENT_MRN")
# Display the rows from t2 that didn't make it into merged_data
rows_not_in_merged_data


same_birth_dates <- (merged_data$STATE_DEATH_DATA_DATE_t1 == merged_data$STATE_DEATH_DATA_DATE_t2) | (is.na(merged_data$STATE_DEATH_DATA_DATE_t1) & is.na(merged_data$STATE_DEATH_DATA_DATE_t2))
same_birth_dates
# Filter the merged dataset to include only individuals with different birth dates
different_birth_dates <- merged_data[!same_birth_dates, ]
different_birth_dates

# Perform anti-join to find rows in t2 that didn't match with t1 based on CURRENT_MRN
rows_not_in_merged_data <- anti_join(t2, merged_data, by = "CURRENT_MRN")
# Display the rows from t2 that didn't make it into merged_data
rows_not_in_merged_data


merged_data %>% dplyr::select(-c(BIRTH_DATE_t1, BIRTH_DATE_t2)) %>%
  View()

dat3$BIRTH_DATE.x==dat3$BIRTH_DATE.y
sum(is.na(dat3$BIRTH_DATE.x)) #0
sum(is.na(dat3$BIRTH_DATE.y)) #5