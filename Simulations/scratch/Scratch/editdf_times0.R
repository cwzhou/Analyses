newdat <- df_times0 %>%
  mutate(obs_time = pmin(Time_Failure1, Time_Failure2, Time_Censor, Time_Tau))

# Create the new variable 'Label' based on 'obs_time'
newdat1 <- newdat %>%
  mutate(Label = case_when(
    obs_time == Time_Failure1 ~ "Failure1",
    obs_time == Time_Failure2 ~ "Failure2",
    obs_time == Time_Censor ~ "Censor",
    obs_time == Time_Tau ~ "Tau",
    TRUE ~ NA_character_
  ))
