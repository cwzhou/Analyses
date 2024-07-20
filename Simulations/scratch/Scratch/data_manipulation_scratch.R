# Set the seed for reproducibility
set.seed(123)

# Number of subjects
N <- 3

# Sample dataset
df <- data.frame(ID = 1:N,
                 Time_Recurrent1 = runif(N)) %>%
  mutate(Time_Recurrent2 = Time_Recurrent1 + runif(N, min=0, max = 0.2),
                 Time_Recurrent3 = runif(N),
                 # Time_Recurrent4 = runif(N),
                 # Time_Recurrent5 = runif(N),
                 # Time_Recurrent6 = runif(N),
                 # Time_Recurrent7 = runif(N),
                 # Time_Recurrent8 = runif(N),
                 # Time_Recurrent9 = runif(N),
                 # Time_Recurrent10 = runif(N),
                 Time_Failure = runif(N)+0.7,
                 Time_Censor = runif(N)+0.3,
                 Time_Tau = tau) # %>%
  # group_by(ID) %>%
  # mutate(Time_minTC = min(Time_Failure, Time_Censor)) %>%
  # ungroup()

# turn into long-format
df_long = df %>%
  pivot_longer(cols = starts_with("Time_"), names_to = "Label", values_to = "Time") %>%
  mutate(Label = gsub("Time_", "", Label)) %>%
  group_by(ID) %>%
  arrange(ID, Time)

# find min(failure, censoring, tau)
df_min = df_long %>% 
  group_by(ID) %>%
  mutate(failure_status = ifelse(Label == "Failure", 1, ifelse(Label == "Censor" | Label == "Tau", 0, 99))) %>%
  filter(Label %in% c("Censor", "Failure", "Tau")) %>% 
  summarise(min = min(Time),
            failure_status_raw = failure_status[which.min(Time)],
            Label = Label[which.min(Time)])

# combine with long dataset a
tmp1 = inner_join(df_long, df_min, by = "ID") %>% 
  # only keep observed data
  filter(Time <= min) %>% 
  as.data.frame() %>% 
  dplyr::select(-Label.y) %>%
  group_by(ID) %>%
  # create dummy var for first row within ID, recurrence, failure, censoring
  mutate(FirstRow = ifelse(row_number() == 1, 1, 0),
         contains_recurrent = (grepl("recurrent", Label.x, ignore.case = TRUE)),
         contains_failure = (grepl("failure", Label.x, ignore.case = TRUE)),
         contains_censor = (grepl("censor", Label.x, ignore.case = TRUE))) %>%
  # use lag to move future times to "L_open"
  mutate(L_open = ifelse(FirstRow ==1, 0, lag(Time)), #open paranthases ("(")
         R_closed = Time, #closed brack ("]")
         IndR = ifelse(contains_recurrent == TRUE, 1, 0),
         IndD = ifelse(contains_failure == TRUE, 1, 0)) %>%
  rename(R_Label = Label.x) %>%
  dplyr::select(ID, L_open, R_closed, IndR, IndD, R_Label); tmp1
  # dplyr::select(-c(Time, min, FirstRow, failure_status_raw, contains_recurrent, contains_failure, contains_censor)); tmp2

# # One-hot encode the Color variable using model.matrix
# df_encoded <- data.frame(model.matrix(~ Label - 1, data = df_long))
# colnames(df_encoded) <- gsub("Label", "Ind", colnames(df_encoded))
# df_full = cbind(df_long, df_encoded);print(df_full)



