library(survival)
View(bladder1)
# https://vincentarelbundock.github.io/Rdatasets/doc/survival/bladder.html
# ?bladder1
# This is the well-known bladder tumor clinical trial conducted by 
# the Veterans Administration Co-operative Urological Research Group (Byar, 1980). 
# A total of 116 patients with superficial bladder tumors were randomly
# assigned to placebo, pyridoxine (vitamin BG), or thiotepa. 
# The goal of the study was to determine the effectiveness of pyridoxine
# and thiotepa in reducing the rate of tumor recurrence.
# By the end of follow-up, there were 87, 57, and 45 recurrences
# among the 47, 31, and 38 patients in the placebo, pyridoxine,
# and thiotepa groups, respectively. Some patients died before
# any tumor recurrence, while others died after recurrences.
# There were 10, 7, and 11 observed deaths in the placebo, pyridoxine,
# and thiotepa groups, respectively.

library(tidyverse)
bladder1 %>% filter(stop == 0)
blad <- bladder1[!bladder1$id %in% c(1, 49), ] %>%
  dplyr::select(id, start, stop, status, treatment, number, size, recur, enum) %>%
  mutate(IndR = ifelse(status == 1, 1, 0),
         IndD = ifelse(status == 2 | status == 3, 1, 0))
# View(blad)

blad %>% 
  group_by(treatment) %>% 
  summarise(
    death = sum(IndD, na.rm = TRUE),  # Sum of deaths
    recurr = sum(IndR, na.rm = TRUE), # Sum of recurrences
    patients = n_distinct(id)        # Count of unique patients
  )

blad1 = blad %>%
  mutate(Trt = ifelse(treatment == "placebo", 0,
                      ifelse(treatment == "pyridoxine", 2,
                             ifelse(treatment == "thiotepa", 1, NA))),
         L_open = start,
         R_closed = stop) %>%
  dplyr::select(-c(treatment, status, start, stop))
# View(blad1)
library(itrSurv)

