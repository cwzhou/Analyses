library(tidycmprsk)
library(tidyverse)
library(ggsurvfit)
data(Melanoma, package = "MASS")
Melanoma <- 
  Melanoma %>% 
  mutate(
    status = as.factor(recode(status, `2` = 0, `1` = 1, `3` = 2))
  )
head(Melanoma)

cuminc(Surv(time, status) ~ 1, data = Melanoma) %>% 
  ggcuminc(outcome = c("1", "2")) +
  ylim(c(0, 1)) + 
  labs(x = "Days")

cuminc(Surv(time, status) ~ ulcer, data = Melanoma) %>% 
  tbl_cuminc(
    times = 1826.25, 
    label_header = "**{time/365.25}-year cuminc**") %>% 
  add_p()

cuminc(Surv(time, status) ~ ulcer, data = Melanoma) %>% 
  ggcuminc() + 
  labs(
    x = "Days"
  ) + 
  add_confidence_interval() +
  add_risktable()

#########################
head(data.df)
cuminc(Surv(obs_time, factor(status)) ~ 1, 
       data = data.df) %>% 
  ggcuminc(outcome = c("1", "2")) +
  ylim(c(0, 1)) + 
  labs(x = "Years")

cuminc(Surv(obs_time, factor(status)) ~ Trt, 
       data = data.df) %>% 
  tbl_cuminc(
    times = 1, 
    label_header = "**{time}-year cuminc**") %>% 
  add_p()

cuminc(Surv(obs_time, factor(status)) ~ Trt, 
       data = data.df) %>% 
  ggcuminc(outcome = c("1", "2")) + 
  labs(x = "Years",
       color = "Trt") #+ 
  # add_confidence_interval() #+
  # add_risktable()

cuminc(Surv(event.time, factor(status)) ~ action, 
       data = rep_czmk) %>% 
  ggcuminc(outcome = c("1", "2")) + 
  labs(x = "Years",
       color = "Trt") #+ 
# add_confidence_interval() #+
# add_risktable()
