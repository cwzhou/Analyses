set.seed(2026)
revision = 0
tol_arg = 0.07 #mean_tol1[1]

## ---------------------------------
## Large Monte Carlo population
## ---------------------------------
n_mc <- 5000#100000 # very large sample

# Generate covariates from true DGM
arg.mc = arg_list
arg.mc$ctype = 99 # no censoring
arg.mc$N <- n_mc
arg.mc$u1 = runif(n_mc)
arg.mc$u2 = runif(n_mc)
arg.mc$u3 = runif(n_mc)

arg_mc_action0 <- arg.mc
arg_mc_action0$policy <- "test_action0"  
arg_mc_action1  <- arg.mc        
arg_mc_action1$policy <- "test_action1"    

## ---------------------------------
## True structural outcome model
## Generate data
## ---------------------------------
set.seed(123)
mc.data_action0<- do.call(gdata_CR, arg_mc_action0);
set.seed(123)
mc.data_action1<- do.call(gdata_CR, arg_mc_action1);

data0 <- mc.data_action0
colnames(data0)[colnames(data0) == "failure_t1"] <- "failure_time_cause1_action0"
colnames(data0)[colnames(data0) == "failure_t2"] <- "failure_time_cause2_action0"
colnames(data0)[colnames(data0) == "status"]     <- "status_action0"
colnames(data0)[colnames(data0) == "event.time"] <- "event_time_action0"

data1 <- mc.data_action1
colnames(data1)[colnames(data1) == "failure_t1"] <- "failure_time_cause1_action1"
colnames(data1)[colnames(data1) == "failure_t2"] <- "failure_time_cause2_action1"
colnames(data1)[colnames(data1) == "status"]     <- "status_action1"
colnames(data1)[colnames(data1) == "event.time"] <- "event_time_action1"

# 3. Select only the outcome columns from data1 to avoid duplicate subj.id, Z1, Z2
data1_out <- data1[, c("event_time_action1",
                               "failure_time_cause1_action1",
                               "failure_time_cause2_action1",
                               "status_action1")]

# 4. Merge side by side
library(dplyr)
mc.data <- bind_cols(
  data0,
  data1 %>%
    select(event_time_action1, status_action1,
           failure_time_cause1_action1, failure_time_cause2_action1)
) %>%
  select(subj.id, Z1, Z2,
         event_time_action0, event_time_action1,
         status_action0, status_action1,
         failure_time_cause1_action0, failure_time_cause2_action0,
         failure_time_cause1_action1, failure_time_cause2_action1)
View(mc.data)

mc.data1 <- mc.data %>%
  mutate(
    max_event_time = pmax(event_time_action0, event_time_action1),
    tolerated_max_event_time =  max_event_time * (1 - tol_arg),
    
    # indicator: check non-max elements
    adv_to_ph2 = if_else(
      # for each row, compare the element that is NOT the max
      (event_time_action0 != pmax(event_time_action0, event_time_action1) & 
         event_time_action0 >= tolerated_max_event_time) |
        (event_time_action1 != pmax(event_time_action0, event_time_action1) & 
           event_time_action1 >= tolerated_max_event_time),
      1, 0
    )
  )

mc.data1 %>% dplyr::select(subj.id, max_event_time, tolerated_max_event_time, 
                     event_time_action0, event_time_action1, 
                     tolerated_max_event_time, adv_to_ph2) %>%
  head(2)


mean(mc.data1$adv_to_ph2)

mc.data2 <- mc.data1 %>%
  mutate(final_trt = case_when(
    adv_to_ph2 == 0 & event_time_action0 == max_event_time ~ 0,
    adv_to_ph2 == 0 & event_time_action1 == max_event_time ~ 1,
    adv_to_ph2 == 1 ~ NA_real_
  ))

# mc.data2 %>% dplyr::select(subj.id, max_event_time, tolerated_max_event_time, 
#                            event_time_action0, event_time_action1, 
#                            tolerated_max_event_time, adv_to_ph2,
#                            final_trt) %>%
#   View()

mc.data3 <- mc.data2 %>%
  # mutate(
  #   auc_cif_action0 = ifelse(
  #     adv_to_ph2 == 1,
  #     ifelse(
  #       failure_time_cause1_action0 <= failure_time_cause2_action0,
  #       failure_time_cause1_action0,  # cause 1 occurred first
  #       0                             # cause 2 occurred first
  #     ),
  #     NA_real_  # set NA if adv_to_ph2 != 1
  #   ),
  #   auc_cif_action1 = ifelse(
  #     adv_to_ph2 == 1,
  #     ifelse(
  #       failure_time_cause1_action1 <= failure_time_cause2_action1,
  #       failure_time_cause1_action1,  # cause 1 occurred first
  #       0                             # cause 2 occurred first
  #     ),
  #     NA_real_  # set NA if adv_to_ph2 != 1
  #   )
  # ) 
mutate(
    # auc_cif_action0 = 
    #   ifelse(
    #     failure_time_cause1_action0 <= failure_time_cause2_action0,
    #     failure_time_cause1_action0,  # cause 1 occurred first
    #     0                             # cause 2 occurred first
    #   ),
  auc_cif_action0 = ifelse(
    failure_time_cause1_action0 <= failure_time_cause2_action0 &
      failure_time_cause1_action0 <= arg.mc$tau,
    arg.mc$tau - failure_time_cause1_action0,
    0
  ),
    # auc_cif_action1 = 
    #   ifelse(
    #     failure_time_cause1_action1 <= failure_time_cause2_action1,
    #     failure_time_cause1_action1,  # cause 1 occurred first
    #     0                             # cause 2 occurred first
    #   )
  auc_cif_action1 = ifelse(
    failure_time_cause1_action1 <= failure_time_cause2_action1 &
      failure_time_cause1_action1 <= arg.mc$tau,
    arg.mc$tau - failure_time_cause1_action1,
    0
  )
  ) %>%
  mutate(
    min_cif_time = pmin(auc_cif_action0, auc_cif_action1)
  )

mc.data4 <- mc.data3 %>%
  mutate(
    final_trt = if_else(
      adv_to_ph2 == 1,
      if_else(
        auc_cif_action0 <= auc_cif_action1,  # pick action with smaller CIF AUC
        0,
        1
      ),
      final_trt  # keep existing value for adv_to_ph2 != 1
    )
  )

mc.data5 <- mc.data4 %>%
  mutate(
    event_time_final  = if_else(final_trt == 0, event_time_action0, event_time_action1),
    cif_time_final  = if_else(final_trt == 0, auc_cif_action0, auc_cif_action1),
    status_final      = if_else(final_trt == 0, status_action0, status_action1),
    failure_time_cause1_final = if_else(final_trt == 0, failure_time_cause1_action0, failure_time_cause1_action1),
    failure_time_cause2_final = if_else(final_trt == 0, failure_time_cause2_action0, failure_time_cause2_action1)
  )

mc.data_final <- mc.data5 %>%
  select(subj.id, Z1, Z2, final_trt, adv_to_ph2,
         event_time_final, cif_time_final, status_final,
         failure_time_cause1_final, failure_time_cause2_final)

mean(mc.data_final$event_time_final)*365.25
mean(mc.data_final$cif_time_final)*365.25
