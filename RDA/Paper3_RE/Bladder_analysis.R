library(purrr);
library(dplyr);
library(survival)
library(randomForestSRC);
library(caret);
library(itrSurv);library(ggplot2); library(tidyverse);
library(MASS);
local = 1
if (local == 1){
  setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper3_RE/")
  source("02.Simulation_Functions_RE.R") # includes F02.ComparatorMethod_Functions.R
  setwd("../../RDA/Paper3_RE/")
} else{
  setwd("/nas/longleaf/home/cwzhou/Dissertation/Analyses/Simulations/Paper3_RE")
  source("02.Simulation_Functions_RE.R") # includes F02.ComparatorMethod_Functions.R
  setwd("/nas/longleaf/home/cwzhou/Dissertation/Analyses/RDA/Paper3_RE")
}
source("RDA_Functions_RE.R")    # All the library and source files
############################################################
######### temporary testing things out parameters ##########
############################################################
criterion_phase1 = "mean"
rule2 = "gray_re"
############################################################
################# methods for comparison ###################
############################################################
rda_methods = c("CZMK", "ZOM", "OBS")
skip_method <- c(!TRUE,!TRUE,!TRUE);
assign_skip_function(rda_methods, skip_method)
skipped_methods <- rda_methods[skip_method]
loop_methods <- rda_methods[!rda_methods %in% c(skipped_methods, "OBS")]
############################################################
#################### 0. parameters ####################
############################################################
tol1_param = c(0.05,0) # c(0.03,0)
t0_crit = 30 #2200 # 6-mo survival
pooled1 = FALSE # stratified = lower nodesize; pooled = can have larger nodesize
tau = 54
K = 5 # number of CV
endpoint = "RE" # endpoint
Tx.nm = "A"
init_seed = 116

## other parameters
nodeSizeSurv = 3
nodeSizeEnd = 3
minEventSurv = round(sqrt(c(nodeSizeSurv)), 0)
minEventEnd = round(sqrt(c(nodeSizeEnd)), 0)
Ntree = 300
ert = FALSE
rs = 0.2 # randomSplit = 0.2

############################################################
##################### Bladder Dataset ######################
############################################################
library(tidyverse)
# # bladder2 dataset generated from SAS code #
# file_path <- "~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Scratch/Project3/bladder2.csv"
# # Read the CSV file into a data frame, including headers
# bladder2 <- read.csv(file_path, header = TRUE)
# colnames(bladder2) = c("id", "TStart", "TStop", "A", "Z1", "Z2", "Visit", "Status", "Gaptime")
# # this dataset doesnt have terminal events so just
# # BS-ing a few in for now to have that variable as we
# # test the package/code up the method
# # Step 1: Create a new variable status_D that is initially a copy of status
# bladder2$Status_D <- bladder2$Status
# # Step 2: For rows where status is 0, randomly set 90% to 1 in status_D
# # Create a logical vector to identify rows where status is 0
# zero_status_rows <- bladder2$Status == 0
# # Randomly sample 90% of these rows and set them to 1 in status_D
# set.seed(123) # Set seed for reproducibility, if needed
# num_to_change <- sum(zero_status_rows) * 0.9
# change_indices <- sample(which(zero_status_rows), num_to_change)
# bladder2$Status_D[change_indices] <- 1
# # Step 3: Set Status_D to 0 where Status is 1
# # Status = 1 for recurrent event and = 0 for censored or death
# # Status_D = 1 for death and = 0 otherwise
# bladder2$Status_D[bladder2$Status == 1] <- 0
# # Check the result
# head(bladder2)
# make it so people who die at time TStart = TStop have an interval where TStart = TStop + 0.001
# bladder2.1 <- bladder2 %>%
#   mutate(
#     TStop = ifelse(TStart == TStop, TStop + 0.001, TStop)  # Avoids TStart = TStop
#   )
# # need to add in censoring rows for the last dataset if status_d = 0
# df_tmp = bladder2.1 %>%
#   arrange(id, Visit) %>%
#   group_by(id) %>%
#   mutate(
#     need_to_add_row = ifelse(!(row_number() == n() & Status_D == 0 & Visit > 1), FALSE,
#                              ifelse(Status == 0 & Status_D == 0, FALSE,
#                                     TRUE))) %>%
#   ungroup() %>%
#   mutate(
#     max_stop_time = max(TStop)
#   )
# new_rows <- df_tmp %>%
#   filter(need_to_add_row) %>%
#   mutate(
#     TStart1 = TStop,
#     TStop1= max_stop_time,
#     Visit1 = Visit + 1,
#     Gaptime1 = max_stop_time - TStop,
#     Status1 = 0
#   ) %>%
#   dplyr::select(id, TStart = TStart1, TStop = TStop1, A, Z1, Z2, Visit=Visit1,
#                 Status = Status1, Gaptime=Gaptime1, Status_D)
# new_rows
# # Combine the original data with the new rows
# df <- bind_rows(bladder2.1, new_rows) %>%
#   group_by(id) %>%
#   arrange(id, Visit) %>%
#   ungroup() %>%
#   as.data.frame()
# df

library(survival)
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
  mutate(Status = ifelse(status == 1, 1, 0),
         Status_D = ifelse(status == 2 | status == 3, 1, 0))
blad %>%
  group_by(treatment) %>%
  summarise(
    death = sum(Status_D, na.rm = TRUE),  # Sum of deaths
    recurr = sum(Status, na.rm = TRUE), # Sum of recurrences
    patients = n_distinct(id)        # Count of unique patients
  )

bladder2.1 = blad %>%
  mutate(A = ifelse(treatment == "placebo", 0,
                      ifelse(treatment == "pyridoxine", 2,
                             ifelse(treatment == "thiotepa", 1, NA))),
         TStart = start,
         TStop = stop) %>%
  dplyr::select(-c(treatment, status, start, stop))
bladder2.1 %>% filter(TStart == TStop)

# need to add in censoring rows for the last dataset if status_d = 0
df_tmp = bladder2.1 %>%
  arrange(id, enum) %>%
  group_by(id) %>%
  dplyr::mutate(
    need_to_add_row = ifelse(!(row_number() == n() & Status_D == 0 & enum > 1), FALSE,
                             ifelse(Status == 0 & Status_D == 0, FALSE,
                                    TRUE))) %>%
  ungroup() %>%
  mutate(
    max_stop_time = max(TStop)
  )

df_tmp1 <- df_tmp %>%
  mutate(
    next_id = lead(id), # Get the next row's id
    need_to_add_row = ifelse(enum == 1 & Status == 1 & id != next_id, TRUE, need_to_add_row)
  ) %>%
  dplyr::select(-next_id) # Remove intermediate column if not needed
# View(df_tmp1)

new_rows <- df_tmp1 %>%
  filter(need_to_add_row) %>%
  mutate(
    TStart1 = TStop,
    TStop1= max_stop_time,
    Visit1 = enum + 1,
    Gaptime1 = max_stop_time - TStop,
    Status1 = 0
  ) %>%
  mutate(Visit = Visit1,
         Status = Status1,
         Gaptime = Gaptime1,
         TStart = TStart1,
         TStop = TStop1) %>%
  dplyr::select(id, TStart, TStop, A, number, size, recur,
                Visit, Gaptime, Status, Status_D)
head(new_rows)
# Combine the original data with the new rows
df <- bind_rows(bladder2.1, new_rows) %>%
  group_by(id) %>%
  arrange(id, enum) %>%
  ungroup() %>%
  as.data.frame()
# head(df)

bladder_surv = df %>%
  filter(Status == 0); #View(bladder_surv)
#checking no. subjects
# bladder2.1 %>% as.data.frame() %>% dplyr::select(id) %>% distinct() %>% nrow()
# bladder_surv %>% as.data.frame() %>% dplyr::select(id) %>% distinct() %>% nrow()
# df %>% as.data.frame() %>% dplyr::select(id) %>% distinct() %>% nrow()

dat1 = bladder_surv # Phase 1 Surv
dat2 = df # Phase 2 Endpoint (RE)

bladder_df = dat2
dataset_name = "bladder"

bladder_df$id %>% unique() %>% length() # 116 patients
sub_ids = bladder_df$id %>% unique()
num_subj = length(sub_ids)

data_to_use = bladder_df %>%
  dplyr::select(-c(enum, Visit, Gaptime))
a1 = data_to_use %>% filter(A == 1) ;
end_a1 = a1 %>%
  dplyr::select(id, number, size, Status, Status_D, TStop)
  # dplyr::select(id, Z1, Z2, Status, Status_D, TStop)
surv_a1 = a1 %>%
  filter(Status == 0) %>%
  dplyr::select(id, number, size, Status, Status_D)
  # dplyr::select(id, Z1, Z2, Status, Status_D)
recurr_a1 = a1 %>%
  filter(Status == 1) %>%
  dplyr::select(id, number, size, Status, Status_D)
  # dplyr::select(id, Z1, Z2, Status, Status_D)
num_unique_people <- n_distinct(recurr_a1$id)

timePointsSurvival = data_to_use %>%
  filter(Status_D == 1) %>%
  dplyr::select(TStop) %>%
  distinct() %>%
  unlist(use.names = F)

timePointsEndpoint = data_to_use %>%
  filter(Status == 1) %>%
  dplyr::select(TStop) %>%
  distinct() %>%
  unlist(use.names = F)

# Get all column names of the dataframe
all_column_names <- colnames(data_to_use)
# Define the column names you want to exclude
excluded_column_names <- c("id", "A", "TStart", "TStop",
                           "Status_D", "Status", "recur")
# Exclude the specified column names
covariate_names <- setdiff(all_column_names, excluded_column_names)

# manipulating dataset to obtain truncated status,D.0,D.1,D.2
dat = data_to_use %>%
  dplyr::select(id, TStart, TStop, A,
                Status_D, Status,
                recur,
                !!!covariate_names) %>%
  as.data.frame()

#################### 1. data preprocessing
# mean(dat$Status_D==1, na.rm = TRUE)
# mean(dat$Status==1, na.rm = TRUE) # 62%
lvls = lapply(Tx.nm, function(x) levels(dat[, x]))
for (x in 1:length(lvls)) {if (is.null(lvls[[x]])) lvls[[x]] = 0:2}

############################################################
###################### 0.2 criterion #######################
############################################################
print(criterion_phase1)
criterion_phase2 = criterion_phase1
if (!(criterion_phase1[1] %in% c("mean", "area"))) {
  criterion_phase1[2] = t0_crit
  rule1 = "logrank_surv"
} else {
  rule1 = "mean_surv"
  criterion_phase1[2] = NA
}

if (!(criterion_phase2[1] %in% c("mean", "area"))) {
  criterion_phase2[2] = t0_crit
} else {
  criterion_phase2[2] = NA
}
if (is.null(criterion_phase1[2]) | is.na(criterion_phase1[2])){
  crit_tmp = ""
} else{
  crit_tmp = round(as.numeric(criterion_phase1[2]), 1)
}

cat("\ncriterion_phase1 =",criterion_phase1[1],
    "\nvalue.os =", criterion_phase1[2],
    "\nvalue.re =", criterion_phase2[2],
    "\ntau =", tau, "\n")

#################### propensity models
modelPr = create_model_formulas("factor(A)",
                                covariate_names,
                                formula = FALSE)
form.weight <- list(modelPr %>% as.formula)

# skeleton
# Exclude "observed" for prop2_values
traintest_methods <- setdiff(rda_methods, "OBS")
prop2_values = c(".PropPhase2")
# train_values = c(".train.terminal", ".train.RE")
test_values = c(".terminal", ".RE")
values_colsnames <- c(
  # unlist(lapply(rda_methods, function(method) paste0(method, train_values))),
  unlist(lapply(rda_methods, function(method) paste0(method, test_values))),
  unlist(lapply(traintest_methods, function(method) paste0(method, prop2_values))),
  "ns1.CZMK",
  "ns2.CZMK",
  "train_cens", "test_cens"
) %>%
  unlist()
# Define a custom sorting function
custom_sort <- function(names) {
  order(
    !grepl("train.terminal", names),
    !grepl("train.RE", names),
    !grepl("terminal", names),
    !grepl("RE", names),
    !grepl("PropPhase2", names),
    !grepl("train_cens", names),
    !grepl("test_cens", names),
    !grepl("^time_", names),
    !grepl("^CZMK_", names),
    !grepl("^ZOM_", names),
    !grepl("^OBS_", names),
    names
  )
}
sorted_values_colsnames <- values_colsnames[custom_sort(c(values_colsnames))]
values <- matrix(NA, K, length(values_colsnames),
                 dimnames = list(1:K, values_colsnames))
train_eval_result_list = train_eval_result_list_tmp = list()

nm = paste0(criterion_phase1[1],
            crit_tmp, "_rule1",
            rule1, "_rule2", rule2, "_tau", tau)
# if (local == 1){
#   fnm = paste0("./output/", endpoint, "/",dataset_name,"/", Sys.Date(), "/") # folder name
# } else{
#   fnm = paste0("/work/users/c/w/cwzhou/Proj3RE/output/", endpoint, "/",dataset_name,"/", Sys.Date(), "/")
# }
# if (!dir.exists(sprintf("./output/%s", endpoint))) dir.create(sprintf("./output/%s", endpoint))
# if (!dir.exists(sprintf("./output/%s/%s", endpoint, dataset_name))) dir.create(sprintf("./output/%s/%s", endpoint, dataset_name))
# if (!dir.exists(sprintf("./figure/%s", endpoint))) dir.create(sprintf("./figure/%s", endpoint))
# if (!dir.exists(sprintf("./figure/%s/%s", endpoint, dataset_name))) dir.create(sprintf("./figure/%s/%s", endpoint, dataset_name))
# if (!dir.exists(fnm)) dir.create(fnm)
if (local == 1) {
  base_dir <- "./output/"
  figure_dir <- "./figure/"
  fnm <- paste0(base_dir, endpoint, "/", dataset_name, "/", Sys.Date(), "/") # folder name
} else {
  base_dir <- "/work/users/c/w/cwzhou/Proj3RE/output/"
  figure_dir <- "/work/users/c/w/cwzhou/Proj3RE/figure/"
  fnm <- paste0(base_dir, endpoint, "/", dataset_name, "/", Sys.Date(), "/") # folder name
}

# Create necessary directories
if (!dir.exists(file.path(base_dir, endpoint))) {
  dir.create(file.path(base_dir, endpoint), recursive = TRUE)
}
if (!dir.exists(file.path(base_dir, endpoint, dataset_name))) {
  dir.create(file.path(base_dir, endpoint, dataset_name), recursive = TRUE)
}
if (!dir.exists(file.path(figure_dir, endpoint))) {
  dir.create(file.path(figure_dir, endpoint), recursive = TRUE)
}
if (!dir.exists(file.path(figure_dir, endpoint, dataset_name))) {
  dir.create(file.path(figure_dir, endpoint, dataset_name), recursive = TRUE)
}
if (!dir.exists(fnm)) {
  dir.create(fnm, recursive = TRUE)
}
rds = paste0("_", Sys.Date(), ".rds")

tab = data.frame(s1 = rep(NA, K),
                 # s2 = NA,
                 total = NA)

############################################################
####################### 0.X weights ########################
############################################################
# weights
weights_bladder <- function(data, lvls1) {
  # treatment was randomized
  s.dat = data %>% filter(Status == 0)
  propensity1 = sapply(lvls1, function(x) mean(s.dat$A == x))
  propensity1 = sapply(s.dat$A, function(s) propensity1[lvls1 == s])
  propensity  = propensity1

  # IPCW
  Sc.hat1 <- coxph(Surv(TStop, 1 - Status_D) ~ number + size, data = s.dat)
  Sc.hat1 <- exp( - predict(Sc.hat1, type = "expected"))
  Sc.hat <- Sc.hat1
  weight.censor = s.dat$Status_D/Sc.hat

  return(list(propensity = propensity, weight.censor = weight.censor))
}
weight0 = weights_bladder(data = dat, lvls1 = lvls[[1]])

############################################################
########################## 0.2 cv ##########################
############################################################
## cross-validation
K = K
dat_surv = dat %>% filter(Status == 0); head(dat_surv)
# Create cross-validation folds for individuals
set.seed(init_seed)
if (K > 1){
  folds <- createFolds(factor(dat_surv$Status_D), k = K, list = TRUE)
} else{
  folds <- createFolds(factor(dat_surv$Status_D), k = 2, list = TRUE)
}
# library(plyr)
# d1 = dat_surv %>% dplyr::select(id, Status_D)
# # Loop through each fold and compute the proportion of Status_D for each fold
# fold_prop_res <- lapply(1:K, function(k) {
#   d1$fold <- ifelse(d1$id %in% sub_ids[folds[[k]]], k, NA)  # Assign fold number
#   # Summarize the proportion of Status_D for each fold
#   fold_summary <- ddply(d1[!is.na(d1$fold), ], 'fold', summarise, prop = mean(Status_D))
#   return(fold_summary)
# })
# # Combine results from all folds
# fold_prop <- do.call(rbind, fold_prop_res); print(fold_prop)

d1 <- dat_surv %>%
  dplyr::select(id, Status_D)  # Select necessary columns
# Loop through each fold and compute the proportion of Status_D for each fold
fold_prop_res <- lapply(1:K, function(k) {
  # Assign fold number to d1
  d1$fold <- ifelse(d1$id %in% sub_ids[folds[[k]]], k, NA)
  # Summarize the proportion of Status_D for each fold using dplyr
  fold_summary <- d1 %>%
    filter(!is.na(fold)) %>%
    group_by(fold) %>%
    summarise(prop = mean(Status_D)) %>%
    ungroup()  # Ensure to ungroup after summarizing

  return(fold_summary)
})
fold_prop <- bind_rows(fold_prop_res);print(fold_prop)


############################################################
####################### STARTING CV ########################
############################################################
for (cv in 1:K) {
  cat(cv, "th cv.\n")
  set.seed(cv)

  # Get in-sample (training) and out-sample (testing) IDs
  test_ids <- sub_ids[folds[[cv]]]
  train_ids <- setdiff(sub_ids, test_ids)
  test_indices = folds[[cv]] # outsample person index
  train_indices = setdiff(c(1:num_subj), test_indices) # insample person index

  # Subset data based on IDs
  train <- dat %>% filter(id %in% train_ids) %>%
    dplyr::select(-c(recur)) %>% as.data.frame()
  # View(train)
  test  <- dat %>% filter(id %in% test_ids) %>% as.data.frame()

  # Debug print statements
  print(paste("Train size:", nrow(train), "with Unique IDs in Train:", length(unique(train$id))))
  print(paste("Test size:", nrow(test), "with Unique IDs in Test:", length(unique(test$id))))

  # Calculate the proportion of censored data conditionally on Status == 0
  train_censored <- train %>% filter(Status == 0)
  test_censored <- test %>% filter(Status == 0)

  values[cv, "train_cens"] <- mean(train_censored$Status_D == 0)
  values[cv, "test_cens"] <- mean(test_censored$Status_D == 0)

  # model1 = "Surv(TStop, Status_D) ~ Z1 + Z2" %>% as.formula()
  # model2 = "Surv(TStart, TStop, Status) ~ Z1 + Z2" %>% as.formula()
  model1 = "Surv(TStop, Status_D) ~ number + size" %>% as.formula()
  model2 = "Surv(TStart, TStop, Status) ~ number + size" %>% as.formula()
  models_RE = list(model1, model2)

  data_surv = train %>% filter(Status == 0)
  missing_in_data_surv <- data_surv$id[!data_surv$id %in% train$id]
  missing_in_train <- train$id[!train$id %in% data_surv$id]
  print(missing_in_data_surv) # IDs in data_surv missing from train
  print(missing_in_train)     # IDs in train missing from data_surv

  # observed policy: training
  # values[cv, "OBS.train.RE"] = mean(train$Status==1, na.rm = T)
  # values[cv, "OBS.train.terminal"] # need to-do: figure out this one

  args.CZMK <- list(data = train,
                    endPoint = "RE",
                    idName = "id",
                    epName = "Status",
                    txName = "A",
                    models = models_RE,
                    timePointsSurvival = timePointsSurvival,
                    timePointsEndpoint = timePointsEndpoint,
                    tau = tau,
                    criticalValue1 = "mean",
                    criticalValue2 = "mean",
                    evalTime = 1,
                    splitRule1 = "mean_surv",
                    splitRule2 = "gray_re",
                    ERT = FALSE,
                    uniformSplit = TRUE,
                    replace = FALSE,
                    randomSplit = 0.2,
                    nTree = Ntree,
                    pooled = FALSE,
                    tol1 = tol1_param,
                    stratifiedSplit = 0.1)

  ############################################################################################################
  #### A. rule estimation ####################################################################################
  ############################################################################################################

  ### A1. The proposed method
  if (skip.CZMK != TRUE){
    print(sprintf("Running CZMK for %s CV", cv))

    values[cv, "ns1.CZMK"] = nodeSizeSurv
    values[cv, "ns2.CZMK"] = nodeSizeEnd
    set.seed(cv)
    CZMK.i <-
      try(do.call(itrSurv::itrSurv,
                  c(args.CZMK, list(nodeSizeSurv = nodeSizeSurv,
                                    nodeSizeEnd = nodeSizeEnd,
                                    minEventSurv = minEventSurv,
                                    minEventEnd = minEventEnd))))
    # values[cv, "CZMK.train.terminal"] = CZMK.i@value[["V1"]][["Et_survival"]]
    values[cv, "CZMK.PropPhase2"] = CZMK.i@value[["V2"]][["PropPhase2"]]
    # values[cv, "CZMK.train.RE"] = CZMK.i@value[["V3"]][["Et_mff"]]
    err.CZMK = class(CZMK.i)[1] == "try-error"
  } else{
    err.CZMK = TRUE
  }

  ### A6. zero-order model
  if (skip.ZOM != TRUE){
    print(sprintf("Running ZOM for %s CV", cv))
    set.seed(cv)
    ZOM.i <-
      try(do.call(itrSurv::itrSurv, c(args.CZMK,
                                      list(nodeSizeEnd = 1e+4,
                                           nodeSizeSurv = 1e+4,
                                           minEventEnd = 1e+4,
                                           minEventSurv = 1e+4))))
    # values[cv, "ZOM.train.terminal"] = ZOM.i@value[["V1"]][["Et_survival"]]
    values[cv, "ZOM.PropPhase2"] = ZOM.i@value[["V2"]][["PropPhase2"]]
    # values[cv, "ZOM.train.RE"] = ZOM.i@value[["V3"]][["Et_mff"]]
    err.ZOM = class(ZOM.i)[1] == "try-error"

    if (!err.ZOM) {
      zom.pred =
        data.frame(
          ZOM.i@phaseResults$FinalOptimalTx_Recc[1])
      tab[cv, 1] = ZOM.i@phaseResults$FinalOptimalTx_Recc[1]
      # cat("ZOM freq table\n"); tab[, 1] %>% table %>% print
    }
  } else{
    err.ZOM = TRUE
  }

  message("Applying rules to test set.")
  ############################################################################################################
  #### B. Rules applied to a test set ########################################################################
  ############################################################################################################
  ### B0. skeletons for the predicted Trt's
  opt.rule.pred = data.frame(Trt = rep(NA, length(unique(test$id))))
  opt.rule.pred[,1] = factor(NA, levels = lvls[[1]])

  for(method in loop_methods) {
    assign(paste("opt.rule.", method, sep = ""),
           opt.rule.pred,
           envir = .GlobalEnv)
  }
  # Create eligibility for each row based on non-missing values in `Tx.nm`
  elig0 <- test %>%
    dplyr::select(id, !!sym(Tx.nm)) %>%
    group_by(id) %>%
    dplyr::mutate(eligibility = all(!is.na(.data[[Tx.nm]]))) %>% # Check if all values are TRUE
    dplyr::mutate(eligibility = ifelse(any(!eligibility), FALSE, eligibility)) %>% # Set all to FALSE if any are FALSE
    ungroup()
  elig = elig0 %>%
    pull(eligibility)
  id_elig = elig0 %>%
    dplyr::select(id, eligibility) %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup() %>%
    pull(eligibility)

  ### B1. the proposed method
  if (!err.CZMK) {
    for (phase in 1:2){
      opt.CZMK_phase =
        itrSurv::predict(CZMK.i,
                newdata = test[elig, ],
                Phase = phase,
                epName1 = "Status",
                endPoint = "RE")
      if (phase == 1){
        action1 <<- opt.CZMK_phase$optimal@optimalTx
        StopatP1 <<- opt.CZMK_phase$optimal@Ratio_Stopping_Ind
      }
      if (phase ==2){
        action2 <<- opt.CZMK_phase$optimal@optimalTx
      }
    }
    tmp_act = as.data.frame(cbind(P1=action1, P2=action2, StopatP1))
    tmp_act = tmp_act %>%
      mutate(Final = ifelse(StopatP1 == 1,
                            P1,
                            P2))
    opt.CZMK = tmp_act$Final
    opt.rule.CZMK[id_elig, ] = factor(opt.CZMK,
                                   levels = lvls[[1]]) # %>% as.numeric()
    print(table(opt.rule.CZMK, useNA = "always"))
  }

  ### B6. zero-order model
  if (!err.ZOM) {
    opt.rule.ZOM[id_elig, ] = rep(zom.pred[1, ], length(id_elig))
    opt.rule.ZOM[id_elig, ] = factor(zom.pred[1, ], levels = lvls[[1]]) # %>% as.numeric()
    print(table(opt.rule.ZOM, useNA = "always"))
  }

  ############################################################################################################
  #### C. Value estimation ###################################################################################
  ############################################################################################################
  test.surv = test %>%
    filter(Status == 0)
  test.tmp = test.surv%>%
    transmute(Trt = test.surv[, Tx.nm])

  if (criterion_phase1[1] != criterion_phase2[1]){
    message("WARNING: CRITERION_PHASE1 AND CRITERION_PHASE2 ARE DIFFERENT - MAKE SURE THAT IS WHAT YOU WANT")
  }
  suffixes = c(".terminal", ".RE")
  for (suffix in 1:length(suffixes)){
    suffix1 = suffixes[suffix]
    if (suffix1 == ".terminal"){
      event_indicator_string = "Status_D"
      criterion = criterion_phase1
      endpoint1 = "survival"
    } else if (suffix1 == ".RE"){
      event_indicator_string = "Status"
      criterion = criterion_phase2
      endpoint1 = "mff"
    }
    testing_dataset = test %>%
      dplyr::select("TStart", "TStop",
                    "recur",
                    Tx.nm,
                    any_of(covariate_names),
                    all_of(event_indicator_string))
    testing_dataset_surv = test %>%
      filter(Status == 0) %>%
      dplyr::select("TStart", "TStop", "recur",
                    Tx.nm,
                    any_of(covariate_names),
                    all_of(event_indicator_string))

    # # weights
    # weight_name <- paste0("weight", suffix1)
    # weight = weights_rda(data = testing_dataset,
    #                      weight.formula.list = form.weight,
    #                      covariates = covariate_names,
    #                      event_indicator_string = event_indicator_string)
    # arg.val_name <- paste0("arg.val.", suffix1)
    # arg.val = list(test = testing_dataset,
    #                actual = test.tmp,
    #                propensity = weight$propensity,
    #                weight.censor = weight$weight.censor,
    #                criterion = criterion,
    #                tau = tau)

    weight_prop = cbind(id = test_ids,
                    weight = weight0$propensity[test_indices]) %>%
      as.data.frame()
    weight_censor = cbind(id = test_ids,
                         weight = weight0$weight.censor[test_indices]) %>%
      as.data.frame()
    # Expand weight to match the number of records for each person in the test dataset
    weight_prop_long <- test %>%
      # Join the weight dataset to the test dataset by 'id'
      left_join(weight_prop, by = "id") %>%
      pull(weight)
    weight_censor_long <- test %>%
      # Join the weight dataset to the test dataset by 'id'
      left_join(weight_censor, by = "id") %>%
      pull(weight)

    arg.val = list(test_dat = testing_dataset_surv,
                   actual = test.tmp,
                   propensity = weight_prop %>%
                     pull(weight), #weight$propensity[test_indices],
                   weight.censor = weight_censor %>%
                     pull(weight), #weight$weight.censor[test_indices],
                   criterion = criterion,
                   tau = tau,
                   endpoint1 = endpoint1)
    for (method in loop_methods) {
      err_method <- paste0("err.", method)
      opt_method <- paste0("opt.rule.", method)
      if (!get(err_method)) {
        # est1 = cbind(id = test_ids,
        #              est = get(opt_method)) %>% as.data.frame()
        # est2 = test %>%
        #   dplyr::select(id) %>%
        #   left_join(est1, by = "id") %>%
        #   pull(Trt) %>%
        #   as.data.frame()
        est2 = get(opt_method)

        # Truncated Mean/Prob
        values[cv, paste0(method, suffix1)] <-
          do.call(getValue,
                  c(arg.val,
                    list(estimated = est2)))
                           #get(opt_method))))
      }
    }
    # if (!err.CSK)  values[cv, "CSK"] = do.call(getValue, c(arg.val, list(estimated = opt.tmp.CSK)))
    # if (!err.zom)  values[cv, "ZOM"] = do.call(getValue, c(arg.val, list(estimated = opt.tmp.zom)))

    # observed
    obs1 <- paste0("OBS", suffix1)
    values[cv, obs1] <- do.call(getValue, c(arg.val, list(estimated = test.tmp)))
  } # end of suffix

  # source("Training_Eval_Results.R")
  message("Printing values")
  print(values[cv, ])
  message("Printing mean across cv iterations")
  print(apply(values, 2, mean, na.rm = TRUE)) # across columns (over all cv iterations)
  # View(values)
  ############################################################################################################
  #### C. saving the results #################################################################################
  ############################################################################################################
  attr(values, "spec") = data.frame(criterion1 = criterion_phase1[1],
                                    criterion2 = criterion_phase2[1],
                                    criterion1.s = as.numeric(criterion_phase1[2]),
                                    criterion2.s = as.numeric(criterion_phase2[2]),
                                    rule1 = rule1,
                                    rule2 = rule2,
                                    tau = tau,
                                    ert = ert,
                                    rs = rs,
                                    ntree = Ntree)
  saveRDS(values, paste0(fnm, "values_cv",cv,"_values_", nm, rds))
} # end of cv
saveRDS(values, paste0(fnm,"Values_", K,"CV_", nm, rds))
saveRDS(tab, paste0(fnm, "tab_ZOM_", nm, rds))

message("End of Bladder_analysis.R")
