library(caret);
library(purrr);
library(dplyr);
library(survival)
library(randomForestSRC);
library(itrSurv);library(ggplot2); library(tidyverse);
library(MASS); library(dplyr);
local = 1
setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper3_RE/")
source("02.Simulation_Functions_RE.R") # includes F02.ComparatorMethod_Functions.R
setwd("../../RDA/Paper3_RE/")
source("RDA_Functions_RE.R")    # All the library and source files
############################################################
######### temporary testing things out parameters ##########
############################################################
criterion_phase1 = "mean"
rule2 = "gray_re"
############################################################
################# methods for comparison ###################
############################################################
rda_methods = c("CZMK", "ZOM", "observed")
skip_method <- !c(TRUE,TRUE,TRUE);
assign_skip_function(rda_methods, skip_method)
skipped_methods <- rda_methods[skip_method]
loop_methods <- rda_methods[!rda_methods %in% c(skipped_methods, "observed")]
############################################################
#################### 0. parameters ####################
############################################################
tol1_param = c(0.1,0,0.3,0.01)
t0_crit = 30 #2200 # 6-mo survival
pooled1 = FALSE # stratified = lower nodesize; pooled = can have larger nodesize
tau = 59
K = 1 #300 # number of CV
endpoint = "RE" # endpoint
Tx.nm = "A"

## other parameters
nodeSizeSurv = 5
nodeSizeEnd = 5
minEventSurv = round(sqrt(c(nodeSizeSurv)), 0)
minEventEnd = round(sqrt(c(nodeSizeEnd)), 0)
Ntree = 1
ert = FALSE
rs = 0.2 # randomSplit = 0.2

############################################################
##################### Bladder Dataset ######################
############################################################
library(tidyverse)
# toy dataset with RE and terminal death (fake added in)

# bladder2 dataset generated from SAS code #
file_path <- "~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Scratch/Project3/bladder2.csv"
# Read the CSV file into a data frame, including headers
bladder2 <- read.csv(file_path, header = TRUE)
colnames(bladder2) = c("id", "TStart", "TStop", "A", "Z1", "Z2", "Visit", "Status", "Gaptime")

# this dataset doesnt have terminal events so just
# BS-ing a few in for now to have that variable as we
# test the package/code up the method
# Step 1: Create a new variable status_D that is initially a copy of status
bladder2$Status_D <- bladder2$Status
# Step 2: For rows where status is 0, randomly set 90% to 1 in status_D
# Create a logical vector to identify rows where status is 0
zero_status_rows <- bladder2$Status == 0
# Randomly sample 90% of these rows and set them to 1 in status_D
set.seed(123) # Set seed for reproducibility, if needed
num_to_change <- sum(zero_status_rows) * 0.9
change_indices <- sample(which(zero_status_rows), num_to_change)
bladder2$Status_D[change_indices] <- 1
# Step 3: Set Status_D to 0 where Status is 1
# Status = 1 for recurrent event and = 0 for censored or death
# Status_D = 1 for death and = 0 otherwise
bladder2$Status_D[bladder2$Status == 1] <- 0
# Check the result
head(bladder2)

# make it so people who die at time TStart = TStop have an interval where TStart = TStop + 0.001
bladder2.1 <- bladder2 %>%
  mutate(
    TStop = ifelse(TStart == TStop, TStop + 0.001, TStop)  # Avoids TStart = TStop
  )


# need to add in censoring rows for the last dataset if status_d = 0
df_tmp = bladder2.1 %>% arrange(id, Visit) %>%
  group_by(id) %>%
  mutate(
    need_to_add_row = ifelse(!(row_number() == n() & Status_D == 0 & Visit > 1), FALSE,
                             ifelse(Status == 0 & Status_D == 0, FALSE,
                                    TRUE))) %>%
  ungroup() %>%
  mutate(
    max_stop_time = max(TStop)
  )
new_rows <- df_tmp %>%
  filter(need_to_add_row) %>%
  mutate(
    TStart1 = TStop,
    TStop1= max_stop_time,
    Visit1 = Visit + 1,
    Gaptime1 = max_stop_time - TStop,
    Status1 = 0
  ) %>%
  dplyr::select(id, TStart = TStart1, TStop = TStop1, A, Z1, Z2, Visit=Visit1,
                Status = Status1, Gaptime=Gaptime1, Status_D)
new_rows
# Combine the original data with the new rows
df <- bind_rows(bladder2.1, new_rows) %>%
  group_by(id) %>%
  arrange(id, Visit) %>%
  ungroup() %>%
  as.data.frame()
df
bladder_surv = df %>%
  filter(Status == 0)
#checking no. subjects
# bladder2.1 %>% as.data.frame() %>% dplyr::select(id) %>% distinct() %>% nrow()
# bladder_surv %>% as.data.frame() %>% dplyr::select(id) %>% distinct() %>% nrow()
# df %>% as.data.frame() %>% dplyr::select(id) %>% distinct() %>% nrow()

dat1 = bladder_surv # Phase 1 Surv
dat2 = df # Phase 2 Endpoint (RE)

bladder_df = dat2
dataset_name = "bladder"

bladder_df$id %>% unique() %>% length() # 12 patients (6-12, 56-60)
sub_ids = bladder_df$id %>% unique()

data_to_use = bladder_df %>%
  dplyr::select(-c(Visit, Gaptime))
a1 = data_to_use %>% filter(A == 1) ;
end_a1 = a1 %>%
  dplyr::select(id, Z1, Z2, Status, Status_D, TStop)
surv_a1 = a1 %>%
  filter(Status == 0) %>%
  dplyr::select(id, Z1, Z2, Status, Status_D)
recurr_a1 = a1 %>%
  filter(Status == 1) %>%
  dplyr::select(id, Z1, Z2, Status, Status_D)
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
                           "Status_D", "Status")
# Exclude the specified column names
covariate_names <- setdiff(all_column_names, excluded_column_names)

# manipulating dataset to obtain truncated status,D.0,D.1,D.2
dat0 = data_to_use %>%
  dplyr::select(id, TStart, TStop, A,
                Status_D, Status,
                !!!covariate_names)

#################### 1. data preprocessing
dat = as.data.frame(dat0)

mean(dat$Status_D==1, na.rm = TRUE)
mean(dat$Status==1, na.rm = TRUE)

lvls = lapply(Tx.nm, function(x) levels(dat[, x]))
for (x in 1:length(lvls)) {if (is.null(lvls[[x]])) lvls[[x]] = 0:1}


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
# Exclude "observed" for train_values and test_values
traintest_methods <- setdiff(rda_methods, "observed")
train_values = c(".train.OS", ".train.PropPhase2", ".train.RE")
test_values = c(".test.OS", ".test.PropPhase2", ".test.RE")
values_colsnames <- c(
  unlist(lapply(traintest_methods, function(method) paste0(method, train_values))),
  unlist(lapply(traintest_methods, function(method) paste0(method, test_values))),
  # paste0(rda_methods, ".OS"),
  # paste0(rda_methods, ".RE"),
  "ns1.CZMK",
  "ns2.CZMK",
  "train_cens", "train_RE", "train_terminal"
)
values <- matrix(NA, K, length(values_colsnames),
                 dimnames = list(1:K, values_colsnames))
train_eval_result_list = train_eval_result_list_tmp = list()

nm = paste0(criterion_phase1[1],
            crit_tmp, "_rule1",
            rule1, "_rule2", rule2, "_tau", tau)
fnm = paste0("./output/", endpoint, "/",dataset_name,"/", Sys.Date(), "/") # folder name
if (!dir.exists(sprintf("./output/%s", endpoint))) dir.create(sprintf("./output/%s", endpoint))
if (!dir.exists(sprintf("./output/%s/%s", endpoint, dataset_name))) dir.create(sprintf("./output/%s/%s", endpoint, dataset_name))
if (!dir.exists(sprintf("./figure/%s", endpoint))) dir.create(sprintf("./figure/%s", endpoint))
if (!dir.exists(sprintf("./figure/%s/%s", endpoint, dataset_name))) dir.create(sprintf("./figure/%s/%s", endpoint, dataset_name))
if (!dir.exists(fnm)) dir.create(fnm)
rds = paste0("_", Sys.Date(), ".rds")

tab = data.frame(s1 = rep(NA, K),
                 # s2 = NA,
                 total = NA)

## cross-validation
K = K
# Create cross-validation folds for individuals
if (K > 1){
  folds <- createFolds(sub_ids, k = K, list = TRUE)
} else{
  folds <- createFolds(sub_ids, k = 2, list = TRUE)
}

for (cv in 1:K) {
  cat(cv, "th cv.\n")
  set.seed(cv)

  # Get in-sample (training) and out-sample (testing) IDs
  test_ids <- sub_ids[folds[[cv]]]
  train_ids <- setdiff(sub_ids, test_ids)

  # Subset data based on IDs
  train <- dat %>% filter(id %in% train_ids) %>% as.data.frame()
  # View(train)
  test  <- dat %>% filter(id %in% test_ids) %>% as.data.frame()

  # Debug print statements
  print(paste("Train size:", nrow(train), "with Unique IDs in Train:", length(unique(train$id))))
  print(paste("Test size:", nrow(test), "with Unique IDs in Test:", length(unique(test$id))))

  values[cv, "train_cens"] = mean(train$Status_D==0 & train$Status == 0); # proportion of people censored
  values[cv, "train_RE"] = mean(train$Status==1, na.rm = T)
  values[cv, "train_terminal"] = mean(train$Status_D==1, na.rm = T)

  ############################################################################################################
  #### A. rule estimation ####################################################################################
  ############################################################################################################

  ### A1. The proposed method
  if (skip.CZMK != TRUE){
    print(sprintf("Running CZMK for %s CV", cv))

    model1 = "Surv(TStop, Status_D) ~ Z1 + Z2" %>% as.formula()
    model2 = "Surv(TStart, TStop, Status) ~ Z1 + Z2" %>% as.formula()
    models_RE = list(model1, model2)

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
                      tol1 = c(0.1,0.1),
                      stratifiedSplit = 0.1)
    values[cv, "ns1.CZMK"] = nodeSizeSurv
    values[cv, "ns2.CZMK"] = nodeSizeEnd
    set.seed(cv)
    CZMK.i <-
      try(do.call(itrSurv::itrSurv,
                  c(args.CZMK, list(nodeSizeSurv = nodeSizeSurv,
                                    nodeSizeEnd = nodeSizeEnd,
                                    minEventSurv = minEventSurv,
                                    minEventEnd = minEventEnd))))
    values[cv, "CZMK.train.OS"] = CZMK.i@value[["V1"]][["Et_survival"]]
    values[cv, "CZMK.train.PropPhase2"] = CZMK.i@value[["V2"]][["PropPhase2"]]
    values[cv, "CZMK.train.RE"] = CZMK.i@value[["V3"]][["Et_mff"]]
    err.CZMK = class(CZMK.i)[1] == "try-error"
  } else{
    err.CZMK = TRUE
  }

  ### A6. zero-order model
  if (skip.ZOM != TRUE){
    print(sprintf("Running ZOM for %s CV", cv))
    ZOM.i <-
      try(do.call(itrSurv::itrSurv, c(args.CZMK,
                                      list(nodeSizeEnd = 1e+4,
                                           nodeSizeSurv = 1e+4,
                                           minEventEnd = 1e+4,
                                           minEventSurv = 1e+4))))
    values[cv, "ZOM.train.OS"] = ZOM.i@value[["V1"]][["Et_survival"]]
    values[cv, "ZOM.train.PropPhase2"] = ZOM.i@value[["V2"]][["PropPhase2"]]
    values[cv, "ZOM.train.RE"] = ZOM.i@value[["V3"]][["Et_mff"]]
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
  opt.rule.pred = data.frame(Trt = rep(NA, dim(test)[1]))
  opt.rule.pred[,1] = factor(NA, levels = lvls[[1]])

  for(method in loop_methods) {
    assign(paste("opt.rule.", method, sep = ""),
           opt.rule.pred,
           envir = .GlobalEnv)
  }
  # Create eligibility for each row based on non-missing values in `Tx.nm`
  elig <- test %>%
    dplyr::select(id, !!sym(Tx.nm)) %>%
    group_by(id) %>%
    mutate(eligibility = !any(is.na(.data[[Tx.nm]]))) %>%
    ungroup() %>%
    pull(eligibility)

  ### B1. the proposed method
  if (!err.CZMK) {
    for (phase in 1:2){
      opt.CZMK_phase =
        predict(CZMK.i,
                newdata = test[elig, ],
                Phase = phase)
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
    opt.rule.CZMK[elig, ] = factor(opt.CZMK,
                                   levels = lvls[[1]]) # %>% as.numeric()
    print(table(opt.rule.CZMK, useNA = "always"))
  }

  ### B6. zero-order model
  if (!err.ZOM) {
    opt.rule.ZOM[elig, ] = zom.pred[1, ]
    opt.rule.ZOM[elig, ] = factor(zom.pred[1, ], levels = lvls[[1]]) # %>% as.numeric()
    print(table(opt.rule.ZOM, useNA = "always"))
  }

  ############################################################################################################
  #### C. Value estimation ###################################################################################
  ############################################################################################################
  test.tmp = test %>% transmute(Trt = test[, Tx.nm])

  if (criterion_phase1[1] != criterion_phase2[1]){
    message("WARNING: CRITERION_PHASE1 AND CRITERION_PHASE2 ARE DIFFERENT - MAKE SURE THAT IS WHAT YOU WANT")
  }
  suffixes = c("OS", "RE")
  for (suffix in 1:length(suffixes)){
    suffix1 = suffixes[suffix]
    if (suffix1 == "OS"){
      event_indicator_string = "Status_D"
      criterion = criterion_phase1
    } else if (suffix1 == "RE"){
      event_indicator_string = "Status"
      criterion = criterion_phase2
    }
    testing_dataset = test %>%
      dplyr::select("TStart", "TStop",
                    Tx.nm,
                    any_of(covariate_names),
                    all_of(event_indicator_string))

    # weights
    weight_name <- paste0("weight.", suffix1)
    weight = weights_rda(data = testing_dataset,
                         weight.formula.list = form.weight,
                         covariates = covariate_names,
                         event_indicator_string = event_indicator_string)
    arg.val_name <- paste0("arg.val.", suffix1)
    arg.val = list(test = testing_dataset,
                   actual = test.tmp,
                   propensity = weight$propensity,
                   weight.censor = weight$weight.censor,
                   criterion = criterion,
                   tau = tau)

    for (method in loop_methods) {
      err_method <- paste0("err.", method)
      opt_method <- paste0("opt.rule.", method)
      if (!get(err_method)) {
        # Truncated Mean/Prob (anything other than OS/PC is counted as censored)
        values[cv, paste0(method, ".", suffix1)] <-
          do.call(getValue,
                  c(arg.val,
                    list(estimated =
                           get(opt_method))))
      }
    }

    # observed
    obs1 <- paste0("observed.", suffix1)
    values[cv, obs1] <- do.call(getValue,
                                c(arg.val[-which(names(arg.val) == "propensity")],
                                  list(estimated = test.tmp,
                                       propensity = 1)))
  } # end of suffix

  # source("Training_Eval_Results.R")

  print(values[cv, ])
  print(apply(values, 2, mean, na.rm = TRUE)) # across columns (over all cv iterations)

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
                                    ntree = Ntree,
                                    nodesize = nodesize,
                                    mindeath = mindeath)
  saveRDS(values, paste0(fnm, "values_cv",cv,"_values_", nm, rds))
} # end of cv
saveRDS(values, paste0(fnm,"Values_", K,"CV_", nm, rds))
saveRDS(tab, paste0(fnm, "tab_ZOM_", nm, rds))
