library(caret); library(purrr); library(dplyr); library(survival)
library(DTRreg); library(randomForestSRC); library(dtrSurv)
library(itrSurv);library(ggplot2); library(tidyverse);
library(MASS); library(dplyr);
source("0_simulations/F03.ComparatorMethod_Functions.R")
source("F03.RDA_Functions.R")    # All the library and source files

K = 1 # number of CV
criterion_phase1 = "mean"
criterion_phase2 = criterion_phase1
endpoint = "CR"
dataset_name = "follic"
data(follic)
# make sure Treatment is "A" and time is "obs_time"
dat0 = follic %>%
  mutate(stage2 = as.numeric(clinstg == 2),
         Trt = as.numeric(ch == "Y"),
         D.0 = ifelse(status > 0, 1, 0),
         D.1 = ifelse(status == 1, 1, 0),
         D.2 = ifelse(status == 2, 1, 0)) %>%
  dplyr::select(A = Trt, stage2, age, hgb, obs_time = time, status, D.0, D.1, D.2)
range(dat0$obs_time)

### 0. parameters
Tx.nm = "A"
covariates = c("stage2", "age", "hgb")
## time points
tau = 32
timepoints = seq(0, sqrt(tau), length.out = 1000)^2

## other parameters
nodesize = 50
mindeath = round(sqrt(c(nodesize)), 0)
Ntree = 300
ert = TRUE; rs = 0.2 # randomSplit = 0.2

if (criterion_phase1 == "mean.prob.combo"){
  dtr_criterion = "surv.mean"
} else if (criterion_phase1 == "prob"){
  dtr_criterion = "surv.prob"
} else{
  dtr_criterion = criterion_phase1
}

### 0.2 criterion
# for (criterion_phase1 in c("mean", "mean.prob.combo")){ # for (dtr_criterion in c("mean", "surv.mean")) {
if (criterion_phase1[1] != "mean") {
  criterion_phase1[2] = 2200 # six-year survival
  dtr_criterion[2] = criterion_phase1[2]
  rule = "logrank"
} else {
  rule = "mean"
  criterion_phase1[2] = NA
}

if (is.null(criterion_phase1[2]) | is.na(criterion_phase1[2])){
  t0_pmcr = tau/2
  crit_tmp = ""
} else{
  t0_pmcr = criterion_phase1[2]
  crit_tmp = round(as.numeric(criterion_phase1[2]), 1)
}

cat(#"all_cause = ", all_cause,
  "\ncriterion_phase1 =",criterion_phase1[1],
  "\nvalue.s =", criterion_phase1[2], "\ntau =", tau, "\n")

## Overall Survival PM models
modelOS = create_model_formulas("Surv(obs_time, D.0)",
                                covariates,
                                formula = FALSE)
form.OS <- modelOS %>% as.formula()

## propensity models
modelPr_PMCR = create_model_formulas("A",
                                     covariates,
                                     formula = TRUE)
modelPr = create_model_formulas("factor(A)",
                                covariates,
                                formula = FALSE)
form.weight <- list(modelPr %>% as.formula)

### 1. data preprocessing
dat = dat0 %>%
  mutate(#time = pmin(time, tau),
    D.0 = ifelse(status > 0, 1, 0),
    D.1 = ifelse(status == 1, 1, 0),
    D.2 = ifelse(status == 2, 1, 0)
  )
#        last.fu = ifelse(time == tau, 1, last.fu),
#        d.1 = ifelse(time == tau & last.fu == 1, 1, d.1),
#        d.2 = ifelse(time == tau, ifelse(last.fu == 1, NA, 1), d.2),
#        delta = ifelse(time == tau, 1, delta))

mean(dat$status==0, na.rm = TRUE)

lvls = lapply(Tx.nm, function(x) levels(dat[, x]))
for (x in 1:length(lvls)) {if (is.null(lvls[[x]])) lvls[[x]] = 0:1}

## cross-validation
K = K
set.seed(100, kind="Mersenne-Twister")
cv.insample  = createDataPartition(dat[, Tx.nm],
                                   p = .8,
                                   list = FALSE,
                                   times = K)
cv.outsample = map_dfc(1:K, ~which(!(1:dim(dat)[1] %in% cv.insample[, .]))) %>% as.matrix

# skeleton
values_colsnames = c("CZMK", "CSK", "PMCR",
                     "ZOM", "CSKZOM", "observed",
                     "ns.CZMK", "ns.CSK",
                     "cens1", "cause1", "cause2")#
values <- matrix(NA, K, length(values_colsnames),
                 dimnames = list(1:K, values_colsnames))

nm = paste0(criterion_phase1[1],
            crit_tmp, "_rule",
            rule, "_tau", tau, all_cause_nm)
fnm = paste0("../3_output/", endpoint, "/",dataset_name,"/", Sys.Date(), "/") # folder name
if (!dir.exists(sprintf("../3_output/%s", endpoint))) dir.create(sprintf("../3_output/%s", endpoint))
if (!dir.exists(sprintf("../3_output/%s/%s", endpoint, dataset_name))) dir.create(sprintf("../3_output/%s/%s", endpoint, dataset_name))
if (!dir.exists(fnm)) dir.create(fnm)
rds = paste0("_", Sys.Date(), ".rds")

tab = data.frame(s1 = rep(NA, K),
                 # s2 = NA,
                 total = NA)

for (cv in 1:K) {
  cat(cv, "th cv.\n")
  set.seed(cv)
  in.cv = cv.insample[, cv]   # insample index
  print(tail(in.cv))
  out.cv = cv.outsample[, cv] # outsample index
  train = dat[in.cv, ]
  test  = dat[out.cv,]
  values[cv, "cens1"] = mean(train$status==0); # proportion of people censored
  values[cv, "cause1"] = mean(train$status==1, na.rm = T)
  values[cv, "cause2"] = mean(train$status==2, na.rm = T)

  ############################################################################################################
  #### A. rule estimation ####################################################################################
  ############################################################################################################

  ### A1. The proposed method
  priority_cause = 1
  priority_vector <- c(0, priority_cause)
  # Use lapply to create a list of sprintf statements
  models_itr <- lapply(priority_vector, function(x) {
    paste0(sprintf("Surv(obs_time, D.%s) ~ ",
                   ifelse(x == 0, "0", as.character(x))),
           gsub(".*~\\s*", "", modelOS)) %>% as.formula
  })
  args.CZMK <- list(data = train,
                    txName = Tx.nm,
                    models = models_itr,
                    tau = tau, timePoints = timepoints,
                    criticalValue1 = criterion_phase1[1],
                    criticalValue2 = criterion_phase2[1],
                    evalTime = as.numeric(criterion_phase1[2]),
                    splitRule = ifelse(criterion_phase1[1] == "mean", "mean", "logrank"),
                    ERT = ert, uniformSplit = ert, replace = !ert,
                    randomSplit = rs, nTree = Ntree, mTry = 6,
                    pooled = FALSE,
                    tol1 = c(0.1,0),
                    stratifiedSplit = FALSE)# actual fitting
  values[cv, "ns.CZMK"] = nodesize
  set.seed(cv)
  CZMK.i <-
    try(do.call(itrSurv,
                c(args.CZMK, list(nodeSize = nodesize,
                                  minEvent = mindeath))))
  err.CZMK = class(CZMK.i)[1] == "try-error"

  ### A2. Cho et al (2022)
  args.CSK <- list(data = train,
                   txName = Tx.nm,
                   models = form.OS,
                   usePrevTime = FALSE, tau = tau, timePoints = timepoints,
                   criticalValue = dtr_criterion[1],
                   evalTime = as.numeric(dtr_criterion[2]),
                   splitRule = ifelse(dtr_criterion[1] == "mean", "mean", "logrank"),
                   ERT = ert, uniformSplit = ert, replace = !ert,
                   randomSplit = rs, nTree = Ntree, mTry = 6,
                   pooled = FALSE, stratifiedSplit = FALSE)

  # actual fitting
  values[cv, "ns.CSK"] = nodesize
  set.seed(cv)
  CSK.i <-
    try(do.call(dtrSurv,
                c(args.CSK, list(nodeSize = nodesize,
                                 minEvent = mindeath))))
  err.CSK = class(CSK.i)[1] == "try-error"

  ### A3. PMCR - 2021
  message("unrestricted regime")
  set.seed(cv)
  temp.unrestr.fit<-PMCR(Time="obs_time",
                         Event="status",
                         formula=modelPr_PMCR, #trt needs to be 0/1
                         data=train,
                         rgenoud=FALSE,
                         Restrict=FALSE,
                         propscore="logistic",
                         t0=t0_pmcr)
  #range of alpha # alpha needs to be tuned
  alps<-c(temp.unrestr.fit$Fbeta2[2],
          seq(round(temp.unrestr.fit$Fbeta2[2],2)+0.01,
              round(temp.unrestr.fit$Fbeta1[2],2)+0.03,
              0.01))
  alp<-alps[2] # alp needs to be tuned based on alps
  message("restricted regime")
  M_pmcr<-1e+5
  temp.restr.fit<-PMCR(Time="obs_time",
                       Event="status", # CR event indicator
                       formula=modelPr_PMCR, # Trt must be 0/1
                       data=train,
                       rgenoud=FALSE,
                       Restrict=TRUE,
                       propscore="logistic",
                       t0=t0_pmcr,
                       alp=alp,
                       M=M_pmcr)
  # dat1 = temp.restr.fit$data
  # cov = dat1 %>%
  #   dplyr::select(-c(id, obs_time, status, pix.hat)) %>%
  #   mutate(int = rep(1, nrow(dat1))) %>%
  #   dplyr::select(int, covariates)
  PMCR.i <- temp.restr.fit$beta3
  err.PMCR = class(PMCR.i)[1] == "try-error"

  # ### A4. Goldberg-Kosorok RF
  # args.GK = list(common.formula = form.GK,
  #                common.Tx.label = Tx.nm, stage.label = 1:2, tau = tau,
  #                data = train, stage.sep = ".", regress.prev.time = T)
  #
  # # actual fitting
  # values[cv, "ns.GK"] = nodesize
  # set.seed(cv)
  # GKRF.i <-
  #   try(do.call(gk.separate, c(args.GK, list(nodesize = nodesize, method = "rf"))))
  # err.GKRF = class(GKRF.i)[1] == "try-error"


  ### A6. zero-order model
  zom.i <-
    try(do.call(itrSurv, c(args.CZMK,
                           list(nodeSize = 1e+4,
                                minEvent = 1e+4))))
  err.zom = class(zom.i)[1] == "try-error"

  if (!err.zom) {
    zom.pred =
      data.frame(
        zom.i@phaseResults$FinalOptimalTx_Recc[1])
    tab[cv, 1] = zom.i@phaseResults$FinalOptimalTx_Recc[1]
    cat("ZOM freq table\n"); tab[, 1] %>% table %>% print
  }

  ### A7. CSKzero-order model
  CSKzom.i <-
    try(do.call(dtrSurv, c(args.CSK,
                           list(nodeSize = 1e+4,
                                minEvent = 1e+4))))
  err.CSKzom = class(CSKzom.i)[1] == "try-error"

  if (!err.CSKzom) {
    CSKzom.pred =
      data.frame(
        CSKzom.i@stageResults[[1]]@optimal@optimalTx[1])
    tab[cv, 1] = CSKzom.i@stageResults[[1]]@optimal@optimalTx[1]
    cat("CSKZOM freq table\n");
    tab[, 1] %>% table %>% print
  }

  ############################################################################################################
  #### B. Rules applied to a test set ########################################################################
  ############################################################################################################
  ### B0. skeletons for the predicted A's
  opt.rule.pred = data.frame(A = rep(NA, dim(cv.outsample)[1]))
  opt.rule.pred[,1] = factor(NA, levels = lvls[[1]])
  opt.rule.CZMK <- opt.rule.CSK <- opt.rule.PMCR <- opt.rule.ZOM <- opt.rule.CSKzom <- opt.rule.pred
  # opt.rule.GKRF <- opt.rule.GKLM <-
  elig = !test[, Tx.nm] %>% is.na

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
  }
  ### B2. Cho et al (2022)
  if (!err.CSK) {
    opt.CSK =
      predict(CSK.i,
              # newdata = test[elig, ],
              newdata = test[elig, ],
              stage = 1)
    opt.rule.CSK[elig, ] = factor(opt.CSK$optimal@optimalTx,
                                  levels = lvls[[1]]) # %>% as.numeric()
  }

  ### B3. PMCR - 2021
  if (!err.PMCR) {
    cov2 = test[elig, ] %>%
      dplyr::select(-c(obs_time, status)) %>%
      mutate(int = rep(1, nrow(test[elig,]))) %>%
      dplyr::select(int, covariates)

    opt.PMCR = ifelse(
      (as.matrix(cov2) %*% as.matrix(PMCR.i)) > 0,
      1, 0) %>% as.vector()

    opt.rule.PMCR[elig, ] = factor(opt.PMCR,
                                   levels = lvls[[1]]) # %>% as.numeric()
  }

  ### B4. Goldberg-Kosorok RF
  # if (!err.GKRF) {
  #   opt.GKRF = predict.opt.sep(GKRF.i$survRF[[1]],
  #                              newdata = test[elig, ] %>% dplyr::select(-V.2, -d.2) %>%
  #                                mutate(prevTime.2 = V.1),
  #                              Tx.label = paste0(Tx.nm, ".", q))
  #   opt.rule.GKRF[elig, q] = lvls[[1]][opt.GKRF$optimal.Tx]
  # }
  # ### B5. Goldberg-Kosorok linear
  # if (!err.GKLM) {
  #   opt.GKLM = predict.opt.sep(GKLM.i$survRF[[1]],
  #                              newdata = test[elig, ] %>% dplyr::select(-V.2, -d.2) %>%
  #                                mutate(prevTime.2 = V.1),
  #                              Tx.label = paste0(Tx.nm, ".", q))
  #   opt.rule.GKLM[elig, q] = lvls[[1]][opt.GKLM$optimal.Tx]
  # }

  ### B6. zero-order model
  if (!err.zom) {
    opt.rule.ZOM[elig, ] = zom.pred[1, ]
    opt.rule.ZOM[elig, ] = factor(zom.pred[1, ], levels = lvls[[1]]) # %>% as.numeric()
  }

  ### B7. CSKzero-order model
  if (!err.CSKzom) {
    opt.rule.CSKzom[elig, ] = CSKzom.pred[1, 1]
    opt.rule.CSKzom[elig, ] = factor(CSKzom.pred[1, 1],
                                     levels = lvls[[1]]) # %>% as.numeric()
  }
  print(table(opt.rule.CZMK, useNA = "always"))
  print(table(opt.rule.CSK, useNA = "always"))
  print(table(opt.rule.PMCR, useNA = "always"))
  print(table(opt.rule.ZOM, useNA = "always"))
  print(table(opt.rule.CSKzom, useNA = "always"))

  ############################################################################################################
  #### C. Value estimation ###################################################################################
  ############################################################################################################
  # weights
  weight = weights_rda(data = test,
                             weight.formula.list = form.weight,
                             covariates = covariates)
  test.tmp = test %>% transmute(A = test[, Tx.nm])

  arg.val = list(test = test,
                 actual = test.tmp,
                 propensity = weight$propensity,
                 weight.censor = weight$weight.censor,
                 criterion = criterion_phase1,
                 tau = tau)
  if (!err.CZMK)  values[cv, "CZMK"] = do.call(getValue,
                                               c(arg.val,
                                                 list(estimated = opt.rule.CZMK)))
  if (!err.CSK)  values[cv, "CSK"] = do.call(getValue,
                                             c(arg.val,
                                               list(estimated = opt.rule.CSK)))
  if (!err.PMCR)  values[cv, "PMCR"] = do.call(getValue,
                                               c(arg.val,
                                                 list(estimated = opt.rule.PMCR)))
  if (!err.zom)  values[cv, "ZOM"] = do.call(getValue,
                                             c(arg.val,
                                               list(estimated = opt.rule.ZOM)))
  if (!err.CSKzom)  values[cv, "CSKZOM"] = do.call(getValue,
                                                   c(arg.val,
                                                     list(estimated = opt.rule.CSKzom)))
  values[cv, "observed"] = do.call(getValue,
                                   c(arg.val[-which(names(arg.val) == "propensity")],
                                     list(estimated = test.tmp,
                                          propensity = 1)))

  print(values[cv, ])
  print(apply(values, 2, mean, na.rm = TRUE)) # across columns (over all cv iterations)

  ############################################################################################################
  #### C. saving the results #################################################################################
  ############################################################################################################

  attr(values, "spec") = data.frame(criterion = criterion_phase1[1],
                                    criterion.s = as.numeric(criterion_phase1[2]),
                                    rule = rule,
                                    tau = tau,
                                    ert = ert,
                                    rs = rs,
                                    ntree = Ntree,
                                    nodesize = nodesize,
                                    mindeath = mindeath)
  if (cv==1) saveRDS(CZMK.i, paste0(fnm, "czmk.cv",cv, "_", dataset_name, "_", nm, "_", cv, rds))
  if (cv==1) saveRDS(CSK.i, paste0(fnm, "csk.cv", cv, "_", dataset_name, "_", nm, "_", cv, rds))
  if (cv==1) saveRDS(PMCR.i, paste0(fnm, "pmcr.cv", cv, "_", dataset_name,"_", nm, "_", cv, rds))
  saveRDS(values, paste0(fnm, "values_cv",cv,"_values_", nm, rds))
} # end of cv
saveRDS(values, paste0(fnm,"OSValues_", K,"CV_", nm, rds))
saveRDS(tab, paste0(fnm, "tab_ZOM_", nm, rds))
# } # end of criterion.value_phase1
