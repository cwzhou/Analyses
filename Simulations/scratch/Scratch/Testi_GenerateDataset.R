source("~/Desktop/UNC BIOS PhD/Dissertation_KosorokLab/Github/Paper_1/Paper1_Sims/Sims_SingleStage/F01.dynamics_SS.R")

#' @param p number of covariates
#' @param Sig covariate error covariance
#' 
#' @return event.time observed time
#' @return gamma 1(death <= treatment)
#' @return delta 1(min(death, treatment) <= censor)
#' @at.risk Availability at the beginning of the stage
#' @covariates(Z1-Zp) covariate values at the beginning of the stage
#' 

# n.sample = arg.csk$n.sample
# predHazardFn = arg.csk$predHazardFn; predCensorFn = arg.csk$predCensorFn;
# predPropensityFn = arg.csk$predPropensityFn
# hidden_data = arg.csk$hidden_data
# printFlag = arg.csk$printFlag
# policy = arg.csk$policy
# at.risk = arg.csk$at.risk
# corr = arg.csk$corr
# rho = omega = covariate = NULL
# surv.previous = rep(0, n.sample)


multiStageDynamics <- 
  function(n.sample = 100, n.stages, tau = 10, tick = 0.01,         # structural parameters
           covariate = NULL,           
           p = 5, ncauses, cause1_prob = NULL,
           Sig = diag(p) + 0.2 - diag(p) * 0.2,              # covariate structure
           predHazardFn, predPropensityFn,                          # list of predictor functions
           hidden_data = FALSE, summary = TRUE,                     # output control
           policy = NULL,                                           # optimal rule (if !is.null, propensity scores are ignored.) for value calculation
           printFlag = TRUE                        # policy is an object of rsf.obj list.
  ) {
    # p = arg.obs$p;
    # ncauses = arg.obs$ncauses
    # caues1_prob = arg.obs$cause1_prob;
    # n.sample = arg.obs$n.sample; n.stages = arg.obs$n.stages; tau = arg.obs$tau;
    # tick = arg.obs$tick; hidden_data = arg.obs$hidden_data; printFlag = arg.obs$printFlag;
    # predHazardFn = arg.obs$predHazardFn; predPropensityFn = arg.obs$predPropensityFn;
    # covariate = NULL;
    # Sig = diag(p) + 0.2 - diag(p) * 0.2
    require(dplyr)
    if (!is.null(covariate)) if (dim(covariate)[1] != n.sample) stop ("length of at.risk should match.")
    time.max = max(tau - 0 )
    
    # initial state
    if (is.null(covariate))     covariate = mvrnorm(n.sample, rep(0, p), Sig)
    if (is.null(colnames(covariate))) colnames(covariate) = paste0("Z", 1:p)
    # skeleton
    tmp <-
      gdata_CR(time.max = 1, tick = tick, 
               cause1_prob = cause1_prob,
               censor_rate = 0.2,
               pred.hazard = 0, 
               full_data = ifelse(hidden_data, 2, 0)) %>% t
    
    output <- 
      matrix(NA, nrow = n.sample, ncol = 2 + dim(tmp)[2] + 2 + p)
    colnames(output) = c("subj.id", "rep.id", colnames(tmp), "action", "at.risk", colnames(covariate))
    output[, "subj.id"] <- 1:n.sample
    output[, "rep.id"] <- 0           # rep.id is reserved for later use (repeated random draws)
    action = rep(NA, n.sample)           # initialize action vector
    
    if (printFlag) cat("endpoint ", endpoint, "\n")
    output[, "at.risk"]           <- at.risk
    output[, colnames(covariate)] <- covariate
    ## action                (t=0)
    if (!is.null(policy)) {
      print("policy is not null")
      #         if (attr(policy, "class") %in% c("Goldberg-lm")) {
      #           vars <- policy[[stage]]$model %>% names
      #           vars <- vars[!grepl("(weights)", vars)]  # exclude (weights)
      #           args <- c(list(policy[[stage]],
      #                          newdata = output[as.logical(at.risk), c("subj.id", vars), stage] %>% as.data.frame,
      #                          return.optimal.Q = FALSE))
      #           action[at.risk != 0] <- do.call(decision.rule.gk, args)$rule.fix
      #         } else if (attr(policy, "class") %in% c("Goldberg-rf")) {
      #   args <- c(list(policy[[stage]],
      #                  newdata = output[as.logical(at.risk), c("subj.id", policy[[stage]]$xvar.names), stage] %>% as.data.frame,
      #                  return.optimal.Q = FALSE))
      #   action[at.risk != 0] <- do.call(decision.rule.gk, args)$rule.fix
      # } else 
      if (attr(policy, "class") %in% c("dwSurv")) {
        # tmp.out <<- output      
        dat.wide <- dwTrans(output, n.stages = n.stages, p = p)
        # vars.tf   <- names(get_all_vars(policy$tf.mod[[stage]], dat.wide))
        # vars.blip <- names(get_all_vars(policy$blip.mod[[stage]], dat.wide))
        # vars <- unique(c(vars.tf, vars.blip))
        args <- c(list(policy,
                       newdata = dat.wide[as.logical(at.risk), ],
                       # output[as.logical(at.risk), c("subj.id", vars), stage] %>% as.data.frame,
                       stage = stage,
                       return.optimal.Q = FALSE))
        action[at.risk != 0] <- do.call(decision.rule.dw, args)$rule.fix
      } else if (attr(policy, "class") %in% c("cho")) {  #### old code
        args <- c(list(policy[[stage]],
                       newdata = output[as.logical(at.risk), c("subj.id", policy[[stage]][[1]]$xvar.names), stage] %>% as.data.frame,
                       return.optimal.s = FALSE, tau = tau),
                  attr(policy[[stage]], "criterion"))
        action[at.risk != 0] <- do.call(decision.rule.list, args)$rule.fix
      } else if (attr(policy, "class") %in% c("DTRSurv")) { #### new package!
        x = output[as.logical(at.risk),,] %>% output2observable()
        x[, paste0(c("T_", "delta_"), stage)] = 0  # dummy values in order for the package not to drop the NA rows.
        x = get_all_vars(policy@stageResults[[stage]]@model, x)
        args <- list(policy, newdata = x, stage = stage)
        action[at.risk != 0] <- do.call(predict, args)$optimal@optimalTx
      } else if (attr(policy, "class") %in% c("GKLM")) { #### new function!
        x = output[as.logical(at.risk),,] %>% output2observable()
        action[at.risk != 0] <- 
          predict.opt.sep(policy[[stage]],  newdata = x, Tx.label = paste0("A_", stage))$optimal.Tx %>% 
          {ifelse(.==1, -1, 1)}
      } else if (attr(policy, "class") %in% c("GKRF")) { #### new function!
        x = output[as.logical(at.risk),,] %>% output2observable()
        x = x[, c(paste0("A_", stage), policy[[stage]][[1]]$xvar.names)]
        x[, paste0("A_", stage)] = 0 # dummy values
        
        action[at.risk != 0] <- 
          predict.opt.sep(policy[[stage]],  newdata = x, Tx.label = paste0("A_", stage))$optimal.Tx %>% 
          {ifelse(.==1, -1, 1)}
      }
      
      # if (stage == 2) stop("")
      action[at.risk == 0] <- NA
    } else {
      propensity = predPropensityFn(covariate = covariate)
      action = suppressWarnings(rbinom(n.sample, 1, propensity) * 2 - 1) # 1 for aggressive and 0 for gentle
      # for NAs in propensity, the action is NA. Thus, the warnings are suppressed.
      action[at.risk == 0] <- NA
    }
    
    ## progress              (0 < t < T_k)
    pred.hazard = predHazardFn(action = action, covariate = covariate)
    tmp3 <<- list(action = action, 
                  covariate = covariate
                  )
    # print(92)
    # tmp3 <<- list(time.max = time.max, tick = tick, rho.plus = rho.plus, omega.plus = omega.plus, 
    #               pred.hazard = pred.hazard, pred.censor = pred.censor, at.risk = at.risk,
    #               corr = 0.3, admin.censor = tau - surv.previous, 
    #               full_data = ifelse(hidden_data, 2, 0),
    #               terminal.stage = (stage == n.stages))
    stage.output <- #this creates event.time
      dynamics.vec(time.max = time.max, tick = tick, 
                   pred.hazard = pred.hazard, 
                   full_data = ifelse(hidden_data, 2, 0)
                   ) %>% t
    
    # tmp3 <<- stage.output
    # print(104)      
    ## bookkeeping and reassigning
    output[, c(colnames(tmp), "action", "at.risk"), stage] <- cbind(stage.output, action, at.risk)
    at.risk <- (stage.output[, "gamma"] == 0)  # Only those with gamma == 0 is at risk. Those with gamma = NA or 1 are not available.
    at.risk[is.na(at.risk)] = 0

    }
    attr(output, "tau") = tau
    attr(output, "p") = p
    print(n.stages)
    if (n.stages >1){
      if (summary) {
        cum.event.time = apply(output[, "event.time",], 1, cumsum) %>% t
        output.summary <- 
          tibble(subj.id = output[, "subj.id", 1],
                 rep.id = output[, "rep.id", 1],
                 terminal.stage = apply(output[, "at.risk", ], 1, sum), 
                 cumulative.event.time = cum.event.time[cbind(1:n.sample, terminal.stage)],
                 censor.status = output[cbind(1:n.sample, "delta", terminal.stage)],
                 actions       = apply(output[, "action", ], 1, 
                                       function(s) gsub("NA", "*", paste0(s, collapse = ""))))
        return(list(output = output, summary = output.summary))
      } else {
        return(output)
      }
    } else if (n.stages == 1){
      print("n.stages = 1")
      if (summary) {
        cum.event.time = output[, "event.time",] %>% t
        output.summary <- 
          tibble(subj.id = output[, "subj.id", 1],
                 rep.id = output[, "rep.id", 1],
                 terminal.stage = output[, "at.risk", ], 
                 cumulative.event.time = output[, "event.time",],
                 censor.status = output[,"delta",],
                 actions       = output[, "action", ])
        return(list(output = output, summary = output.summary))
      } else {
        return(output)
      }
    } else{
      print("ERROR: n.stages < 1")
    }
  }


if (0) {
  
  ## 0. basic settings
  n = 300
  p = 5
  
  beta.propensity1 = c(0, -1, 3, -1, rep(1, p))   # log(surv.prev + 1), rho, omega, covariate (1~5)
  beta.hazard1 = c(0, -1, rep(1, p))              # log(surv.prev + 1), covariate (1~5)
  # longer survived has less hazard, and milder treatment
  beta.censor1 = c(2, 1.5, rep(0.3, p))               # log(surv.prev + 1), covariate (1~5)
  # longer survived is more likely to be censored
  # predPropensityFn1 <- function(surv.previous, rho, omega, covariate) {
  #   plogis(cbind(1, log(surv.previous + 1), rho, omega, covariate) %*% beta.propensity1)
  # }
  # predHazardFn1 <- function(surv.previous, rho, omega, covariate) {
  #   # rho, omega are there for compatibility only.
  #   cbind(1, log(surv.previous + 1), covariate) %*% beta.hazard1
  # }
  # predCensorFn1 <- function(surv.previous, rho, omega, action, covariate) {
  #   # rho, omega are there for compatibility only.
  #   cbind(1, log(surv.previous + 1), covariate) %*% beta.censor1 * ifelse(action > 0, 1, 0.5) - 5  # trt interaction
  # }
  # noCensorFn1 <- function(surv.previous, rho, omega, action, covariate) {
  #   matrix(-10, nrow = length(surv.previous), ncol = 1)
  # }
  # predCensorFn <- predHazardFn <- predPropensityFn <- noCensorFn <- vector("list", 3)
  # for (i in 1:3) {
  #   predPropensityFn[[i]] <- predPropensityFn1
  #   predHazardFn[[i]]     <- predHazardFn1
  #   predCensorFn[[i]]     <- predCensorFn1
  #   noCensorFn[[i]]       <- noCensorFn1
  # }
  # 
  
  ## initial state         (t=0)
  set.seed(1)
  a.tmp <- 
    multiStageDynamics (
      at.risk = 1,      # initial state vector
      n.sample = n, n.stages = 3, tau = 10, tick = 0.01,  # structural parameters
      p = p, corr = -0.5,                              # cor of two error processes
      predHazardFn = predHazardFn, predPropensityFn = predPropensityFn, 
      predCensorFn = predCensorFn,             # list of predictor functions
      hidden_data = TRUE)                     # output control
  a.tmp$output
  a.tmp$summary
  1 - mean(a.tmp$summary$censor.status)
  a.tmp$summary$terminal.stage %>% table / n
  a.tmp$output[1:10, c("failure.time", "treatment.time", "censor.time"), 1]
  a.tmp$output[, c("failure.time"), 1] %>% hist
  a.tmp$output[, c("treatment.time"), 1] %>% hist
  a.tmp$output[, c("censor.time"), 1] %>% hist
  ggplot(a.tmp$summary, aes(terminal.stage, cumulative.event.time, col = actions)) +
    geom_point()
}


