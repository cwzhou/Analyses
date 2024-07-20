require(dplyr)
if (length(at.risk) == 1) at.risk = rep(at.risk, n.sample)
if (length(surv.previous) != n.sample)  stop ("length of surv.previous should match.")
if (length(at.risk) != n.sample)  stop ("length of at.risk should match.")
if (!is.null(omega)) if (length(omega) != n.sample)  stop ("length of omega should match.")
if (!is.null(covariate)) if (dim(covariate)[1] != n.sample) stop ("length of at.risk should match.")
if (!is.null(rho)) if(length(rho) != n.sample)  stop ("length of covariate should match.")

time.max = max(tau - surv.previous)

# initial state
if (is.null(rho))           rho = runif(n.sample) * 0.5 + 0.5
# if (is.null(rho))           rho = rep(1, n.sample)
if (is.null(omega))         omega = runif(n.sample) * 0.5 + 0.5
if (is.null(surv.previous)) surv.previous = rep(0, n.sample)
if (is.null(covariate))     covariate = mvrnorm(n.sample, rep(0, p), Sig)
if (is.null(colnames(covariate))) colnames(covariate) = paste0("Z", 1:p)
# skeleton
tmp <-
  dynamics.vec(time.max = 1, tick = tick, rho.plus = rep(0.5, 2), omega.plus = rep(0.5, 2), 
               pred.hazard = 0, pred.censor = 0, admin.censor = Inf,
               at.risk = 1, corr = corr, full_data = ifelse(hidden_data, 2, 0)) %>% t

output <- 
  array(NA, dim = c(n.sample, 2 + dim(tmp)[2] + 6 + p, n.stages), 
        dimnames = list(1:n.sample, c("subj.id", "rep.id", colnames(tmp), "action", "at.risk", 
                                      "surv.previous", "lB", "rho.0", "omega.0", colnames(covariate)), 
                        1:n.stages))
output[, "subj.id", ] <- 1:n.sample
output[, "rep.id", ] <- 0           # rep.id is reserved for later use (repeated random draws)
action = rep(NA, n.sample)           # initialize action vector

# tmp.tmp <<- list()    
for (stage in 1:n.stages) {
  # print(paste0("stage ", stage))
  if (printFlag) cat("stage ", stage, "\n")
  output[, "at.risk", stage]           <- at.risk
  output[, "surv.previous", stage]     <- surv.previous
  output[, "lB", stage]                <- log(surv.previous + 1)
  output[, colnames(covariate), stage] <- covariate
  output[, "rho.0", stage]             <- rho
  output[, "omega.0", stage]           <- omega
  # tmp.out <<- output 
  ## action                (t=0)
  if (!is.null(policy)) {
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
    print("hi")
    propensity = predPropensityFn[[stage]](surv.previous = surv.previous, rho = rho, 
                                           omega = omega, covariate = covariate)
    action = suppressWarnings(rbinom(n.sample, 1, propensity) * 2 - 1) # 1 for aggressive and 0 for gentle
    # for NAs in propensity, the action is NA. Thus, the warnings are suppressed.
    action[at.risk == 0] <- NA
  }
  
  ## instantaneous state   (t=0+)
  rho.plus = rho / omega / ifelse (action > 0, 10, 4)
  omega.plus = omega - ifelse(action > 0, 0.5, 0.25)
  
  ## progress              (0 < t < T_k)
  pred.hazard = predHazardFn[[stage]](surv.previous = surv.previous, rho = rho, 
                                      omega = omega, action = action, covariate = covariate)
  pred.censor = predCensorFn[[stage]](surv.previous = surv.previous, rho = rho, 
                                      omega = omega, action = action, covariate = covariate)
  tmp3 <<- list(surv.previous = surv.previous, rho = rho, 
                omega = omega, action = action, covariate = covariate,
                pred.censor = pred.censor,
                cens.fn = predCensorFn[[stage]])
  # print(92)
  # tmp3 <<- list(time.max = time.max, tick = tick, rho.plus = rho.plus, omega.plus = omega.plus, 
  #               pred.hazard = pred.hazard, pred.censor = pred.censor, at.risk = at.risk,
  #               corr = 0.3, admin.censor = tau - surv.previous, 
  #               full_data = ifelse(hidden_data, 2, 0),
  #               terminal.stage = (stage == n.stages))
  stage.output <- 
    dynamics.vec(time.max = time.max, tick = tick, rho.plus = rho.plus, omega.plus = omega.plus, 
                 pred.hazard = pred.hazard, pred.censor = pred.censor, at.risk = at.risk,
                 corr = corr, admin.censor = tau - surv.previous, 
                 full_data = ifelse(hidden_data, 2, 0),
                 terminal.stage = (stage == n.stages)) %>% t
  # tmp3 <<- stage.output
  # print(104)      
  ## bookkeeping and reassigning
  output[, c(colnames(tmp), "action", "at.risk"), stage] <- cbind(stage.output, action, at.risk)
  at.risk <- (stage.output[, "gamma"] == 0)  # Only those with gamma == 0 is at risk. Those with gamma = NA or 1 are not available.
  at.risk[is.na(at.risk)] = 0
  
  rho   <- stage.output[, "rho.x"  ]
  omega <- stage.output[, "omega.x"]
  
  ## covariate for the next stage. (actually not needed for the terminal stage)
  covariate = covariate * 0.5 + mvrnorm(n.sample, rep(0, p), Sig) * 0.5
  covariate[, 1] <- covariate[, 1] * 0.2 + rho * 0.8
  covariate[, 2] <- covariate[, 2] * 0.2 + omega * 0.8
  surv.previous <- surv.previous + stage.output[, "event.time"] # updating for the next stage
}
attr(output, "tau") = tau
attr(output, "p") = p
if (summary) {
  print(output)
  if (n.stages == 1) {
    cum.event.time <- output[, "event.time",]
    terminal.stage = output[,"at.risk",]
    cumulative.event.time = cum.event.time[cbind(1:n.sample, terminal.stage)]
  } else {
    cum.event.time <- t(apply(output[, "event.time",], 1, cumsum))
    terminal.stage = apply(output[, "at.risk", ], 1, sum)
    cumulative.event.time = cum.event.time[cbind(1:n.sample, terminal.stage)]
  }
  print(cum.event.time)
  print(terminal.stage)
  output.summary <- 
    tibble(subj.id = output[, "subj.id", 1],
           rep.id = output[, "rep.id", 1],
           terminal.stage = terminal.stage, 
           cumulative.event.time = cumulative.event.time,
           censor.status = output[cbind(1:n.sample, "delta", terminal.stage)],
           actions       = apply(output[, "action", ], 1, 
                                 function(s) gsub("NA", "*", paste0(s, collapse = ""))))
  return(list(output = output, summary = output.summary))
} else {
  return(output)
}
