pooled1 = FALSE
nodesize = 100
criterion_phase1 = "mean"
t0_crit = 19#2200 # six-year survival
### 0.2 criterion
# for (criterion_phase1 in c("mean", "mean.prob.combo")){ # for (dtr_criterion in c("mean", "surv.mean")) {
if (criterion_phase1[1] != "mean") {
  criterion_phase1[2] = t0_crit
  dtr_criterion[2] = criterion_phase1[2]
  rule = "logrank"
} else {
  rule = "mean"
  criterion_phase1[2] = NA
}
if (criterion_phase2[1] != "mean") {
  criterion_phase2[2] = t0_crit
} else {
  criterion_phase2[2] = NA
}
if (is.null(criterion_phase1[2]) | is.na(criterion_phase1[2])){
  t0_pmcr = tau/2
  crit_tmp = ""
} else{
  t0_pmcr = criterion_phase1[2]
  crit_tmp = round(as.numeric(criterion_phase1[2]), 1)
}
criterion_phase2 = criterion_phase1
skip.PMCR = TRUE

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
                  pooled = pooled1,
                  tol1 = c(0.1,0),
                  stratifiedSplit = FALSE)# actual fitting
values[cv, "ns.CZMK"] = nodesize
# set.seed(cv)
CZMK.i <-
  try(do.call(itrSurv,
              c(args.CZMK, list(nodeSize = nodesize,
                                minEvent = mindeath))))
err.CZMK = class(CZMK.i)[1] == "try-error"
