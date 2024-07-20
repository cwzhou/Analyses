n.sample = 3; n.stages=1; tau = 10; tick = 0.01;          # structural parameters
rho = NULL; omega = NULL; surv.previous = rep(0, n.sample); # initial state vector
covariate = NULL; at.risk = rep(1, n.sample);             
p = 5; Sig = diag(p) + 0.2 - diag(p) * 0.2;               # covariate structure
corr = -0.5;                                               # cor of two error processes
predHazardFn; predPropensityFn; predCensorFn;             # list of predictor functions
hidden_data = FALSE; summary = TRUE;                      # output control
policy = NULL;                                            # optimal rule (if !is.null, propensity scores are ignored.) for value calculation
printFlag = TRUE   
