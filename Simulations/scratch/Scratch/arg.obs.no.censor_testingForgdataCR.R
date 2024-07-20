# gdata_CR in functions

N = arg.obs.no.censor$N; ncov = arg.obs.no.censor$ncov; mass_p = arg.obs.no.censor$mass_p;
ctype = arg.obs.no.censor$ctype
cparam=2;censor_min=0;censor_max=NULL;
ztype = 0; zparam = 0.5; num_A = 2
tau = arg.obs.no.censor$tau; policy = NULL; seed1 = arg.obs.no.censor$seed1
u1 = arg.obs.no.censor$u1


# multiphasedynamicsCR
censor_time = trunc_cens_time
N = arg.obs.no.censor$N
at.risk = arg.obs.no.censor$at.risk; ncov = arg.obs.no.censor$ncov; corr = arg.obs.no.censor$corr;
at.risk = rep(1, N)
N = arg.obs.no.censor$N; n.phases = arg.obs.no.censor$n.phases; tau = arg.obs.no.censor$tau;
tick = arg.obs.no.censor$tick; hidden_data = arg.obs.no.censor$hidden_data; printFlag = arg.obs.no.censor$printFlag;
predHazardFn = arg.obs.no.censor$predHazardFn; predPropensityFn = arg.obs.no.censor$predPropensityFn;
predCensorFn = arg.obs.no.censor$predCensorFn; surv.previous = rep(0, N)
rho = NULL; omega = NULL; covariate = NULL;
Sig = diag(ncov) + 0.2 - diag(ncov) * 0.2
