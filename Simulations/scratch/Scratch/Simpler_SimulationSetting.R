# Generate 2 exponential RV which depend on covariates, 
# one for the competing event of interest $E_1$ 
# and one for other causes, $E_2$.If $E_1<E_2$, 
# then the individual dies of the event of interest at time $E_1$; 
# otherwise the person dies from the other event at time $E_2$.
# Hopefully get more direct control over the marginal cumulative incidence curves.

beta1 = list(beta1.hazard0 = c(0, 0.1 * c(1,1)),         # int, covariate (1~2) #cause1 treatment1
             beta1.hazard1 = c(0, 5 * c(1,1)),      # int, covariate (1~2) #cause1 treatment0
             beta2.hazard0 = c(0, 1 * c(1,1)),         # int, covariate (1~2)
             beta2.hazard1 = c(0, 1 * c(1,1))) 