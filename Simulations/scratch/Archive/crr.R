# simulated data to test 
set.seed(10)
ftime <- rexp(200)
fstatus <- sample(0:2,200,replace=TRUE)
cov <- matrix(runif(600),nrow=200)
dimnames(cov)[[2]] <- c('x1','x2','x3')
trt = sample(0:1,200,replace=TRUE)
d = cbind(id = 1:200, ftime,fstatus,cov)
View(as.data.frame(d) %>%
       filter(fstatus == 1) %>%
       dplyr::select(ftime) %>%
       distinct() %>%
       arrange(ftime))

print(z <- crr(ftime,fstatus,cov))
summary(z)
z.p <- predict(z,
               rbind(c(.1,.5,.8),
                     c(.1,.5,.2)
                     )); head(z.p)
dim(z.p) # unique failure times from cause 1 x ncov
# z.p[,1] is unique failure times from cause 1 in dataset
# z.p[,2] is x
# Returns a matrix with the unique type 1 failure times 
# in the first column, and the other columns giving 
# the estimated subdistribution function 
# corresponding to the covariate combinations 
# in the rows of cov1 and cov2, 
# at each failure time 
# (the value that the estimate jumps to 
#   at that failure time).

#plots the subdistribution functions 
# estimated by predict.crr, 
# by default using a different line type for each curve
plot(z.p) 


z2 = crr(ftime,fstatus,cov,failcode=2)
z.p2 <- predict(z2,
               rbind(c(.1,.5,.8),
                     c(.1,.5,.2)))
plot(z.p2,lty=1,color=2:3)
