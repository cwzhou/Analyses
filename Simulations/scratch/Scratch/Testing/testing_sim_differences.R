df = as.data.frame(cbind(ID = rep1$subj.id, CZMK = rep1$event.time, TRT1 = rep1$action, ZOM = rep2$event.time, TRT2 = rep2$action))
climb <- df[df[,"CZMK"] != df[,"ZOM"], ]

climbczmk <- climb[climb[,"CZMK"] >= climb[,"ZOM"],]
climbzom <- climb[climb[,"CZMK"] < climb[,"ZOM"],]
check_ids = climbzom$ID
good_ids = climbczmk$ID

dim(df); dim(climb); dim(climbczmk); dim(climbzom)

check_ids1 = head(check_ids,20)
good_ids1 = head(good_ids)
for (id in check_ids1){
  idczmk1 = rep1 %>% filter(subj.id == id)
  idzom1 = rep2 %>% filter(subj.id == id)
  message("czmk")
  print(idczmk1)
  message("zom")
  print(idzom1)
}

for (id in good_ids1){
  idczmk1 = rep1 %>% filter(subj.id == id)
  idzom1 = rep2 %>% filter(subj.id == id)
  message("czmk")
  print(idczmk1)
  message("zom")
  print(idzom1)
}


ggplot(d, aes(y = mean, x = factor(sim, levels = c(1:10)), color = trt)) +
  geom_point() +
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 1) +
  facet_grid(~method) +
  ggtitle("Scatter Plot")

d$mean = as.numeric(d$mean)
d$within_10_percent <- (d$mean[d$trt == "A"] >= (1-0.1) * d$mean[d$trt == "B"])
  
ggplot(d, aes(y = mean, x = factor(sim, levels = c(1:10)), color = within_10_percent, shape = trt)) +
  geom_point() +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red")) + 
  scale_shape_manual(values = c("A" = 2, "B" = 0)) +  # 0 represents a square, 1 represents a circle
  facet_grid(~method) +
  ggtitle("Scatter Plot")



long_data = d 
library(ggplot2)

# Assuming 'long_data' is your long format dataset
# Create a sample long format dataset
set.seed(123)
long_data <- data.frame(
  subject = rep(1:50, each = 20),
  trt = rep(c("A", "B"), each = 500),
  sim = rep(1:10, each = 100, times = 2),
  mean = rnorm(1000)
)

# Calculate mean and confidence interval
summary_data <- aggregate(mean ~ trt + sim, data = long_data, FUN = function(x) c(mean = mean(x), ci = t.test(x)$conf.int))

# Plot means with error bars using ggplot2
ggplot(summary_data, aes(x = mean, color = trt)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = ci[1], ymax = ci[2], group = trt), position = position_dodge(width = 0.5), width = 0.2) +
  facet_grid(~trt) +
  ggtitle("Mean Plot with Error Bars by Treatment and Simulation")
