# OS means similar between treatments BUT
# trt0 better than trt1 (higher survival)
one = rep_obs %>% group_by(action) %>% summarize(mean = mean(event.time), n = n())
one1 = ifelse(one$mean[2]<one$mean[1],"TRUE","FALSE")

# cause1 better survival than cause2 (but lower in number aka more people with cause2) AND
two = rep_obs %>% group_by(status) %>% summarize(mean = mean(event.time), n = n()) %>%
  filter(status != 0)
# & two$n[1]<=two$n[2]
two2 = ifelse(two$mean[1]>=two$mean[2],"TRUE","FALSE")

# WITIHIN cause: trt0 < trt1 for cause1; trt0 > trt1 for cause2
three = rep_obs %>% group_by(status,action) %>% summarize(mean = mean(event.time), n = n()) %>% ungroup() %>%
  filter(status != 0)
three3 = ifelse(three$mean[1] <= three$mean[2] & three$mean[3] >= three$mean[4],
               "TRUE",
               "FALSE")


final = ifelse(one1 == "FALSE" | two2 == "FALSE" | three3 == "FALSE", "at least one is false", "all true")
if (final == "all true"){
  message("HOORAY CONDITIONS MET!")
} else{
  message(final)
  print(sprintf("OS means similar between treatments BUT trt0 better than trt1 (higher survival): %s",one1))
  print(sprintf("cause1 better survival than cause2 (but lower in number aka more people with cause2): %s",two2))
  print(sprintf("WITIHIN cause: trt0 < trt1 for cause1; trt0 > trt1 for cause2: %s",three3))
}
