train_truth = merge(train0.data, train1.data, by = "subj.id") %>%
  dplyr::select(-c(at.risk.x, at.risk.y, Z1.x,Z1.y,#Z2.x,Z2.y,
                   failure_t1.x, failure_t2.x, failure_t1.y, failure_t2.y)) %>%
  mutate(best.action = ifelse(event.time.x > event.time.y, action.x, action.y),
         best.status = ifelse(event.time.x > event.time.y, status.x, status.y),
         best.time = ifelse(event.time.x > event.time.y, event.time.x, event.time.y))  %>%
  dplyr::select(-c(action.x, action.y)) %>%
  mutate(obs.time0 = round(event.time.x,2),
         obs.time1 = round(event.time.y,2),
         status0 = status.x,
         status1 = status.y
  ) %>%
  dplyr::select(subj.id, best.time, best.action, best.status,
                obs.time0, obs.time1, status0, status1)

test_truth = merge(test0.data.rep, test1.data.rep, by = "subj.id") %>%
  dplyr::select(-c(at.risk.x, at.risk.y, Z1.x,Z1.y,#Z2.x,Z2.y,
                   failure_t1.x, failure_t2.x, failure_t1.y, failure_t2.y)) %>%
  mutate(best.action = ifelse(event.time.x > event.time.y, action.x, action.y),
         best.status = ifelse(event.time.x > event.time.y, status.x, status.y),
         best.time = ifelse(event.time.x > event.time.y, event.time.x, event.time.y))  %>%
  dplyr::select(-c(action.x, action.y)) %>%
  mutate(obs.time0 = round(event.time.x,2),
         obs.time1 = round(event.time.y,2),
         status0 = status.x,
         status1 = status.y
  ) %>%
  dplyr::select(subj.id, best.time, best.action, best.status,
                obs.time0, obs.time1, status0, status1)


# View(train_truth); View(test_truth)

message('overall survival')
train_truth %>%
  ungroup() %>%
  group_by(best.action) %>%
  summarise(mean = mean(best.time)) %>%
  print()

message('by cause')
train_truth %>%
  ungroup() %>%
  group_by(best.status,best.action) %>%
  summarise(mean = mean(best.time)) %>%
  print()

message('by cause')
test_truth %>%
  ungroup() %>%
  group_by(best.status,best.action) %>%
  summarise(mean = mean(best.time)) %>%
  print()

message("\nnumber of people with actions equaling true best.action")
message("czmk:",round(sum(rep_czmk$action == test_truth$best.action),2))
message("csk:",round(sum(rep_csk$action == test_truth$best.action),2))
message("zom:",round(sum(rep_zom$action == test_truth$best.action),2))
message("obs:",round(sum(rep_obs$action == test_truth$best.action),2))

checking = cbind(test_truth$best.time,
      rep_czmk$event.time,
      rep_csk$event.time,
      rep_zom$event.time,
      rep_obs$event.time,
      ifelse(test_truth$best.time==rep_czmk$event.time,1,0),
      ifelse(test_truth$best.time==rep_csk$event.time,1,0),
      ifelse(test_truth$best.time==rep_zom$event.time,1,0),
      ifelse(test_truth$best.time==rep_obs$event.time,1,0))
colnames(checking) = c("true", "czmk", "csk", "zom", "obs",
                       "ind_czmk", "ind_csk", "ind_zom", "ind_obs")
# View(checking)

message("\n\nproportion of people from CZMK with true best time: ",
        round(mean(checking[,"ind_czmk"]),2))
message("proportion of people from CSK with true best time: ",
        round(mean(checking[,"ind_csk"]),2))
message("proportion of people from ZOM with true best time: ",
        round(mean(checking[,"ind_zom"]),2))
message("proportion of people from OBS with true best time: ",
        round(mean(checking[,"ind_obs"]),2))

c0czmk = test_truth[(rep_czmk$action==test_truth$best.action),]
cczmk = rep_czmk[(rep_czmk$action==test_truth$best.action),]
c0csk = test_truth[(rep_csk$action==test_truth$best.action),]
ccsk = rep_csk[(rep_csk$action==test_truth$best.action),]
c0zom = test_truth[(rep_zom$action==test_truth$best.action),]
czom = rep_zom[(rep_zom$action==test_truth$best.action),]
c0obs = test_truth[(rep_obs$action==test_truth$best.action),]
cobs = rep_obs[(rep_obs$action==test_truth$best.action),]

# View(c0csk)
# View(ccsk)
# View(c0czmk)
# View(cczmk)

# Looking at how many MATCHING TO TRUTH actions overlap between CZMK and CSK
# Subject ID vectors
vector1 <- cczmk$subj.id
vector2 <- ccsk$subj.id
vector3 <- czom$subj.id
vector4 <- cobs$subj.id
length(vector1)
length(vector2)
length(vector3)
length(vector4)

# Find overlapping elements
overlap <- vector1 %in% vector2
overlapping_elements <- vector1[overlap]
print(overlapping_elements);message("length:",length(overlapping_elements))
message("overlap out of czmk:", round(length(overlapping_elements)/length(vector1),2))
message("overlap out of csk:", round(length(overlapping_elements)/length(vector2),2))


# looking at the mean over survival times when czmk action EQ or NEQ truth
mean(rep_czmk[(rep_czmk$action == test_truth$best.action),]$event.time)
# mean(test_truth[(rep_czmk$action == test_truth$best.action),]$best.time)
mean(rep_csk[(rep_czmk$action == test_truth$best.action),]$event.time)
mean(rep_zom[(rep_czmk$action == test_truth$best.action),]$event.time)

mean(rep_czmk[!(rep_czmk$action == test_truth$best.action),]$event.time)
mean(test_truth[!(rep_czmk$action == test_truth$best.action),]$best.time)
mean(rep_csk[!(rep_czmk$action == test_truth$best.action),]$event.time)
mean(rep_zom[!(rep_czmk$action == test_truth$best.action),]$event.time)
mean(rep_obs[!(rep_czmk$action == test_truth$best.action),]$event.time)
#
# look = rep_czmk[!(rep_czmk$action == test_truth$best.action),] %>%
#   dplyr::select(-at.risk) %>%
#   mutate(abs_diff = abs(failure_t1-failure_t2))
# range(look$abs_diff)
#
# look %>% filter(abs_diff > 100) %>% head()
source("Scratch/checking_prop_of_truth_testingset.R")
