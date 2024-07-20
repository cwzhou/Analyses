saving_dataset = FALSE
# LINE 126: trying to get more people in open group.

# smb://ad.unc.edu/med/tracs/groups/research/CDWH/McGinigle_Kate_IRB18-1153/Analysis/4.GrCURRENT_MRN_creation/grCURRENT_MRN2.xlsx
# setwd("/Volumes/med/tracs/Groups/Research/CDWH/McGinigle_Kate_IRB18-1153")
setwd("/Volumes/McGinigle_Kate_IRB18-1153")

### 0. library
library(readxl)
library(tidyverse)
library(dplyr)

### 1. Raw data read-in
proc0 <- read_excel("Analysis/4.Grid_creation/procedureLat.xlsx")
death0 <- read_excel("Analysis/5.Outcome_creation/mortality.xlsx")
comor0 <- read_excel("Analysis/4.Grid_creation/grid2.xlsx")
demo0 <- read_excel("_IRB18-1153_PROJECT_DATA/IRB18-1153_MASTER_PATIENT_WITH_ATA_PTA_TOE_PRESSURE.xlsx")
colnames(demo0) <- demo0[3, ]
demo0 = demo0[-(1:3), ]
rownames(demo0) <- NULL

# Dimension of the raw data
dim(proc0) # 2420 x 11
dim(death0) # 601 x 15
dim(comor0) # 601 x 44
dim(demo0) # 2029 x 20

### 2. Variables

# procedure codes (trt group and major amputations)
proc1 <- proc0 %>%
  # only keep subjects in either treatment group or who had major amputation
  filter(procGroup == "Open" | procGroup == "Endovascular" | procGroup == "Major amputation") %>%
  dplyr::select(-c(PROC_CODE, episode, woundClass, PROC_NAME, infectionClass,
                   DAYS_FROM_CLI_DX, CLI_DX_DATE)) %>%
  arrange(CURRENT_MRN, PROCEDURE_DATE)

# comorbidities
comor1 <- comor0 %>%
  arrange(CURRENT_MRN) %>%
  # if value is -50 then set equal to 0
  mutate(InflowDisease = ifelse(InflowDisease == -50, 0, InflowDisease),
         OutflowDisease = ifelse(OutflowDisease == -50, 0, OutflowDisease),
         RunoffDisease = ifelse(RunoffDisease == -50, 0, RunoffDisease),
  ) %>%
  dplyr::select(-c(hemodynamicLaterality, ptSide, STUDY_DATE...5, STUDY_DATE...44))
colnames(comor1)

# demographics
demo1 = demo0 %>%
  dplyr::select(ORIGINAL_FILE_MRN, CURRENT_MRN, MRN_SOURCE,
                # BIRTH_DATE, EPIC_DEATH_DATE, STATE_DEATH_DATA_DATE,
                GENDER, ETHNICITY, PATIENT_RACE_1) %>%
  filter(MRN_SOURCE == "CURRENT_MRN_IS_SAME") %>%
  dplyr::select(-c(ORIGINAL_FILE_MRN, MRN_SOURCE)) %>%
  arrange(CURRENT_MRN)

# death
death1 = death0 %>%
  dplyr::select(CURRENT_MRN,
                BIRTH_DATE,
                # EPIC_DEATH_DATE, STATE_DEATH_DATA_DATE,
                deathDate) %>%
  mutate(CURRENT_MRN = as.character(sprintf("%012s", CURRENT_MRN))) %>%
  arrange(CURRENT_MRN, deathDate) %>%
  mutate(BIRTH_DATE = format(BIRTH_DATE, "%Y-%m-%d"),
         # EPIC_DEATH_DATE = format(EPIC_DEATH_DATE, "%Y-%m-%d"),
         # STATE_DEATH_DATA_DATE = format(STATE_DEATH_DATA_DATE, "%Y-%m-%d"),
         deathDate = format(deathDate, "%Y-%m-%d")
  )


### 2B. description of the key variables
## (1) proc
# ID
# C7_FUTIME : follow-up time for MI or fatal CHD
# FUTIMED   : follow-up time for death (longer than or equal to C7_FUTIME)
# DEAD19    : status of death as of 2019-12-31
# UCOD      : underlying cause of death (coded)
# C7_SOURCINC : the source of event when C7_INC_BY19=1  ‘MI’ if the event is definite or probable MI
# C7_SOURCIP : the source of event when C7_IN_BY19P=1   ‘MI’ (‘FATCHD’) if there is no cardiac procedure or cardiac procedure is after the MI.
# C7_SOURCISP : the source of event when C7_IN_19SP=1   ‘MI’ (‘FATCHD’) if the event is definite or probable MI (definite fatal CHD)
# PREVHF01  : prevalent heart failure at visit 1.
# PRVCHD05  : prevalent coronary heart disease at visit 1.
# TIAB01    : prevalent stroke at visit 1.

## (2) death
# CENTER    : center four centers
# KNWNDEADBYVISIT21 : known to be dead at visit 2 => delta_1
# KNWNDEADBYVISIT31 : known to be dead at visit 3 => delta_2
# KNWNDEADBYVISIT41 : known to be dead at visit 4 ...
# KNWNDEADBYVISIT51 : known to be dead at visit 5
# KNWNDEADBYVISIT61 : known to be dead at visit 6
# KNWNDEADBYVISIT71 : known to be dead at visit 7 => delta_6


## (3) comor

## (4) demo
# ID, CENTER,
# ASPIRINCODE01   : Aspirin use in the past 2 weeks based on 2004 medication codes
# STATINCODE01    : Statin use in the past 2 weeks based on 2004 medication codes
# ANTICOAGCODE01  : Used Anticoagulates (At Visit 1) Last 2 Weeks (0=no, 1=yes) Based On 2004 Med Code
# PREVHF01, PRVCHD05 (overlaps with inc)
# RACEGRP, GENDER, V1AGE01
# V1DATE01,
# CIGTYR01    : Cigarette years of smoking
# BMI01       :
# WSTHPR01    : Waist-To-Hip Ratio
# DRNKR01     : Drinker Status (1: current drinker, 2: former, 3: never, 4: unknown)
# HYPERT04    : Hypertension, definition 4; replaces HYPERT01
# GLUCOS01
#
# GLEFH01 FAST1202 TCHSIU01 WORK_I02 SPRT_I02 ELEVEL01 ELEVEL02 HYPTMD01 ...
# OCCUPN01 MENOPS01 INTPLQ01 HYPTMDCODE01 MOMHISTORYSTR DADHISTORYSTR MOMHISTORYCHD DADHISTORYCHD MOMHISTORYDIA DADHISTORYDIA


### 3. Remove IDs with more than one amputable limb
dups_death = death1$CURRENT_MRN[duplicated(death1$CURRENT_MRN)]; length(dups_death)
dups_comor = comor1$CURRENT_MRN[duplicated(comor1$CURRENT_MRN)]; length(dups_comor)
dups_demo = demo1$CURRENT_MRN[duplicated(demo1$CURRENT_MRN)]; length(dups_demo)

# remove people who have duplicates aka multiple leg segments
comor2 = comor1 %>% filter(! CURRENT_MRN %in% dups_comor) %>% arrange(CURRENT_MRN)
death2 = death1 %>% filter(! CURRENT_MRN %in% dups_death) %>% arrange(CURRENT_MRN)
# sum(duplicated(comor2$CURRENT_MRN)) # want this to be 0
if (all(dups_comor == dups_death)){
  dupes = dups_comor
} else{
  dups = unique(c(dups_comor, dups_death))
}
# only select the people who do not have multiple limb options
proc2.0 = proc1 %>% filter(! CURRENT_MRN %in% dupes) %>% arrange(CURRENT_MRN)
demo2 = demo1 %>% filter(! CURRENT_MRN %in% dupes) %>% arrange(CURRENT_MRN)
# this includes major amputation laterality which uses proc2.0 to be created.
proc2_majoramplat = read_excel("Analysis/4.Grid_creation/major_amputation_laterality.xlsx") %>%
  mutate(CURRENT_MRN = as.character(sprintf("%012s", CURRENT_MRN)))
# check that they are the same
checking_procdate_proc2 = all(proc2.0$PROCEDURE_DATE == proc2_majoramplat$PROCEDURE_DATE)
if (checking_procdate_proc2){
  proc2 = proc2_majoramplat %>% dplyr::select(-Notes)
} else{
  stop("ERROR: PROC2.0 AND PROC2_MAJORAMPLAT ARE NOT THE SAME")
}

# checking IDs for comor and death and make sure they're the same.
matching_CURRENT_MRNs <- death2$CURRENT_MRN %in% comor2$CURRENT_MRN
# Check if there are any matching CURRENT_MRNs
all(matching_CURRENT_MRNs) # want to be true
mean(matching_CURRENT_MRNs) # want to be 1
dim(comor2);dim(death2);dim(proc2);dim(demo2)
# 543; 543; 1094; 1975
length(unique(proc2$CURRENT_MRN)) #349

### 4. Creating Treatment Groups and Identifying Event Dates
# modifying procedure data
proc2.5 = proc2 %>%
  # only keep subjects who have amputations on the diseased limb
  mutate(ind = if_else(is.na(`Confirmed Lat`) | laterality == `Confirmed Lat`, 1, 0)) %>%
  filter(ind == 1) %>%
  dplyr::select(CURRENT_MRN, PROCEDURE_DATE, procGroup, laterality)

proc3 <- proc2.5 %>%
  mutate(open_date = case_when(
    procGroup == "Open" ~ PROCEDURE_DATE,
    TRUE ~ NA
  ),
  endo_date = case_when(
    procGroup == "Endovascular" ~ PROCEDURE_DATE,
    TRUE ~ NA
  ),
  major_amp_date = case_when(
    procGroup == "Major amputation" ~ PROCEDURE_DATE,
    TRUE ~ NA
  )
  ) %>%
  dplyr::select(-c(PROCEDURE_DATE, procGroup))

proc4 <- proc3 %>%
  group_by(CURRENT_MRN) %>%
  # below assumes we just take the first date to occur
  summarize_all(.funs = list(~first(na.omit(.)))) %>%
  # only keep rows with missing amputation date, or where amputation is on or AFTER each treatment procedure date
  filter(is.na(major_amp_date) | major_amp_date >= open_date | major_amp_date >= endo_date) %>%
  mutate(open_ind0 = ifelse(!is.na(open_date) & is.na(endo_date), 1,
                            ifelse(!is.na(open_date) & !is.na(endo_date) & open_date < endo_date, 1, 0)),
         endo_ind0 = ifelse(!is.na(endo_date) & is.na(open_date), 1,
                            ifelse(!is.na(open_date) & !is.na(endo_date) & open_date >= endo_date, 1, 0))) %>%
  # CURRENTLY, ASSIGNING OPEN IF IT OCCURS 1 WEEK OR LESS AFTER ENDO!!!! OTHERWISE INTENTION-TO-TREAT (doing this so we get more people in OPEN treatment group)
  mutate(endo_week = endo_date + weeks(1)) %>%
  mutate(open2 = ifelse(is.na(endo_week), 1,
                        ifelse(is.na(open_date), 0,
                               as.integer(open_date <= endo_week)))) %>%
  mutate(open_ind = ifelse(open_ind0 == 1 | open2 == 1, 1, 0),
         # if open2 is 1 then we change that person's endo to be 0
         endo_ind = ifelse(open2 == 1, 0, endo_ind0))
dim(proc4) #325 x 10
# View(proc4)
# proc4 %>% mutate(sum = open_ind+endo_ind) %>% dplyr::select(sum) %>% filter(sum != 1) # should always be 0

proc5 <- proc4 %>%
  mutate(Trt = ifelse(open_ind == 1, "open", ifelse(endo_ind == 1, "endo", NA)),
         trt_date = pmin(open_date,endo_date,na.rm = T)) %>%
  dplyr::select(-c(open_date,endo_date,endo_week,
                   open_ind0,endo_ind0,open2,
                   open_ind,endo_ind)) %>%
  mutate(end_date = trt_date + years(1)) %>%
  mutate(amp_ind_tmp = ifelse(!is.na(major_amp_date) & major_amp_date <= end_date, 1, 0))
dim(proc5) #325 x 11
table(proc5$Trt) # endo: 252; open: 73

# modifying patients for other datasets
proc_ids = unique(proc5$CURRENT_MRN);length(proc_ids)
comor3 = comor2 %>%
  filter(CURRENT_MRN %in% proc_ids)
death3 = death2 %>%
  filter(CURRENT_MRN %in% proc_ids)
demo3 = demo2 %>%
  filter(CURRENT_MRN %in% proc_ids)
# View(demo3);View(proc5);View(comor3);View(death3)
dim(demo3); length(unique(proc5$CURRENT_MRN)); dim(comor3); dim(death3);
# now we end up with 320,325,325,325

### 5. Combining Datasets
dat1 <- left_join(death3, comor3, by = "CURRENT_MRN")
dat2 <- left_join(dat1, demo3, by = "CURRENT_MRN"); dim(dat2)
dat3 <- left_join(proc5, dat2, by = "CURRENT_MRN")
dat4 = dat3 %>%
  mutate(birth_date = as.Date(BIRTH_DATE),
         trt_date = as.Date(trt_date)) %>%
  # CALCULATE AGE AT TREATMENT
  mutate(age_at_treatment = as.numeric(difftime(trt_date,
                                                birth_date,
                                                units = "days") / 365.25)) %>%
  dplyr::select(-c(BIRTH_DATE, birth_date)) %>%
  mutate(death_ind_tmp = ifelse(is.na(deathDate),
                                0,
                                ifelse(deathDate >= trt_date & deathDate <= end_date,
                                       1,
                                       0)))  %>%
  mutate(event_tmp = ifelse(death_ind_tmp == 1 & amp_ind_tmp == 0, 1,
                            ifelse(death_ind_tmp == 0 & amp_ind_tmp == 1, 2,
                                   ifelse(death_ind_tmp == 0 & amp_ind_tmp == 0, 0, 99)))) %>%
  mutate(event_tmp = ifelse(event_tmp == 99,
                            ifelse(deathDate <= major_amp_date, 1, 2),
                            event_tmp)) %>%
  mutate(observed_time_tmp = case_when(
               event_tmp == 1 ~ as.numeric(difftime(deathDate, trt_date, units = "days")),
               event_tmp == 2 ~ as.numeric(difftime(major_amp_date, trt_date, units = "days")),
               event_tmp == 0 ~ as.numeric(difftime(end_date, trt_date, units = "days")),
               TRUE ~ NA_real_  # If none of the conditions are met, assign NA
             )) %>%
  mutate(ID = 1:nrow(dat3),
         obs_time = as.numeric(observed_time_tmp),
         status = as.numeric(event_tmp),
         D.0 = ifelse(event_tmp == 1 | event_tmp == 2, 1, 0),
         D.1 = ifelse(event_tmp == 1, 1, 0),
         D.2 = ifelse(event_tmp == 2, 1, 0)) %>%
  # choose the order i want
  dplyr::select(ID, Trt, trt_date, end_date,
                obs_time, status,
                event_tmp, observed_time_tmp,
                death_ind_tmp, deathDate, D.0, D.1,
                amp_ind_tmp, major_amp_date, D.2,
                everything()) %>%
  # remove the unneeded variables
  dplyr::select(-c(CURRENT_MRN, major_amp_date, deathDate,
                   trt_date, end_date,
                   observed_time_tmp, event_tmp,
                   death_ind_tmp, amp_ind_tmp,
                   laterality,
                   absTimeToImageStudy))

### 6. Covariates
# clean up covariates
dat5 = dat4 %>%
  mutate(ethnicity = ifelse(ETHNICITY == "UNKNOWN", NA, ETHNICITY),
         race = ifelse(PATIENT_RACE_1 == "UNKNOWN", NA,
                       ifelse(PATIENT_RACE_1 == "WHITE OR CAUCASIAN", "WHITE OR CAUCASIAN",
                              ifelse(PATIENT_RACE_1 == "BLACK OR AFRICAN AMERICAN", "BLACK OR AFRICAN AMERICAN",
                                     "OTHER"
                                     )))
         ) %>%
  dplyr::select(-c(ETHNICITY, PATIENT_RACE_1))
# View(dat5)
table(dat5$Trt)
table(dat5$status) #censored: 231; death: 38; major amputation: 56
table(dat5$GENDER) # female: 131; male: 189
table(dat5$race)
table(dat4$PATIENT_RACE_1)
table(dat5$ethnicity)

df_tmp0 = dat5 %>%
  mutate(Diabetes = ifelse(DIABETES == 1 | `DIABETES-COMPLICATED` == 1 | `DIABETES-UNCOMPLICATED` == 1, 1, 0),
         Hypertension = ifelse(`HYPERTENSION-COMPLICATED` == 1 | `HYPERTENSION-UNCOMPLICATED` == 1, 1, 0),
         Renal = ifelse(`RENAL DISEASE-CKD` == 1 | `RENAL DISEASE-COMPLICATED` == 1 | `RENAL DISEASE-ESRD` == 1| RENAL_DISEASE == 1, 1, 0),
         Hyperlipidemia = ifelse(`HYPERLIPIDEMIA-COMPLICATED` == 1 | `HYPERLIPIDEMIA-NA` == 1, 1, 0),
         Wound_Comorbidity = ifelse(`WOUND-COMPLICATED` == 1 | `WOUND-EXTREMITY` == 1 | WOUND  == 1, 1, 0),
         Cerebrovascular = ifelse(`CEREBROVASCULAR DISEASE-COMPLICATED` == 1 | CEREBROVASCULAR == 1, 1, 0)) %>%
  dplyr::select(-c(DIABETES, `DIABETES-COMPLICATED`, `DIABETES-UNCOMPLICATED`,
                   `HYPERTENSION-COMPLICATED`, `HYPERTENSION-UNCOMPLICATED`,
                   `RENAL DISEASE-CKD`, `RENAL DISEASE-COMPLICATED`, `RENAL DISEASE-ESRD`, RENAL_DISEASE,
                   `HYPERLIPIDEMIA-COMPLICATED`, `HYPERLIPIDEMIA-NA`,
                   `WOUND-COMPLICATED`, `WOUND-EXTREMITY`, WOUND,
                   `CEREBROVASCULAR DISEASE-COMPLICATED`, CEREBROVASCULAR
  )) %>%
  mutate(Male = ifelse(GENDER == "MALE", 1, 0), # 1 = male; 0 = female
         Ethnicity = ifelse(ethnicity == "NOT HISPANIC OR LATINO", 1, 0) # 1 = "NOT HISPANIC OR LATINO"; 0 = "HISPANIC OR LATINO"
  ) %>%
  mutate(ID = factor(ID),
         Trt = factor(Trt, levels = c("open", "endo"), labels = c("Open", "Endovascular")),
         woundClass = factor(woundClass, levels = c(0,1,2,3)),
         maxRutherfordClass = factor(maxRutherfordClass, levels = c(4,5,6)),
         Race = factor(race, levels = c("BLACK OR AFRICAN AMERICAN",
                                        "WHITE OR CAUCASIAN",
                                        "OTHER"),
                       labels = c("BLACK OR AFRICAN AMERICAN",
                                  "WHITE OR CAUCASIAN",
                                  "OTHER RACE")),
         ischemia = factor(ischemia, levels = c(0,1,2,3))
  ) %>%
  # dplyr::select(-c(GENDER,
  #                  race,
  #                  ethnicity)) %>%
  # # below these have values of 0
  # dplyr::select(-c(ANEMIA,
  #                  infectionClass)) %>%
  as.data.frame()

# Get column names into a vector
column_names <- colnames(df_tmp0)
# Convert column names into a comma-separated string
column_names_string <- paste(column_names, collapse = ", ")
# Print the resulting string
print(column_names_string)

df_tmp1 = df_tmp0 %>%
  dplyr::select(ID, # factor
                Trt, # factor
                obs_time, status, D.0, D.1, D.2,
                ischemia, #(4-level factor)
                woundClass, #(4-level factor)
                maxRutherfordClass, #(3-level factor)
                Race, #(3-level factor)
                Ethnicity, #indicator
                Male, #indicator
                InflowDisease, OutflowDisease, RunoffDisease,
                stentBPG, CAD,
                CHRONIC_PULM_DISEASE, COAGULOPATHY,
                `CONGESTIVE HEART FAILURE-COMPLICATED`, CHF,
                `DEFICIENCY ANEMIA-NA`, `DEMENTIA-NA`, MI, OBESITY,
                SMOKING, VENOUS_INSUFFICIENCY, WEIGHT_LOSS, VTE,
                age_at_treatment,
                Diabetes, Hypertension, Renal, Hyperlipidemia,
                Wound_Comorbidity, Cerebrovascular,
                ); dim(df_tmp1)
# convert all indicators to factors except D.0,D.1,D.2
df_tmp <- df_tmp1 %>%
  mutate_at(vars(-one_of("D.0", "D.1", "D.2")), ~ if(all(n_distinct(na.omit(.)) == 2)) as.factor(.) else .)

# below is for complete case which we aren't doing.
#length(unique(df_tmp$ID)) # N = 325 but 4 people don't have demographics
# which(! proc_ids %in% demo4$CURRENT_MRN) # the 4 people without demographics
# # Check for missing values in each column
# missing_values <- colSums(is.na(df_tmp))
# # Identify variables with missing values
# variables_with_missing <- names(missing_values[missing_values > 0])
# print(variables_with_missing)
# rows_with_missing <- !complete.cases(df_tmp)
# if (any(is.na(df_tmp))) {
#   print("There is missing data in the dataframe.")
# } else {
#   print("There is no missing data in the dataframe.")
# }
df_tmpcc = df_tmp[complete.cases(df_tmp), ] #303 people (311)
dim(df_tmpcc) #303 x 37
table(df_tmpcc$Trt) # endo: 236; open: 67
# Proportions within each level of Trt
prop.table(table(df_tmp$Trt, df_tmp$status), margin = 1)
prop.table(table(df_tmpcc$Trt, df_tmpcc$status), margin = 1)
# View(df_tmpcc)

# # see which vars are factors
# factor_vars <- pad_df %>%
#   select_if(is.factor) %>%
#   names()
# # Print the names of factor variables
# print(factor_vars)

# # seeing how many 1's in the indicators
# numeric_vars <- sapply(pad_df, is.numeric)
# # Apply sum() function only to numeric variables
# col_sums <- sapply(pad_df[, numeric_vars], sum)
# # Print the column sums
# print(col_sums)

pad_df = df_tmpcc;dim(pad_df)

if (saving_dataset == TRUE){
  print("saving dataset")
  write.csv(pad_df, file = "Analysis/pad_df.csv", row.names = FALSE)
}

