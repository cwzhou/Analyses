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
  View(proc4)
  # proc4 %>% mutate(sum = open_ind+endo_ind) %>% dplyr::select(sum) %>% filter(sum != 1) # should always be 0

  proc5 <- proc4 %>%
    mutate(Trt = ifelse(open_ind == 1, "open", ifelse(endo_ind == 1, "endo", NA)),
           trt_date = pmin(open_date,endo_date,na.rm = T)) %>%
    dplyr::select(-c(open_date,endo_date,open_ind,endo_ind)) %>%
    mutate(end_date = trt_date + years(1)) %>%
    mutate(amp_ind_tmp = ifelse(!is.na(major_amp_date) & major_amp_date <= end_date, 1, 0))
  dim(proc5) #325 x 11
  table(proc5$Trt) # endo: 252; open: 73

  ## 9. Final set.
  dat.clean =
    outcome %>%
    ##    9.1 Removes kd.?, date.?, new.date.?, cum.?, next.fu, date.next.fu.
    dplyr::select(ID, T.0, DEAD19, UCOD, last.fu, d.1:d.7, V.1:V.7) %>%
    ##    9.2 Adding necessary information back
    left_join(inc %>% select(ID, PREVHF01, PRVCHD05, TIAB01), by = "ID") %>%
    left_join(sta %>% select(ID, CENTER), by = "ID") %>%
    left_join(visit1 %>% select(-date.1), by = "ID") %>%
    left_join(visit2 %>% select(-date.2), by = "ID") %>%
    left_join(visit3 %>% select(-date.3), by = "ID") %>%
    left_join(visit4 %>% select(-date.4), by = "ID") %>%
    left_join(visit5 %>% select(-date.5), by = "ID") %>%
    left_join(visit6 %>% select(-date.6), by = "ID") %>%
    left_join(visit7 %>% select(-date.7), by = "ID")


  saveRDS(dat.clean, "../data_ARIC/01.aric.clean.rds")



  ### STEP 2.
  ### Subsetting the latest visits only Visit 5 (2011 ~ 2012) and after. (Visits 4 were in 1996-1997, visits 1 were 1986-1988)
  all_cause = TRUE

  ### 1. trimmer(): stages 1:7 to 5:7
  trimmer = function(data, full = 1:7, reduced = 4:6,
                     stage.length = "V", status = "d",
                     sep = ".",
                     total.time = "T.0", last.fu = "last.fu") {
    require(dplyr)
    # reduced should be continuous (no intermittent missing).

    # 0. defining the variables
    vars = names(data)
    # reduced = reduced         # 4, 5, 6
    # to.be = 1:length(reduced) # 1, 2, 3

    # V.full = paste0(stage.length, sep, full)
    # d.full = paste0(status, sep, full)
    V.new = paste0(stage.length, sep, reduced)
    # d.new = paste0(status, sep, reduced)

    # if reduced does not end with the final stage, the last stage length is augmented.
    if (max(reduced) < max(full)) {
      V.residual = paste0(stage.length, sep, max(reduced) : max(full))     # if reduced = 2:4, obtain c("V.4", "V.5", "V.6", "V.7").
      data[, V.residual[1]] = apply(data[, V.residual], 1, sum, na.rm = TRUE)
    }

    # 1. updating the total time
    data[, total.time] = apply(data[, V.new], 1, sum, na.rm = TRUE)

    new.i = 0
    for (i in full) {
      if (i %in% reduced) {
        new.i = new.i + 1
        vars = gsub(paste0("\\.", i), paste0("\\.NEW", new.i), vars)
      } else {
        vars = gsub(paste0("\\.", i), paste0("\\.DEPRECATED", i), vars)
      }
    }
    names(data) = vars

    ## 2. Drop the unnecessary variables, and update the last.fu
    # drop the unused visits
    data = data %>%
      dplyr::select(-contains("DEPRECATED"))
    # fix the names by dropping "NEW"
    names(data) = gsub("\\.NEW([0-9]+)", ".\\1", names(data))

    # updating the last.fu
    data[, last.fu] = pmin(data[, last.fu], max(reduced))
    data[, last.fu] = match(data[, last.fu], reduced)

    # updating the statuses
    for (i in 1:length(reduced)) {
      data[, sprintf("%s%s%s", status, sep, i)] =
        ifelse(data[, last.fu] > i, 1, ifelse(data[, last.fu] == i, data[, "DEAD19"], NA))
    }

    # Removing those NA's (who failed or dropped out)
    data = data %>%
      dplyr::filter(!is.na(last.fu))

    return(data)
  }

  ### 2. specific.failure: Treat failure from other causes as censoring.
  ##     See https://health.mo.gov/data/documentation/death/death-icd10.php for the list of the codes
  specific.failure = function(data, cause = "^I.*") {
    fail = grepl(cause, data$UCOD)
    cens = data$UCOD == ""
    cens.other.cause = !(fail | cens)
    data[cens.other.cause, "DEAD19"] = 0
    for (i in 1:max(data$last.fu, na.rm = TRUE)) {
      matched = data[cens.other.cause, "last.fu"] == i
      after.matched = data[cens.other.cause, "last.fu"] < i
      data[cens.other.cause, paste0("d.", i)][matched] = 0         # Change death = 1 to = 0
      data[cens.other.cause, paste0("d.", i)][after.matched] = NA  # Change delta = NA afterwards
    }
    return(data)
  }

  ### dat.567 has only seven observed deaths. So we move the window to 4:6.
  ### dat.456 still has only seven observed deaths. So we move the window to 5:6.
  # dat.567 = dat.clean %>% trimmer(reduced = 5:7) %>% specific.failure %>% filter(PRVCHD05 == 1| PREVHF01 == 1) # 92% censoring
  # dat.567 = dat.clean %>% trimmer(reduced = 5:7) %>% specific.failure %>% dplyr::filter(hf.1 == 1) # 81% censoring
  # dat.567 = dat.clean %>% trimmer(reduced = 5:7) %>% specific.failure %>% dplyr::filter(hf.1 == 1 | PRVCHD05 == 1| PREVHF01 == 1) # 82% censoring
  # dat.456 = dat.clean %>% trimmer(reduced = 4:6) %>% specific.failure %>% dplyr::filter( PRVCHD05 == 1| PREVHF01 == 1) # 76% censoring
  dat.56 = dat.clean %>% trimmer(reduced = 5:6) %>% {if (all_cause) . else specific.failure(.)} %>% dplyr::filter(hf.1 == 1 | PRVCHD05 == 1| PREVHF01 == 1) # 82% censoring
  # dat.567 %>% filter(DEAD19 == 1) %>% {table(.$last.fu)} # 147, 20, 7 (2011, 2017, 2019) # Too few deaths at the end
  # dat.456 %>% filter(DEAD19 == 1) %>% {table(.$last.fu)} # 187, 37, 8 (1996, 2011, 2017) # Too few deaths at the end
  dat.56 %>% filter(DEAD19 == 1) %>% {table(.$last.fu)} # 147, 27 (2011, 2017) / 324, 65 for all cause deaths.
  # dat.45 %>% filter(DEAD19 == 1) %>% {table(.$last.fu)} # Too old (1996, 2011)

  # combineTx = function(vec1, vec2) {ifelse(is.na(vec1) | is.na(vec2), NA, paste0(vec1, vec2))}
  combineTx = function(vec1, vec2) {vec1 + vec2}
  # dat.456 = dat.456 %>%
  dat.56 = dat.56 %>%
    dplyr::filter(!race %in% c("A", "I")) %>% mutate(race = race %>% droplevels) %>%   # Removing the only one "A" and no "I"s from analysis.
    rename(hf = hf.1) %>%               # Baseline feature, only available at stage 5.  #hf is not available at visit 4.
    rename(delta = DEAD19)              # Replace DEAD19 with delta
  saveRDS(dat.56, sprintf("../data_ARIC/02.aric.56%s.rds", if (all_cause) "_allcause" else "" ))

  # 1. raw data:  n=15,760, % non-censoring = 55%
  dim(dat.clean) ; mean(dat.clean$DEAD19)

  # 2. visits 5-7 only:  n=5,890, % non-censoring = 23%
  dim(trimmer(dat.clean, reduced = 5:6)) ; mean(trimmer(dat.clean, reduced = 5:6)$DEAD19)

  # 3. circular death only:  n=5,890, % non-censoring = 7.6%
  dim(specific.failure(trimmer(dat.clean, reduced = 5:6))) ; mean(specific.failure(trimmer(dat.clean, reduced = 5:6))$DEAD19)

  # 3. the cases at visit 1: n = 945,  % non-censoring = 18% / 41% for all cause mortality
  dim(dat.56) ; mean(dat.56$delta)



  ### STEP 3.
  ### Handling the missing values.
  library(mice)  # multiple imputation
  all_cause = TRUE
  dat.56 = readRDS(sprintf("../data_ARIC/02.aric.56%s.rds", if (all_cause) "_allcause" else "" ))
  m = 5 # Five imputation sets

  ### Composition of the missing values
  ##  type A. Intrinsic missing values -- For those who were censored or failed at a previous stage, the current stage values are inherently missing.
  ##  type B. Undue missing values -- Lack of information due to any other reasons than type A.
  ### Approches
  ##  Only type B missings are imputed.
  ##  For the remaining, multiple imputation beginging from baseline, stage 1, 2, and 3.
  ##  The multiple imputation is done for the first stage (baseline), and for each imputed set, the next stage imputation is done once.

  ### 1. Staging the variables
  names(dat.56)
  var.y = c("ID", "T.0", "delta", "last.fu", "UCOD", "d.1", "V.1", "d.2", "V.2")
  var.0 = c("PREVHF01", "PRVCHD05", "TIAB01", "CENTER", "gender", "race", "age", "cig_yrs")
  var.1 = c("AA.1", "AS.1", "AC.1", "bmi.1", "wth.1", "drink.1", "hypert.1", "glucose.1", "smoke.1", "hdl.1", "fast.1", "lfast.1", "hf")
  var.2 = c("AA.2", "AS.2", "AC.2", "bmi.2", "wth.2", "drink.2", "hypert.2", "glucose.2", "smoke.2", "hdl.2", "fast.2", "lfast.2")
  # var.3 = c("AA.3", "AS.3", "AC.3", "bmi.3", "wth.3", "drink.3", "hypert.3", "glucose.3", "smoke.3", "hdl.3", "fast.3", "lfast.3")
  var.0a = c("ID", var.0); var.1a = c("ID", var.1); var.2a = c("ID", var.2); #var.3a = c("ID", var.3)

  ### 2. Imputation
  ## 2.0 Stage 0 -- outcomes and the baseline variables
  dat.56[!complete.cases(dat.56[, var.y]), var.y] # complete

  # All less than 3% missing. PREVHF01 (1.3%), PRVCHD05 (1.6%), TIAB01 (2.2%), cig_yrs (1.0%)
  #  missing fractions
  apply(dat.56[, var.0], 2, function(x) mean(is.na(x)))

  set.seed(0)
  imp.0 = mice(dat.56[, var.0a], m = m, block = list("ID", var.0))  # Using block=, "ID" is not used to guess the other variables.
  dat.comp.list0 =
    lapply(1:m, function(i) complete(imp.0, action = i))
  saveRDS(dat.comp.list0, sprintf("../data_ARIC/03.tmp0%s.rds", if (all_cause) "_allcause" else "" ))

  ## 2.1 Stage 1
  #  missing fractions
  apply(dat.56[, var.1], 2, function(x) mean(is.na(x)))

  set.seed(1)
  dat.comp.list1 =
    lapply(1:m, function(i) {
      augmented.i = left_join(dat.comp.list0[[i]], dat.56[, var.1a], by = "ID")
      imp.i = mice(augmented.i, m = 1, block = list("ID", c(var.0, var.1))) # For each previously imputed set, imputation is done once.
      comp.i = complete(imp.i, action = 1)
      comp.i
    })
  saveRDS(dat.comp.list1, sprintf("../data_ARIC/03.tmp1%s.rds", if (all_cause) "_allcause" else "" ))

  ## 2.2 Stage 2
  #  missing fractions
  avail2 = !dat.56$d.2 %>% is.na      # Those who are available (no failure or dropout so far)
  apply(dat.56[avail2, var.2], 2, function(x) mean(is.na(x)))

  set.seed(2)
  dat.comp.list2 =
    lapply(1:m, function(i) {
      augmented.i =
        left_join(dat.comp.list1[[i]], dat.56[, var.2a], by = "ID") %>%
        dplyr::filter(avail2)    ### Consider only those who were observed!
      imp.i = mice(augmented.i, m = 1, block = list("ID", c(var.0, var.1, var.2))) # For each previously imputed set, imputation is done once.
      comp.i = complete(imp.i, action = 1)
      left_join(dat.comp.list1[[i]], comp.i[, var.2a], by = "ID")
    })
  saveRDS(dat.comp.list2, sprintf("../data_ARIC/03.tmp2%s.rds", if (all_cause) "_allcause" else "" ))



  dat.comp.list =
    dat.comp.list2 %>%
    lapply(function(i) {
      i %>%
        left_join(dat.56[, var.y], ., by = "ID") %>%
        mutate(#ACFS.1 = combineTx(AC.1 == 1, lfast.1 == 1) %>% factor,
          #ACFS.2 = combineTx(AC.2 == 1, lfast.2 == 1) %>% factor,
          #ACFS8.1 = combineTx(AC.1 == 1, fast.1 == 1) %>% factor,
          #ACFS8.2 = combineTx(AC.2 == 1, fast.2 == 1) %>% factor,
          AAprev.2 = AA.1,
          ACprev.2 = AC.1,
          ASprev.2 = AS.1,
          fastprev.2 = fast.1,
          lfastprev.2 = lfast.1,
          fast3prev.2 = ifelse(is.na(lfast.1), NA, ifelse(lfast.1, 2, fast.1)) %>% factor,
          ACFS.1 = paste0(AC.1, lfast.1) %>% {ifelse(. == "NANA", NA, .)} %>% factor,
          ACFS.2 = paste0(AC.2, lfast.2) %>% {ifelse(. == "NANA", NA, .)} %>% factor)
    })
  saveRDS(dat.comp.list, sprintf("../data_ARIC/dat56/03.aric.comp%s.rds",  if (all_cause) "_allcause" else "" ))
  saveRDS(dat.comp.list, sprintf("../data_ARIC/03.aric.comp%s.rds",  if (all_cause) "_allcause" else "" ))

