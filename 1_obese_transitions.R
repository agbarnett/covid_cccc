# 1_obese_transitions.R
# called by 1_obese.Rmd
# create transitions for "full" model with 4 states using imputed data
# May 2021

### work out patients to exclude because of missing data in their dates ###

### Section 1: ECMO patient data ###

### start by generating transitions with no intermediate state ###
trans.data = time_transitions_no_intermediate(indata = patients, 
                                              start_state = 'date_ecmo', # starting date/state
                                              discharge_date = 'date_hospital_discharge', # ICU (date_discharge) or hospital discharge (date_hospital_discharge)
                                              check = TRUE, 
                                              censor_day = censor_day) %>%
  rename('entry' = 'start', # use these names to suit mstate function below
         'exit' = 'end')
# check transfers, look fine, no need to apply transfers
# filter(trans.data, !is.na(date_transfer)) %>% select(pin, date_ecmo, date_discharge, date_transfer, last_date)
# now remove excluded patients with too much missing date data
excluded = patients$pin[patients$pin %in% trans.data$pin ==FALSE]
patients = filter(patients, !pin %in% excluded)
# check exclusions, all have something badly wrong with date
#f = filter(patients, pin==excluded[1]) %>% select(pin, date_ecmo, date_discharge, date_death, date_hospital_discharge, last_date, date_transfer)

# now repeat for imputed data
trans.data.imputed = list()
for (k in 1:N.impute){
  trans.data.imputed[[k]] = time_transitions_no_intermediate(indata = patients_imputed[[k]], 
                                                start_state = 'date_ecmo', # starting date/state
                                                discharge_date = 'date_hospital_discharge', # ICU (date_discharge) or hospital discharge (date_hospital_discharge)
                                                check = TRUE, 
                                                censor_day = censor_day) %>%
    rename('entry' = 'start', # use these names to suit mstate function below
           'exit' = 'end')
  trans.data.imputed[[k]] = filter(trans.data.imputed[[k]], !pin %in% excluded)
}

### Section 2: repeat above for MV patients (can be used for time to ECMO model and full model) ###
trans.data.mv = time_transitions_string(
  indata = mv_patients,
  dates = c('date_mechanical_ventilation','date_ecmo','date_transfer','date_hospital_discharge','date_death'),
  int_state = 'date_ecmo', # intermediate
  start_state = 'date_mechanical_ventilation', # starting date/state
  discharge_date = 'date_hospital_discharge', # hospital discharge 
  check = TRUE, 
  censor_day = censor_day)
# transfers are censored
trans.data.mv = mutate(trans.data.mv,
                       to = ifelse(to=='date_transfer','censored', to)) 
# patients with too little data:
mv_excluded = mv_patients$pin[mv_patients$pin %in% trans.data.mv$pin ==FALSE]
# now remove from patients
mv_patients = filter(mv_patients, !pin %in% mv_excluded)

# now repeat for imputed data
trans.data.mv.imputed = list()
for (k in 1:N.impute){
  trans.data.mv.imputed[[k]] = time_transitions_string(
    indata = mv_patients_imputed[[k]],
    dates = c('date_mechanical_ventilation','date_ecmo','date_transfer','date_hospital_discharge','date_death'),
    int_state = 'date_ecmo', # intermediate
    start_state = 'date_mechanical_ventilation', # starting date/state
    discharge_date = 'date_hospital_discharge', # hospital discharge 
    check = TRUE, 
    censor_day = censor_day) 
  # transfers are censored
  trans.data.mv.imputed[[k]] = mutate(trans.data.mv.imputed[[k]],
                         to = ifelse(to=='date_transfer','censored', to)) 
  # remove patients with too little data
  trans.data.mv.imputed[[k]] = filter(trans.data.mv.imputed[[k]], !pin %in% mv_excluded)
}

## Section 3: additional stuff for mstate ##

# Set transition matrix for mstate
# three states (ECMO only)
tra3 <- transMat(x = list(c(2,3), c(), c()), 
                names = c("ECMO", "Discharge","Death"))
# four states
tra <- transMat(x = list(c(2,3,4), c(3,4), c(),c()), # no ECMO ->MV, add 1 to second set if this is needed
                names = c("MV", "ECMO", "Discharge","Death"))

# change transition names to numbers
trans.data.mv = mutate(trans.data.mv,
                         to = case_when(
                           to == 'censored' ~ 0,
                           to == 'date_int_state_start' ~ 2,
                           to == 'discharge' ~ 3,
                           to == 'death' ~ 4),
                         from = case_when(
                           from =='date_mechanical_ventilation' ~ 1,
                           from =='date_int_state_start' ~ 2
                         ))

# Add transition vectors
trans.data.mv$trans <- NA
for (i in 1:nrow(tra)){
  for (j in 1:ncol(tra))  {
    trans.data.mv$trans[which(trans.data.mv$from == i & trans.data.mv$to ==j)] <- tra[i, j]  
  }
}

# Add status vector, indicates observed transition
trans.data.mv$status <- 1
# Status vector for censored observations set to '0'
trans.data.mv$status[trans.data.mv$to == 0] <- 0
# Rename 'to' == 0  to 'cens' for use in 'ext_mstate' function
trans.data.mv$to[trans.data.mv$to == 0] <- 'cens'
trans.data.mv$from = as.character(trans.data.mv$from) # needed for function below
#my.data$id = my.data$pin #needed for function below
# Create data frame with all possible transitions for when a patient is at risk
trans.data.mv_ext = ext_mstate(trans.data.mv, tra, id.var='pin')



