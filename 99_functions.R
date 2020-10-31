# 99_functions.R
# range of functions used by multiple programs 
# October 2020

# function for rounding numbers with zeros kept
roundz = function(x, digits){
  dformat = paste('%.', digits, 'f', sep='')
  x = sprintf(dformat, round(x, digits))
  return(x)
}

# function for making nice tables, used by table1.markdown.Rmd
make_nice = function(intab, rows, pvalue=NULL){
  t1 = intab$cross_table[, c(3,1,2)] # re-order columns
  p1 = round(100*intab$proportions[, c(3,1,2)]) # re-order columns and round to percents
  t1 = t1[rownames(t1) %in% rows, ] # get rows
  p1 = p1[rownames(p1) %in% rows, ]
  tab = NULL
  if(length(rows)==1){
    c1 = c(rows, paste(t1, ' (', p1, ')', sep='')) # combine n and percent
    tab = rbind(tab, c1)
  }
  if(length(rows)>1){
    names = rownames(t1) # to get the correct ordering
    for (k in 1:length(rows)){
      c1 = c(names[k], paste(t1[k,], ' (', p1[k,], ')', sep='')) # combine n and percent
      tab = rbind(tab, c1)
    }
  }
  # add pvalue (just to top row)
  if(is.null(pvalue)==FALSE){
    tab = cbind(tab, rep('', nrow(tab)))
    tab[1,5] =  pvalue
  }
  #
  return(tab)
}

# standard error of the mean
sem <- function(x, na.rm=TRUE){
  x = x[is.na(x)==FALSE] # remove missing
  res = sqrt(var(x)/length(x))
  return(res)
}

# for colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# function to compare blood gas data before and after the ECMO date
blood_gases = function(before_var='ecmo_worst_pa_o2_6hr_before',
                       after_var = 'pa_o2'){
  # cheat to make variables
  patients$before = select(patients, before_var)[,1]
  daily$after = select(daily, after_var)[,1]
  # 
  to_compare = filter(patients, 
                      !is.na(before), # must have reading and ECMO data
                      !is.na(date_ecmo)) %>%
    select(pin, date_ecmo, before)
  to_merge_from_daily = select(daily, date_daily, pin, after) %>%
    filter(!is.na(date_daily),
           !is.na(after)) # just non missing
  for_table = left_join(to_compare, to_merge_from_daily, by='pin') %>%
    mutate(diff = as.numeric(date_ecmo - date_daily)) %>%
    arrange(pin, diff) %>%
    filter(diff >= 1) %>% # must be at least one day after ECMO date
    group_by(pin) %>%
    slice(1) %>% # select daily date closest to ECMO date
    ungroup()
  # data for plot
  to_plot = gather(for_table, `before`, `after`, key='time', value='result') %>%
    mutate(xaxis = ifelse(str_detect(string=time, 'before'), 1, 2),
           xaxis = factor(xaxis, levels=1:2, labels=c('Worst before ECMO','After ECMO')))
  # data for average change line
  av_line = group_by(to_plot, xaxis) %>%
    summarise(result = mean(result)) %>%
    mutate(pin = -99) %>%
    ungroup()
  
  # test difference using paired t-test (after - before)
  one = log(for_table$before) # log transform
  two = log(for_table$after)
  test = t.test(y=one, x=two, paired=TRUE, alternative='two.sided')
  # percent difference
  test$pdiff = 100*(exp(test$estimate) - 1)
  test$conf.int = 100*(exp(test$conf.int) - 1)
  
  # return
  to_return = list()
  to_return$av_line = av_line
  to_return$to_plot = to_plot
  to_return$test = test
  return(to_return)
} # End of function


### function to make transitions times for mechanical ventilation as intermediate state
time_transitions = function(indata, check=TRUE, censor_day = NA){
  
  # State 0: Censored
  # State 1: ICU, non-MV
  # State 2: MV (intermediate state)
  # State 3: Discharge
  # State 4: Death
  
  # exclude those with no event dates (may have ICU entry date)
  index = is.na(indata$date_mechanical_ventilation) & 
    is.na(indata$date_death) &
  is.na(indata$date_discharge) &
  is.na(indata$last_date)&
    is.na(indata$date_mech_vent_discontinued)
  n_excluded = sum(index)
  indata = indata[!index,]

  # if date mechanical ventilation stopped is same as death or discharge then delete MV-discontinued date
  indata = mutate(indata,
      date_mech_vent_discontinued = ifelse(date_mech_vent_discontinued==date_discharge, NA, date_mech_vent_discontinued),
      date_mech_vent_discontinued = ifelse(date_mech_vent_discontinued==date_death, NA, date_mech_vent_discontinued),
      date_mech_vent_discontinued = as.Date(date_mech_vent_discontinued, origin='1970-01-01'))
  
# a) times to death
to_death = filter(indata, is.na(date_death)==FALSE) %>% # must have date of death
  mutate(start_date = ifelse(is.na(date_mechanical_ventilation)==FALSE, date_mechanical_ventilation, date_icu), # start from MV or admission ...
         start_date = ifelse(is.na(date_mech_vent_discontinued)==FALSE, date_mech_vent_discontinued, start_date), # ... update to MV end if it's there
         start_date = as.Date(start_date, origin='1970-01-01'), 
         from = ifelse(is.na(date_mechanical_ventilation)==FALSE, 2, 1), # MV, non-MV
         from = ifelse(is.na(date_mech_vent_discontinued)==FALSE, 1, from), # from non-MV if stop-MV date
         to = 4, # death
         start = as.numeric(start_date - date_icu), # start time in days
         end = as.numeric(date_death - date_icu)) # end time in days (to death)

# b) times to discharge
to_discharge = filter(indata, is.na(date_death)==TRUE, is.na(date_discharge)==FALSE) %>% # must have no death date but discharge date
  mutate(start_date = ifelse(is.na(date_mechanical_ventilation)==FALSE, date_mechanical_ventilation, date_icu), # start from MV or admission
         start_date = ifelse(is.na(date_mech_vent_discontinued)==FALSE, date_mech_vent_discontinued, start_date), # ... update to MV end if it's there
         start_date = as.Date(start_date, origin='1970-01-01'), 
         from = ifelse(is.na(date_mechanical_ventilation)==FALSE, 2, 1), # MV, non-MV
         from = ifelse(is.na(date_mech_vent_discontinued)==FALSE, 1, from), # from non-MV if stop-MV date
         to = 3, # Discharge
         start = as.numeric(start_date - date_icu),
         end = as.numeric(date_discharge - date_icu)) # end time in days (to discharge)

# c) times to censored
to_censored = filter(indata, is.na(date_death)==TRUE, is.na(date_discharge)==TRUE) %>%
  mutate(start_date = ifelse(is.na(date_mechanical_ventilation)==FALSE, date_mechanical_ventilation, date_icu),
         start_date = as.Date(start_date, origin='1970-01-01'), 
         from = ifelse(is.na(date_mechanical_ventilation)==FALSE, 2, 1), # MV, non-MV
         to = 0, # 'Censored',
         start = as.numeric(start_date - date_icu),
         end = as.numeric(last_date - date_icu)) # days to censored

# d) times to MV (from ICU)
to_MV = filter(indata, is.na(date_mechanical_ventilation)==FALSE) %>%
  mutate(start_date = date_icu,
         start_date = as.Date(start_date, origin='1970-01-01'), 
         from = 1, # Non-MV
         to = 2, # MV
         start = 0, 
         end = as.numeric(date_mechanical_ventilation - date_icu), # days to MV
         end = ifelse(end < 0, 0, end)) # if negative assume ventilated on admission

# e) times to non-MV (from MV)
to_non_MV = filter(indata, is.na(date_mech_vent_discontinued)==FALSE) %>%
  mutate(start_date = date_mechanical_ventilation,
         end_date = date_mech_vent_discontinued,
         from = 2, # MV
         to = 1, # non-MV
         start = as.numeric(start_date - date_icu), # days to MV start
         end = as.numeric(end_date - date_icu), # days to MV end
         start = ifelse(start<0, 0, start)) # assume ventilated on arrival if negative

# combine five transitions 
for_model = bind_rows(to_death, to_discharge, to_censored, to_MV, to_non_MV) %>%
  mutate(start = ifelse(start < 0, 0, start)) # assume negative dates mean ventilated on arrival

## fixes because of same dates
# if non-MV to MV at time zero then assume ventilated on entry (so no transition)
index = with(for_model, from==1 & to==2 & start==0 & end==0)
for_model = for_model[!index,]
# if MV to non-MV at time zero then assume ventilated on entry (so no transition)
index = with(for_model, from==2 & to==1 & start==0 & end==0)
for_model = for_model[!index,]

## add half days
# add half to discharge if start and end are same
index = with(for_model, from==1 & to==2 & start==end)
for_model$end[index] = for_model$end[index] + 0.5
# add half to discharge if start and end are same
index = with(for_model, from==2 & to==1 & start==end)
for_model$end[index] = for_model$end[index] + 0.5
# add half to discharge if start and end are same
index = with(for_model, to==3 & start==end)
for_model$end[index] = for_model$end[index] + 0.5
# add half to death if start and end are same
index = with(for_model, to==4 & start==end)
for_model$end[index] = for_model$end[index] + 0.5
# add half to censored if start and end are same
index = with(for_model, to==0 & start==end)
for_model$end[index] = for_model$end[index] + 0.5

# remove impossible patients with obvious errors in dates
exclude = filter(for_model, start > end) 
for_model = filter(for_model,
                   !pin %in% exclude$pin) %>%
  arrange(pin, start, end)

#remove any remaining NAs
index = with(for_model,is.na(from)|is.na(to)|is.na(start)|is.na(end))
for_model = for_model[!index,]

# censor at a specified time (if it exists)
if(is.na(censor_day) >= FALSE){
  for_model = filter(for_model, start < censor_day) %>% # only transitions before censor day
    mutate(beyond = end > censor_day,
           to = ifelse(beyond==TRUE, 0, to), # change state to censored
           end = ifelse(beyond==TRUE, censor_day, end))
}

# quick check
if(check==TRUE){
  t1 = filter(for_model, from==1, to==2) %>%
    sample_n(1)
  t2 = filter(for_model, from==1, to==3) %>%
    sample_n(1)
  t3 = filter(for_model, from==1, to==4) %>%
    sample_n(1)
  t4 = filter(for_model, from==2, to==3) %>%
    sample_n(1)
  t5 = filter(for_model, from==2, to==4) %>%
    sample_n(1)
  t6 = filter(for_model, to==0) %>%
    sample_n(1)
  t7 = filter(for_model, to==1) %>%
    sample_n(1)
  test = bind_rows(t1, t2, t3, t4, t5, t6, t7)
  ids = unique(test$pin)
  cat('Original data:\n')
  patient_test = filter(indata, pin %in% ids) %>%
    select(pin, date_icu, starts_with('date_mech'), date_discharge, date_death, last_date) %>%
    arrange(pin)
  print(patient_test)
  cat('Transitions:\n')
  transition_test = filter(for_model, pin %in% ids) %>%
    select(pin, from, to, start, end)
  print(transition_test)
  
}

return(for_model)
} # end of function

#
### function to make transitions times for ECMO as intermediate state
time_transitions_ecmo = function(indata, check=TRUE, censor_day = NA){
  
  # State 0: Censored
  # State 1: ICU, non-ECMO
  # State 2: ECMO (intermediate state)
  # State 3: Discharge
  # State 4: Death
  
  # exclude those with no event dates (may have ICU entry date)
  index = is.na(indata$date_ecmo) & 
    is.na(indata$date_death) &
    is.na(indata$date_discharge) &
    is.na(indata$last_date)&
    is.na(indata$date_ecmo_discontinued)
  n_excluded = sum(index)
  indata = indata[!index,]
  
  # if date ECMO stopped is same as death or discharge then delete ECMO discontinued date
  indata = mutate(indata,
                  date_ecmo_discontinued = ifelse(date_ecmo_discontinued==date_discharge, NA, date_ecmo_discontinued),
                  date_ecmo_discontinued = ifelse(date_ecmo_discontinued==date_death, NA, date_ecmo_discontinued),
                  date_ecmo_discontinued = as.Date(date_ecmo_discontinued, origin='1970-01-01'))
  
  # a) times to death
  to_death = filter(indata, is.na(date_death)==FALSE) %>% # must have date of death
    mutate(start_date = ifelse(is.na(date_ecmo)==FALSE, date_ecmo, date_icu), # start from ECMO or admission ...
           start_date = ifelse(is.na(date_ecmo_discontinued)==FALSE, date_ecmo_discontinued, start_date), # ... update to ECMO end if it's there
           start_date = as.Date(start_date, origin='1970-01-01'), 
           from = ifelse(is.na(date_ecmo)==FALSE, 2, 1), # ECMO, non-ECMO
           from = ifelse(is.na(date_ecmo_discontinued)==FALSE, 1, from), # from non-ECMO if stop-ECMO date
           to = 4, # death
           start = as.numeric(start_date - date_icu), # start time in days
           end = as.numeric(date_death - date_icu)) # end time in days (to death)
  
  # b) times to discharge
  to_discharge = filter(indata, is.na(date_death)==TRUE, is.na(date_discharge)==FALSE) %>% # must have no death date but discharge date
    mutate(start_date = ifelse(is.na(date_ecmo)==FALSE, date_ecmo, date_icu), # start from ECMO or admission
           start_date = ifelse(is.na(date_ecmo_discontinued)==FALSE, date_ecmo_discontinued, start_date), # ... update to ECMO end if it's there
           start_date = as.Date(start_date, origin='1970-01-01'), 
           from = ifelse(is.na(date_ecmo)==FALSE, 2, 1), # ECMO, non-ECMO
           from = ifelse(is.na(date_ecmo_discontinued)==FALSE, 1, from), # from non-ECMO if stop-ECMO date
           to = 3, # Discharge
           start = as.numeric(start_date - date_icu),
           end = as.numeric(date_discharge - date_icu)) # end time in days (to discharge)
  
  # c) times to censored
  to_censored = filter(indata, is.na(date_death)==TRUE, is.na(date_discharge)==TRUE) %>%
    mutate(start_date = ifelse(is.na(date_ecmo)==FALSE, date_ecmo, date_icu),
           start_date = as.Date(start_date, origin='1970-01-01'), 
           from = ifelse(is.na(date_ecmo)==FALSE, 2, 1), # ECMO, non-ECMO
           to = 0, # 'Censored',
           start = as.numeric(start_date - date_icu),
           end = as.numeric(last_date - date_icu)) # days to censored
  
  # d) times to ECMO (from ICU)
  to_ECMO = filter(indata, is.na(date_ecmo)==FALSE) %>%
    mutate(start_date = date_icu,
           start_date = as.Date(start_date, origin='1970-01-01'), 
           from = 1, # Non-ECMO
           to = 2, # ECMO
           start = 0, 
           end = as.numeric(date_ecmo - date_icu), # days to ECMO
           end = ifelse(end < 0, 0, end)) # if negative assume ventilated on admission
  
  # e) times to non-ECMO (from ECMO)
  to_non_ECMO = filter(indata, 
                       is.na(date_ecmo)==FALSE, # must have both dates
                       is.na(date_ecmo_discontinued)==FALSE) %>%
    mutate(start_date = date_ecmo,
           end_date = date_ecmo_discontinued,
           from = 2, # ECMO
           to = 1, # non-ECMO
           start = as.numeric(start_date - date_icu), # days to ECMO start
           end = as.numeric(end_date - date_icu), # days to ECMO end
           start = ifelse(start<0, 0, start)) # assume ventilated on arrival if negative
  
  # combine five transitions 
  for_model = bind_rows(to_death, to_discharge, to_censored, to_ECMO, to_non_ECMO) %>%
    mutate(start = ifelse(start < 0, 0, start)) # assume negative dates mean ventilated on arrival
  
  ## fixes because of same dates
  # if non-ECMO to ECMO at time zero then assume ventilated on entry (so no transition)
  index = with(for_model, from==1 & to==2 & start==0 & end==0)
  for_model = for_model[!index,]
  # if ECMO to non-ECMO at time zero then assume ventilated on entry (so no transition)
  index = with(for_model, from==2 & to==1 & start==0 & end==0)
  for_model = for_model[!index,]
  ## add half days
  # add half to discharge if start and end are same
  index = with(for_model, from==1 & to==2 & start==end)
  for_model$end[index] = for_model$end[index] + 0.5
  # add half to discharge if start and end are same
  index = with(for_model, from==2 & to==1 & start==end)
  for_model$end[index] = for_model$end[index] + 0.5
  # add half to discharge if start and end are same
  index = with(for_model, to==3 & start==end)
  for_model$end[index] = for_model$end[index] + 0.5
  # add half to death if start and end are same
  index = with(for_model, to==4 & start==end)
  for_model$end[index] = for_model$end[index] + 0.5
  # add half to censored if start and end are same
  index = with(for_model, to==0 & start==end)
  for_model$end[index] = for_model$end[index] + 0.5
  
  # remove impossible patients with obvious errors in dates
  exclude = filter(for_model, start > end) 
  for_model = filter(for_model,
                     !pin %in% exclude$pin) %>%
    arrange(pin, start, end)
  
  # censor at a specified time (if it exists)
  if(is.na(censor_day) >= FALSE){
    for_model = filter(for_model, start < censor_day) %>% # only transitions before censor day
      mutate(beyond = end > censor_day,
             to = ifelse(beyond==TRUE, 0, to), # change state to censored
             end = ifelse(beyond==TRUE, censor_day, end))
  }
  
  # check one patient
  # check = filter(patients, pin=='OX_00634-0028') %>% select(pin, date_icu, date_ecmo, date_ecmo_discontinued, date_discharge, date_death)
  # check missing start dates
  # check = filter(patients, is.na(date_ecmo), !is.na(date_ecmo_discontinued)) %>% select(pin, date_icu, date_ecmo, date_ecmo_discontinued, date_discharge, date_death)
  
  # quick check
  if(check==TRUE){
    t1 = filter(for_model, from==1, to==2) %>%
      sample_n(1)
    t2 = filter(for_model, from==1, to==3) %>%
      sample_n(1)
    t3 = filter(for_model, from==1, to==4) %>%
      sample_n(1)
    t4 = filter(for_model, from==2, to==3) %>%
      sample_n(1)
    t5 = filter(for_model, from==2, to==4) %>%
      sample_n(1)
    t6 = filter(for_model, to==0) %>%
      sample_n(1)
    t7 = filter(for_model, to==1) %>%
      sample_n(1)
    test = bind_rows(t1, t2, t3, t4, t5, t6, t7)
    ids = unique(test$pin)
    cat('Original data:\n')
    patient_test = filter(indata, pin %in% ids) %>%
      select(pin, date_icu, starts_with('date_mech'), date_discharge, date_death, last_date) %>%
      arrange(pin)
    print(patient_test)
    cat('Transitions:\n')
    transition_test = filter(for_model, pin %in% ids) %>%
      select(pin, from, to, start, end)
    print(transition_test)
    
  }
  
  return(for_model)
} # end of function

### function to make transitions times for ECMO with prone as an intermediate state
# addition option with prone_time_dependent = TRUE to add prone as a cumulative exposire
time_transitions_ecmo_prone = function(indata, 
                                       check = TRUE, 
                                       censor_day = NA, 
                                       prone_time_dependent = FALSE){
  
  # State 0: Censored
  # State 1: ECMO & supine
  # State 2: ECMO & prone (intermediate state)
  # State 3: Discharge
  # State 4: Death
  
  # exclude those with absolutely no event dates, must be on ECMO
  index = 
    is.na(indata$date_ecmo) & 
    is.na(indata$date_death) &
    is.na(indata$date_discharge) &
    is.na(indata$last_date) & 
    is.na(indata$prone_start)
  n_excluded = sum(index)
  indata = indata[!index,]
  
  # date fixes
  indata = mutate(indata,
          # if prone end before ECMO then blank prone start and end (prone happened before ECMO started)
          prone_blank = ifelse(prone_end < date_ecmo, 1, 0),
          prone_start = ifelse(prone_blank == 1, NA, prone_start),
          prone_end = ifelse(prone_blank == 1, NA, prone_end),
          # if prone start before ECMO (and prone end after) then change to ECMO date (i.e., assume prone on arrival)
          prone_start = ifelse(prone_start < date_ecmo & !is.na(prone_start) & (prone_end >= date_ecmo | is.na(prone_end)), date_ecmo, prone_start),
          # if last date is before ecmo start then move last date
          last_date  = ifelse(last_date < date_ecmo, date_ecmo, last_date),
          # if date prone end is same as death or discharge then delete prone end date
          prone_end = ifelse(prone_end==date_discharge & !is.na(prone_end), NA, prone_end),
          prone_end = ifelse(prone_end==date_death & !is.na(prone_end), NA, prone_end),
          # if date ECMO end is same as death or discharge then delete ECMO end date
          ecmo_end = ifelse(ecmo_end==date_discharge & !is.na(ecmo_end), NA, ecmo_end),
          ecmo_end = ifelse(ecmo_end==date_death& !is.na(ecmo_end), NA, ecmo_end),
          # reformat dates
          prone_start = as.Date(prone_start, origin='1970-01-01'),
          prone_end = as.Date(prone_end, origin='1970-01-01'),
          last_date  = as.Date(last_date, origin='1970-01-01'),
          ecmo_end = as.Date(ecmo_end, origin='1970-01-01')
          )
  
  ### versions without time-dependent prone ###
  if(prone_time_dependent == FALSE){
                    
  # a) times to death
  to_death = filter(indata, is.na(date_death)==FALSE) %>% # must have date of death
    mutate(start_date = ifelse(is.na(prone_start)==FALSE, prone_start, date_ecmo), # start from Prone or supine ...
           start_date = ifelse(is.na(prone_end)==FALSE, prone_end, start_date), # ... update to Prone end if it's there
           start_date = as.Date(start_date, origin='1970-01-01'), 
           from = ifelse(is.na(prone_start)==FALSE, 2, 1), # if non-missing start date then Prone, ECMO
           from = ifelse(is.na(prone_end)==FALSE, 1, from), # from ECMO & supine if non-missing stop-Prone date
           to = 4, # death
           start = as.numeric(start_date - date_ecmo), # start time in days
           end = as.numeric(date_death - date_ecmo)) # end time in days (to death)
  
  # b) times to discharge
  to_discharge = filter(indata, is.na(date_death)==TRUE, is.na(date_discharge)==FALSE) %>% # must have no death date but discharge date
    mutate(start_date = ifelse(is.na(prone_start)==FALSE, prone_start, date_ecmo), # start from Prone or supine
           start_date = ifelse(is.na(prone_end)==FALSE, prone_end, start_date), # ... update to Prone end if it's there
           start_date = as.Date(start_date, origin='1970-01-01'), 
           from = ifelse(is.na(prone_start)==FALSE, 2, 1), # Prone or supine
           from = ifelse(is.na(prone_end)==FALSE, 1, from), # from supine if stop-Prone date
           to = 3, # Discharge
           start = as.numeric(start_date - date_ecmo),
           end = as.numeric(date_discharge - date_ecmo)) # end time in days (to discharge)
  
  # c) times to censored
  to_censored = filter(indata, is.na(date_death)==TRUE, is.na(date_discharge)==TRUE) %>%
    mutate(start_date = ifelse(is.na(prone_start)==FALSE, prone_start, date_ecmo), # start from Prone or supine
           start_date = ifelse(is.na(prone_end)==FALSE, prone_end, start_date), # ... update to Prone end if it's there
           start_date = as.Date(start_date, origin='1970-01-01'), 
           from = ifelse(is.na(prone_start)==FALSE, 2, 1), # Prone, non-Prone
           from = ifelse(is.na(prone_end)==FALSE, 1, from), # from supine if stop-Prone date
           to = 0, # 'Censored',
           start = as.numeric(start_date - date_ecmo),
           end = as.numeric(last_date - date_ecmo)) # days to censored
  
  # d) times to Prone (from ECMO)
  to_Prone = filter(indata, is.na(prone_start)==FALSE) %>%
    mutate(start_date = date_ecmo,
           start_date = as.Date(start_date, origin='1970-01-01'), 
           from = 1, # supine
           to = 2, # Prone
           start = 0, 
           end = as.numeric(prone_start - date_ecmo), # days to Prone
           end = ifelse(end < 0, 0, end)) # if negative assume prone on admission
  
  # e) times to non-Prone (from Prone)
to_non_Prone = filter(indata, is.na(prone_end)==FALSE) %>%
mutate(start_date = prone_start,
end_date = prone_end,
from = 2, # Prone
to = 1, # non-Prone (supine)
start = as.numeric(start_date - date_ecmo), # days to Prone start
end = as.numeric(end_date - date_ecmo), # days to Prone end
start = ifelse(start<0, 0, start)) # assume ventilated on arrival if negative

  
  # combine five transitions 
  for_model = bind_rows(to_death, to_discharge, to_censored, to_Prone, to_non_Prone) %>%
    mutate(start = ifelse(start < 0, 0, start)) # assume negative dates mean ECMO on arrival - to check
  
  
  } # end of time-dependent if
  
  ### version with prone ###
  if(prone_time_dependent == TRUE){
    
    # extra date fix for time-dependent prone
    indata = mutate(indata,
                    # if prone start and end the same add one day otherwise daily output does not work
                    prone_end = ifelse(prone_start==prone_end & !is.na(prone_start) & !is.na(prone_end), prone_end+1, prone_end),
                    prone_end = as.Date(prone_end, origin='1970-01-01'))
    
    # a1) times to death in non-prone
    to_death_non_prone = filter(indata, is.na(date_death)==FALSE, is.na(prone_start)==TRUE) %>% # must have date of death
      mutate(start_date = date_ecmo, # start from ecmo ...
             start_date = as.Date(start_date, origin='1970-01-01'), 
             prone_cum = 0, # cumulative days of prone, always none for this group
             from = 1, # ECMO and supine
             to = 4, # death
             start = as.numeric(start_date - date_ecmo), # start time in days
             end = as.numeric(date_death - date_ecmo)) # end time in days (to death)

    # a2) times to death in prone
    prone_and_dead = filter(indata, is.na(date_death)==FALSE, is.na(prone_start)==FALSE) # must have date of death
    # have to loop
    to_death_prone = NULL
    for (p in 1:nrow(prone_and_dead)){
      this_patient = prone_and_dead[p,]
      for(daily in this_patient$date_ecmo:this_patient$date_death){
        out = mutate(this_patient, 
                     start = as.numeric(as.Date(daily, origin='1970-01-01') - date_ecmo), 
                     end = start + 1, # just moving in days
                     prone_cum_start = (as.Date(daily, origin='1970-01-01') - (prone_start-1)) * as.numeric(as.Date(daily, origin='1970-01-01') > (prone_start-1)), # cumulative days of prone
                     prone_cum_end = (as.Date(daily, origin='1970-01-01') - (prone_end)) * as.numeric(as.Date(daily, origin='1970-01-01') >= (prone_end)), # cumulative days of prone
                     prone_cum_end = ifelse(is.na(prone_end)==TRUE, 0, prone_cum_end), # replace missing with none
                     prone_cum = prone_cum_start - prone_cum_end, # residual effect of prone ...
                     prone_cum = ifelse(as.Date(daily, origin='1970-01-01') >= prone_end & !is.na(prone_end), 0, prone_cum), # ... no residual effect of prone 
                     from = ifelse(daily >= prone_start, 2, 1), # prone if beyond prone date
                     from = ifelse(daily >= prone_end & !is.na(prone_end), 1, from), # supine if beyond prone end date
                     to = ifelse(daily == (prone_start-1) , 2, 0), # to prone if delayed start
                     to = ifelse(daily == (prone_end-1) & !is.na(prone_end), 1, to), # to supine if prone end  
                     to = ifelse(daily == date_death, 4, to)) # death if last day, otherwise censored
        to_death_prone = bind_rows(to_death_prone, out)
      }
    }
    
    # b1) times to discharge in non-prone
    to_discharge_non_prone = filter(indata, is.na(date_death)==TRUE, is.na(date_discharge)==FALSE, is.na(prone_start)==TRUE) %>% # must have date of discharge
      mutate(start_date = date_ecmo, # start from ecmo ...
             start_date = as.Date(start_date, origin='1970-01-01'), 
             prone_cum = 0, # cumulative days of prone, always none for this group
             from = 1, # ECMO & supine 
             to = 3, # discharge
             start = as.numeric(start_date - date_ecmo), # start time in days
             end = as.numeric(date_discharge - date_ecmo)) # end time in days (to death)
    
    # b2) times to discharge in prone
    prone_and_discharge = filter(indata, is.na(date_death)==TRUE, is.na(date_discharge)==FALSE, is.na(prone_start)==FALSE) # must have date of discharge
    # have to loop
    to_discharge_prone = NULL
    for (p in 1:nrow(prone_and_discharge)){
      this_patient = prone_and_discharge[p,]
      for(daily in this_patient$date_ecmo:this_patient$date_discharge){
        out = mutate(this_patient, 
                     start = as.numeric(as.Date(daily, origin='1970-01-01') - date_ecmo), 
                     end = start + 1, # just moving in days
                     prone_cum_start = (as.Date(daily, origin='1970-01-01') - (prone_start-1)) * as.numeric(as.Date(daily, origin='1970-01-01') > (prone_start-1)), # cumulative days of prone
                     prone_cum_end = (as.Date(daily, origin='1970-01-01') - (prone_end)) * as.numeric(as.Date(daily, origin='1970-01-01') >= (prone_end)), # cumulative days of prone
                     prone_cum_end = ifelse(is.na(prone_end)==TRUE, 0, prone_cum_end), # replace missing with none
                     prone_cum = prone_cum_start - prone_cum_end, # residual effect of prone ...
                     prone_cum = ifelse(as.Date(daily, origin='1970-01-01') >= prone_end & !is.na(prone_end), 0, prone_cum), # ... no residual effect of prone 
                     from = ifelse(daily >= prone_start, 2, 1), # prone if beyond prone date
                     from = ifelse(daily >= prone_end & !is.na(prone_end), 1, from), # supine if beyond prone end date
                     to = ifelse(daily == (prone_start-1) , 2, 0), # to prone if delayed start
                     to = ifelse(daily == (prone_end-1) & !is.na(prone_end), 1, to), # to supine if prone end  
                     to = ifelse(daily == date_discharge, 3, to)) # discharge if last day, otherwise censored
        to_discharge_prone = bind_rows(to_discharge_prone, out)
      }
    }
    # check:
    # View(select(to_discharge_prone, pin, date_ecmo, date_discharge, prone_start, prone_end, prone_cum_start, prone_cum_end, prone_cum, start, end, from, to))
    
    # c1) times to censored in non-prone
    to_censored_non_prone = filter(indata, is.na(date_death)==TRUE, is.na(date_discharge)==TRUE, is.na(prone_start)==TRUE) %>%
      mutate(start_date = date_ecmo, # start from ecmo ...
             start_date = as.Date(start_date, origin='1970-01-01'), 
             prone_cum = 0, # cumulative days of prone, always none for this group
             from = 1, # ECMO and supine
             to = 0, # 'Censored',
             start = as.numeric(start_date - date_ecmo),
             end = as.numeric(last_date - date_ecmo)) # days to censored
    
    # c2) times to censored in prone
    prone_and_censored = filter(indata, is.na(date_death)==TRUE, is.na(date_discharge)==TRUE, is.na(prone_start)==FALSE) # 
    # have to loop
    to_censored_prone = NULL
    for (p in 1:nrow(prone_and_censored)){
      this_patient = prone_and_censored[p,]
      for(daily in this_patient$date_ecmo:this_patient$last_date){
        out = mutate(this_patient, 
                     start = as.numeric(as.Date(daily, origin='1970-01-01') - date_ecmo), 
                     end = start + 1, # just moving in days
                     prone_cum_start = (as.Date(daily, origin='1970-01-01') - (prone_start-1)) * as.numeric(as.Date(daily, origin='1970-01-01') > (prone_start-1)), # cumulative days of prone
                     prone_cum_end = (as.Date(daily, origin='1970-01-01') - (prone_end)) * as.numeric(as.Date(daily, origin='1970-01-01') >= (prone_end)), # cumulative days of prone
                     prone_cum_end = ifelse(is.na(prone_end)==TRUE, 0, prone_cum_end), # replace missing with none
                     prone_cum = prone_cum_start - prone_cum_end, # residual effect of prone ...
                     prone_cum = ifelse(as.Date(daily, origin='1970-01-01') >= prone_end & !is.na(prone_end), 0, prone_cum), # ... no residual effect of prone 
                     from = ifelse(daily >= prone_start, 2, 1), # prone if beyond prone date
                     from = ifelse(daily >= prone_end & !is.na(prone_end), 1, from), # supine if beyond prone end date
                     to = 0, # default is censored
                     to = ifelse(daily == (prone_start-1) , 2, to), # to prone if delayed start
                     to = ifelse(daily == (prone_end-1) & !is.na(prone_end), 1, to)) # to supine if prone end  
        to_censored_prone = bind_rows(to_censored_prone, out)
      }
    }
    # View(select(to_censored_prone, pin, date_ecmo, last_date, prone_start, prone_end, prone_cum_start, prone_cum_end, prone_cum, start, end, from, to))
    
  # combine six transitions 
  for_model = bind_rows(to_death_prone, to_death_non_prone, 
                        to_discharge_prone, to_discharge_non_prone,
                        to_censored_prone, to_censored_non_prone) %>%
    mutate(start = ifelse(start < 0, 0, start)) # assume negative dates mean ECMO on arrival - to check
  } # end of time-dependent if
  
  ## add half days
  # add half to prone if start and end are same
  index = with(for_model, from==1 & to==2 & start==end)
  for_model$end[index] = for_model$end[index] + 0.5
  # add half to supine if start and end are same
  index = with(for_model, from==2 & to==1 & start==end)
  for_model$end[index] = for_model$end[index] + 0.5
  # add half to discharge if start and end are same
  index = with(for_model, to==3 & start==end)
  for_model$end[index] = for_model$end[index] + 0.5
  # add half to death if start and end are same
  index = with(for_model, to==4 & start==end)
  for_model$end[index] = for_model$end[index] + 0.5
  # add half to censored if start and end are same
  index = with(for_model, to==0 & start==end)
  for_model$end[index] = for_model$end[index] + 0.5
  
  ## fixes because of same dates
  # if non-Prone to Prone at time zero then assume prone on entry (so no transition)
  index = with(for_model, from==1 & to==2 & start==0 & end <= 0.5)
  for_model = for_model[!index,]
  # if Prone to non-Prone at time zero then assume non-prone on entry (so no transition)
  index = with(for_model, from==2 & to==1 & start==0 & end <= 0.5)
  for_model = for_model[!index,]

  # remove impossible patients with obvious errors in dates
  exclude = filter(for_model, start > end) 
  for_model = filter(for_model,
                     !pin %in% exclude$pin) %>%
    arrange(pin, start, end)
  
  # remove any remaining NAs
  index = with(for_model,is.na(from)|is.na(to)|is.na(start)|is.na(end))
  for_model = for_model[!index,]
  
  # censor at a specified time (if it exists)
  if(is.na(censor_day) >= FALSE){
    for_model = filter(for_model, start < censor_day) %>% # only transitions before censor day
      mutate(beyond = end > censor_day,
             to = ifelse(beyond==TRUE, 0, to), # change state to censored
             end = ifelse(beyond==TRUE, censor_day, end))
  }
  
  # quick check
  if(check==TRUE){
    t1 = filter(for_model, from==1, to==2) %>%
      sample_n(1)
    t2 = filter(for_model, from==1, to==3) %>%
      sample_n(1)
    t3 = filter(for_model, from==1, to==4) %>%
      sample_n(1)
    t4 = filter(for_model, from==2, to==3) %>%
      sample_n(1)
    t5 = filter(for_model, from==2, to==4) %>%
      sample_n(1)
    t6 = filter(for_model, to==0) %>%
      sample_n(1)
    t7 = filter(for_model, to==1) %>%
      sample_n(1)
    test = bind_rows(t1, t2, t3, t4, t5, t6, t7)
    ids = unique(test$pin)
    cat('Original data:\n')
    patient_test = filter(indata, pin %in% ids) %>%
      select(pin, date_ecmo, prone_start, prone_end, date_discharge, date_death, last_date) %>%
      arrange(pin)
    print(patient_test)
    cat('Transitions:\n')
    transition_test = filter(for_model, pin %in% ids) %>%
      select(pin, from, to, start, end)
    print(transition_test)
    
  }
  
  return(for_model)
} # end of function



# function to change 0/1 to no/yes
change_factor_lvls <- function(x){
  return(factor(x,levels=0:1,labels=c('No','Yes')))
}

## convert additional date variables
convert_dates <- function(x){
  out = ifelse(x!='',as.Date(x,format='%Y-%m-%d'),NA)
  return(as.Date(out, origin='1970-01-01'))
}

## nice rename for tables
nice_rename = function(invar){
  invar = ifelse(str_detect(string=invar, pattern='^age$')==TRUE, 'Age', invar)
  invar = ifelse(str_detect(string=invar, pattern='age/5')==TRUE, 'Age (+5 Years)', invar)
  invar = ifelse(str_detect(string=invar, pattern='sexMale')==TRUE, 'Male', invar)
  invar = ifelse(str_detect(string=invar, pattern='fromSupine')==TRUE, 'Supine', invar)
  invar = ifelse(str_detect(string=invar, pattern='comorbidity_obesity')==TRUE, 'Comorbidity = obesity', invar)
  invar = ifelse(str_detect(string=invar, pattern='obese1')==TRUE, 'Obesity = Yes', invar)
  invar = ifelse(str_detect(string=invar, pattern='obese2')==TRUE, 'Obesity = Missing', invar)
  invar = ifelse(str_detect(string=invar, pattern='bmi_catLow')==TRUE, 'BMI = Below median', invar)
  invar = ifelse(str_detect(string=invar, pattern='bmi_catMissing')==TRUE, 'BMI = Missing', invar)
  invar = ifelse(str_detect(string=invar, pattern='apache_catLow')==TRUE, 'APACHE II = Below median', invar)
  invar = ifelse(str_detect(string=invar, pattern='apache_catMissing')==TRUE, 'APACHE II = Missing', invar)
  return(invar)
}

## extract results from WingBugs models
bugs_extract = function(inbugs, 
                        var, # variable to extract, e.g., `beta`
                        exp = FALSE, # for ratios, e.g., OR, HR, RR
                        names = NULL, # vector of names to add to rows
                        make.ci = FALSE, # make a confidence interval
                        digits = NULL, # round estimates to this many digits
                        arrange = NULL # arrange estimates from high to low mean
                        ){
  rows = data.frame(inbugs$summary[grep(var, rownames(inbugs$summary)), c(1,3,7)])
  names(rows) = c('mean','lower','upper')
  rows  = mutate(rows,
                 name = rownames(rows))
  type = str_detect(string=rows$name[1], pattern=',') # does it have rows and columns?
  brackets = str_detect(string=rows$name[1], pattern='\\[') # does it have brackets
  if(type == TRUE & brackets==TRUE){
    rows = tidyr::separate(rows, col='name', sep=',', into=c('Var','row','col','junk')) %>% # extract row and column number
      mutate(row = as.numeric(row),
             col = as.numeric(col)) %>%
      arrange(Var, row, col) %>%
      select(Var, row, col, mean, lower, upper) # arrange columns
  }
  if(type == FALSE & brackets==TRUE){ # just rows
    rows = tidyr::separate(rows, col='name', into=c('Var','row','junk')) %>%
      mutate(row = as.numeric(row)) %>%
      arrange(Var, row) %>%
      select(Var, row, mean, lower, upper) # arrange columns
  }
  if(brackets==FALSE){ # no rows or columns
    rows = tidyr::separate(rows, col='name', into=c('Var','junk')) %>%
      arrange(Var) %>%
      select(Var, mean, lower, upper) # arrange columns
  }
  # exponentiate for odds ratio, rate ratio, etc
  if(exp==TRUE){
    rows = mutate(rows,
                  mean = exp(mean),
                  lower = exp(lower),
                  upper = exp(upper))
  }
  # add variable names
  if(is.null(names)==FALSE){
    rows$variable = names
    rows = select(rows, variable, mean, lower, upper)
  }
  # add variable names
  if(is.null(arrange)==FALSE){
    rows = select(rows, -row) %>% # remove old row
      arrange(mean) %>%
      mutate(row = 1:n()) # add new row
  }
  # round mean and intervals
  if(is.null(digits)==FALSE){
    rows = mutate(rows,
                  mean = roundz(mean, digits),
                  lower = roundz(lower, digits),
                  upper = roundz(upper, digits))
  }
  # make CI 
  if(make.ci == TRUE){
    rows = mutate(rows,
                  CI = paste(lower, ' to ', upper, sep='')) %>%
      select(-lower, -upper) # drop intervals
  }
  #
  return(rows)
}

## make box diagram for multi-state models (uses diagram library)
multistate_diagram = function(states, links, pos, arr.pos=0.5){
  par(mai=rep(0.01, 4))
  N = length(states)
  M <- matrix(nrow = N, ncol = N, byrow = TRUE, data = 0)
  for (i in 1:nrow(links)){ # add links between states
    M[links[i,1], links[i,2]] = "' '"
  }
  plotmat(M, pos = pos, name = states, lwd = 1,
          box.lwd = 2, cex.txt = 0.8, box.size = 0.12, curve = 0,
          box.type = "square", box.prop = 0.5, arr.pos=arr.pos)
  
}

### Descriptive tables by group. # group_var must be "Yes" / "No"
# a) for continuous variables
desc_table_cont = function(indata, vars, group_var, heading){
    # simple rename
    index = names(indata) == group_var
    names(indata)[index] = 'group_var'
    # make table by group
    tab = select(indata, pin, group_var, all_of(vars)) %>%
    tidyr::gather(key='Variable', value='Value', -pin, -group_var) %>%
    group_by(Variable, group_var) %>%
    summarise(Missing = sum(is.na(Value)),
              Median = roundz(median(Value, na.rm=TRUE),0),
              Q1 = roundz(quantile(Value, probs=0.25, na.rm=TRUE),0),
              Q3 = roundz(quantile(Value, probs=0.75, na.rm=TRUE),0),
              IQR = paste(Q1 , ' to ', Q3, sep='')) %>%
    select(group_var, Variable, Missing, Median, IQR)
    # total column
    tab_total = select(indata, pin, all_of(vars)) %>%
      tidyr::gather(key='Variable', value='Value', -pin) %>%
      group_by(Variable) %>%
      summarise(Missing = sum(is.na(Value)),
                Median = roundz(median(Value, na.rm=TRUE),0),
                Q1 = roundz(quantile(Value, probs=0.25, na.rm=TRUE),0),
                Q3 = roundz(quantile(Value, probs=0.75, na.rm=TRUE),0),
                IQR = paste(Q1 , ' to ', Q3, sep='')) %>%
      select(Variable, Missing, Median, IQR) %>%
      mutate(group_var = 'Total')
    tab = bind_rows(tab, tab_total)
    # columns by group
  columns = c('Variable','No_Missing','No_Median','No_IQR','Yes_Missing','Yes_Median','Yes_IQR','Total_Missing','Total_Median','Total_IQR')
  ctab = tidyr::gather(tab, variable, value, -(group_var:Variable)) %>% # move stats into one column
    tidyr::unite(temp, group_var, variable) %>%
    tidyr::spread(temp, value) %>%
    select(columns)
  # table header
  table_header <- data.frame(
    col_keys = columns,
    h1 = c('', rep(heading[1], 3), rep(heading[2], 3), rep('Total', 3)),
    h2 = c('Variable','Missing','Median','IQR','Missing','Median','IQR','Missing','Median','IQR'),
    stringsAsFactors = FALSE )
  ftab = flextable(ctab) %>% 
    set_header_df(mapping = table_header, key = "col_keys" ) %>% # add header
    merge_h(i=1, part = "header") %>% # merge "To" column headers
    theme_box() %>%
    autofit()
  return(ftab)
}
# b) for categorical variables
desc_table_cat = function(indata, vars, group_var, heading){
  # simple rename
  index = names(indata) == group_var
  names(indata)[index] = 'group_var'
  # make table by group
  tab = select(indata, pin, group_var, all_of(vars)) %>%
    tidyr::gather(key='Variable', value='Value', -pin, -group_var) %>%
    mutate( # explicit missing:
      Value = ifelse(Value=='', 'Missing', Value),
      Value = ifelse(is.na(Value)==TRUE, 'Missing', Value),
      # nice variable names:
      Variable = ifelse(Variable=='ecmo_prone_before', 'Prone Position Before ECMO', Variable),
      Variable = ifelse(Variable=='pre_ecmo_vent', 'Ventilatory Mode Before Start Of ECMO', Variable)) %>%
    group_by(Variable, group_var, Value) %>%
    tally() %>%
    group_by(Variable, group_var) %>%
    mutate(Percent = roundz(prop.table(n)*100,0)) %>% # calculate percent
    ungroup()
  # make table for totals
  tab_total = select(indata, pin, all_of(vars)) %>%
    tidyr::gather(key='Variable', value='Value', -pin) %>%
    mutate( # explicit missing:
      Value = ifelse(Value=='', 'Missing', Value),
      Value = ifelse(is.na(Value)==TRUE, 'Missing', Value),
      # nice variable names:
      Variable = ifelse(Variable=='ecmo_prone_before', 'Prone Position Before ECMO', Variable),
      Variable = ifelse(Variable=='pre_ecmo_vent', 'Ventilatory Mode Before Start Of ECMO', Variable)) %>%
    group_by(Variable, Value) %>%
    tally() %>%
    group_by(Variable) %>%
    mutate(Percent = roundz(prop.table(n)*100,0)) %>% # calculate percent
    ungroup() %>%
    mutate(group_var = 'Total')
  tab = bind_rows(tab, tab_total) # add totals
  # columns for prone vs any prone
  columns = c('Variable','Value','No_n',"No_Percent",'Yes_n','Yes_Percent','Total_n','Total_Percent')
  ctab = tidyr::gather(tab, `n`, `Percent`, key='stat', value='result', -Variable, -group_var) %>% # move all stats to one column
    tidyr::unite(temp, group_var, stat) %>%
    group_by(Variable, Value) %>%
    tidyr::spread(temp, result) %>%
    select(columns) %>%
    ungroup() %>%
    mutate_if(is.character, tidyr::replace_na, replace = 0)
  # table header
  table_header <- data.frame(
    col_keys = columns,
    h1 = c('', '', rep(heading[1], 2), rep(heading[2], 2), rep('Total', 2)),
    h2 = c('Variable','Value','N','%','N','%','N','%'), stringsAsFactors = FALSE )
  ftab = flextable(ctab) %>% 
    set_header_df(mapping = table_header, key = "col_keys" ) %>% # add header
    merge_h(i=1, part = "header") %>% # merge "To" column headers
    merge_v(j = c("Variable")) %>% # merge rows
    theme_box() %>%
    autofit()
  return(ftab)
}
