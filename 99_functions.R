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
multistate_diagram = function(states, links, pos, arr.pos=0.5, box.size = 0.12){
  par(mai=rep(0.01, 4))
  N = length(states)
  M <- matrix(nrow = N, ncol = N, byrow = TRUE, data = 0)
  for (i in 1:nrow(links)){ # add links between states
    M[links[i,1], links[i,2]] = "' '"
  }
  plotmat(M, pos = pos, name = states, lwd = 1,
          box.lwd = 2, cex.txt = 0.8, box.size = box.size, curve = 0,
          box.type = "square", box.prop = 0.5, arr.pos=arr.pos)
  
}

### Descriptive tables by group. # group_var must be "Yes" / "No"
# a) for continuous variables
desc_table_cont = function(
  indata, 
  vars, 
  group_var, 
  var_levels,
  heading,
  autofit = TRUE,
  include_total = TRUE # add total column or not
){
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
    if(include_total==TRUE){
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
    }
    # nice variable names (continuous):
    tab = mutate(tab,
           Variable = ifelse(Variable=='duration_ecmo', 'Duration of ECMO', Variable),
           Variable = ifelse(Variable=='duration_mv', 'Duration of MV', Variable),
           Variable = ifelse(Variable=='days_hosp_stroke', 'Days from hospital admission to stroke', Variable),
           Variable = ifelse(Variable=='days_hosp', 'Days in hospital', Variable),
           Variable = ifelse(Variable=='days_symptoms_hosp', 'Days between symptoms and hospital', Variable),
           Variable = str_replace_all(Variable, pattern='_', replacement=' ')) 
    # columns in the right order:
    columns = c('Variable', paste(rep(var_levels, each=3), c('_Missing','_Median','_IQR'), sep=''))
    if(include_total == TRUE){columns = c(columns, 'Total_Missing','Total_Median','Total_IQR')}
    ctab = tidyr::gather(tab, variable, value, -(group_var:Variable)) %>% # move stats into one column
    tidyr::unite(temp, group_var, variable) %>%
    tidyr::spread(temp, value) %>%
    select(columns)
    # table header
    if(include_total == FALSE){
      table_header <- data.frame(
        col_keys = columns,
        h1 = c('', rep(heading, each=3)),
        h2 = c('Variable', rep(c('Missing','Median','IQR'), length(heading))))
    }
    if(include_total == TRUE){
      table_header <- data.frame(
        col_keys = columns,
        h1 = c('', rep(heading, each=3), rep('Total', 3)),
        h2 = c('Variable', rep(c('Missing','Median','IQR'), length(heading)+1)))
    }
  ftab = flextable(ctab) %>% 
    set_header_df(mapping = table_header, key = "col_keys" ) %>% # add header
    merge_h(i=1, part = "header") %>% # merge "To" column headers
    theme_box() 
  if(autofit==TRUE){ftab = autofit(ftab)}
  return(ftab)
}

# b) for categorical variables
desc_table_cat = function(indata, 
                          vars, 
                          group_var, 
                          var_levels = c('Yes','No'), # levels of the group variable
                          include_total = TRUE, # add total column or not
                          heading){
  # simple rename
  index = names(indata) == group_var
  names(indata)[index] = 'group_var'
  # make table by group
  tab = select(indata, pin, group_var, all_of(vars)) %>%
    tidyr::gather(key='Variable', value='Value', -pin, -group_var) %>%
    mutate( # explicit missing:
      Value = ifelse(Value=='', 'Missing', Value),
      Value = ifelse(is.na(Value)==TRUE, 'Missing', Value)) %>%
    group_by(Variable, group_var, Value) %>%
    tally() %>%
    group_by(Variable, group_var) %>%
    mutate(Percent = roundz(prop.table(n)*100,0)) %>% # calculate percent
    ungroup()
  # make table for totals
  if(include_total == TRUE){
  tab_total = select(indata, pin, all_of(vars)) %>%
    tidyr::gather(key='Variable', value='Value', -pin) %>%
    mutate( # explicit missing:
      Value = ifelse(Value=='', 'Missing', Value),
      Value = ifelse(is.na(Value)==TRUE, 'Missing', Value)) %>%
    group_by(Variable, Value) %>%
    tally() %>%
    group_by(Variable) %>%
    mutate(Percent = roundz(prop.table(n)*100,0)) %>% # calculate percent
    ungroup() %>%
    mutate(group_var = 'Total')
  tab = bind_rows(tab, tab_total) # add totals
  }
  # nice variable names (categorical):
  tab = mutate(tab,
  Variable = str_remove_all(Variable, pattern='comorbidity_'),
  Variable = ifelse(Variable=='ecmo_prone_before', 'Prone Position Before ECMO', Variable),
  Variable = ifelse(Variable=='pre_ecmo_vent', 'Ventilatory Mode Before Start Of ECMO', Variable))
  # columns in the right order:
  columns = c('Variable','Value', paste(rep(var_levels, each=2), c('_n','_Percent'), sep=''))
  if(include_total == TRUE){columns = c(columns, 'Total_n','Total_Percent')}
  ctab = tidyr::gather(tab, `n`, `Percent`, key='stat', value='result', -Variable, -group_var) %>% # move all stats to one column
    tidyr::unite(temp, group_var, stat) %>%
    group_by(Variable, Value) %>%
    tidyr::spread(temp, result) %>%
    select(columns) %>%
    ungroup() %>%
    mutate_if(is.character, tidyr::replace_na, replace = 0)
  # table header
  if(include_total == FALSE){
  table_header <- data.frame(
    col_keys = columns,
    h1 = c('', '', rep(heading, each=2)),
    h2 = c('Variable','Value', rep(c('N','%'), length(heading))))
  }
  if(include_total == TRUE){
    table_header <- data.frame(
      col_keys = columns,
      h1 = c('', '', rep(heading, each=2), rep('Total', 2)),
      h2 = c('Variable','Value', rep(c('N','%'), length(heading)+1)))
  }
  ftab = flextable(ctab) %>% 
    set_header_df(mapping = table_header, key = "col_keys" ) %>% # add header
    merge_h(i=1, part = "header") %>% # merge "To" column headers
    merge_v(j = c("Variable")) %>% # merge rows
    theme_box() %>%
    autofit()
  return(ftab)
}

### function used by time_transitions to randomly sample data to check
sample_transitions = function(indata, in_from, in_to){
  transitions = filter(indata, from %in% in_from, to %in% in_to) 
  if(nrow(transitions) > 0){
    transitions = sample_n(transitions, 1)
  }
  if(nrow(transitions) ==0){
    transitions = NULL
  }
  return(transitions)
}



### function to make transitions times with a general intermediate state
time_transitions = function(indata, 
                            start_state = 'date_icu', # starting date/state
                            int_state = 'date_mechanical_ventilation', # name of intermediate state
                            int_state_end = NULL, # name of intermediate state end (optional)
                            discharge_date = 'date_hospital_discharge', # ICU (date_discharge) or hospital discharge (date_hospital_discharge)
                            check = TRUE, 
                            censor_day = NA){
  
  # State 0: Censored
  # State 1: starting state (could be ICU admission or ECMO start)
  # State 2: intermediate state
  # State 3: Discharge (ICU or hospital)
  # State 4: Death
  
  # simple rename for starting state
  index = names(indata) == start_state
  names(indata)[index] = 'date_start'
  if(sum(index)==0){cat('Warning, check name starting date\n')}
  # simple rename for intermediate state
  index = names(indata) == int_state
  names(indata)[index] = 'date_int_state'
  if(sum(index)==0){cat('Warning, check name intermediate date\n')}
  if(is.null(int_state_end) == FALSE){
    index = names(indata) == int_state_end
    names(indata)[index] = 'date_int_state_end'
    if(sum(index)==0){cat('Warning, check name intermediate end date\n')}
    
  }
  # simple rename for discharge
  index = names(indata) == discharge_date
  names(indata)[index] = 'discharge_date'
  if(sum(index)==0){cat('Warning, check name discharge date\n')}
  
  # exclude those with no event dates 
  if(is.null(int_state_end) == TRUE){
    index1 = is.na(indata$date_int_state) & 
      is.na(indata$date_death) &
      is.na(indata$discharge_date) &
      is.na(indata$last_date)
  }
  if(is.null(int_state_end) == FALSE){
    index1 = is.na(indata$date_int_state) & 
      is.na(indata$date_death) &
      is.na(indata$discharge_date) &
      is.na(indata$last_date)&
      is.na(indata$date_int_state_end)
  }
  index2 = is.na(indata$date_start) # must have start date
  index = index1 | index2
  n_excluded = sum(index)
  indata = indata[!index,]
  
  # if intermediate state stopped is same as death or discharge then delete intermediate end date
  if(is.null(int_state_end) == FALSE){
    indata = mutate(indata,
                    date_int_state_end = ifelse(date_int_state_end == discharge_date & !is.na(discharge_date), NA, date_int_state_end),
                    date_int_state_end = ifelse(date_int_state_end == date_death & !is.na(date_death), NA, date_int_state_end),
                    date_int_state_end = as.Date(date_int_state_end, origin='1970-01-01'))
  }
  
  # a) times to death
  to_death = filter(indata, is.na(date_death) == FALSE) # must have date of death
  if(is.null(int_state_end) == TRUE){ # if no intermediate state
    to_death = mutate(to_death, 
                      start_date = ifelse(is.na(date_int_state) == FALSE, date_int_state, date_start), # start from intermediate or admission
                      from = ifelse(is.na(date_int_state) == FALSE, 2, 1)) # if intermediate ended then back to 1
  }
  if(is.null(int_state_end) == FALSE){ # if intermediate state
    to_death = mutate(to_death, 
                      start_date = ifelse(is.na(date_int_state) == FALSE, date_int_state, date_start),# start from intermediate or admission ...
                      start_date = ifelse(is.na(date_int_state_end) == FALSE, date_int_state_end, start_date), # ... update to intermediate end if it's there
                      from = case_when(
                        is.na(date_int_state) == TRUE & is.na(date_int_state_end) == TRUE ~ 1, # not intermediate
                        is.na(date_int_state) == FALSE & is.na(date_int_state_end) == TRUE ~ 2, # still in intermediate
                        is.na(date_int_state) == TRUE & is.na(date_int_state_end) == FALSE ~ 1, # back to start
                        is.na(date_int_state) == FALSE & is.na(date_int_state_end) == FALSE ~ 1 # back to start
                      ))
  }
  to_death = mutate(to_death, 
                    start_date = as.Date(start_date, origin='1970-01-01'), 
                    to = 4, # death
                    start = as.numeric(start_date - date_start), # start time in days
                    end = as.numeric(date_death - date_start)) # end time in days (to death)
  
  # b) times to discharge
  to_discharge = filter(indata, is.na(date_death) == TRUE, is.na(discharge_date) == FALSE)  # must have no death date but discharge date
  if(is.null(int_state_end) == TRUE){ # if no intermediate state
    to_discharge = mutate(to_discharge, 
                          start_date = ifelse(is.na(date_int_state) == FALSE, date_int_state, date_start), # start from intermediate or admission
                          from = ifelse(is.na(date_int_state) == FALSE, 2, 1)) # if intermediate ended then back to 1
  }
  if(is.null(int_state_end) == FALSE){ # if intermediate state
    to_discharge = mutate(to_discharge, 
                          start_date = ifelse(is.na(date_int_state) == FALSE, date_int_state, date_start),# start from intermediate or admission ...
                          start_date = ifelse(is.na(date_int_state_end) == FALSE, date_int_state_end, start_date), # ... update to intermediate end if it's there
                          from = case_when(
                            is.na(date_int_state) == TRUE & is.na(date_int_state_end) == TRUE ~ 1, # not intermediate
                            is.na(date_int_state) == FALSE & is.na(date_int_state_end) == TRUE ~ 2, # still in intermediate
                            is.na(date_int_state) == TRUE & is.na(date_int_state_end) == FALSE ~ 1, # back to start
                            is.na(date_int_state) == FALSE & is.na(date_int_state_end) == FALSE ~ 1 # back to start
                          ))
  }
  to_discharge = mutate(to_discharge, 
                        start_date = as.Date(start_date, origin='1970-01-01'), 
                        to = 3, # discharge
                        start = as.numeric(start_date - date_start), # start time in days
                        end = as.numeric(discharge_date - date_start)) # end time in days (to discharge)
  
  # c) times to censored
  to_censored = filter(indata, is.na(date_death) == TRUE, is.na(discharge_date) == TRUE) 
  if(is.null(int_state_end) == TRUE){ # if no intermediate state
    to_censored = mutate(to_censored, 
                         start_date = ifelse(is.na(date_int_state) == FALSE, date_int_state, date_start), # start from intermediate or admission
                         from = ifelse(is.na(date_int_state) == FALSE, 2, 1)) # if intermediate ended then back to 1
  }
  if(is.null(int_state_end) == FALSE){ # if intermediate state
    to_censored = mutate(to_censored, 
                         start_date = ifelse(is.na(date_int_state) == FALSE, date_int_state, date_start),# start from intermediate or admission ...
                         start_date = ifelse(is.na(date_int_state_end) == FALSE, date_int_state_end, start_date), # ... update to intermediate end if it's there
                         from = case_when(
                           is.na(date_int_state) == TRUE & is.na(date_int_state_end) == TRUE ~ 1, # not intermediate
                           is.na(date_int_state) == FALSE & is.na(date_int_state_end) == TRUE ~ 2, # still in intermediate
                           is.na(date_int_state) == TRUE & is.na(date_int_state_end) == FALSE ~ 1, # back to start
                           is.na(date_int_state) == FALSE & is.na(date_int_state_end) == FALSE ~ 1 # back to start
                         ))
  }
  to_censored = mutate(to_censored, 
                       start_date = as.Date(start_date, origin='1970-01-01'), 
                       to = 0, # censored
                       start = as.numeric(start_date - date_start), # start time in days
                       end = as.numeric(last_date - date_start)) # end time in days (to discharge)
  
  # d) times to intermediate (from starting state)
  to_intermediate = filter(indata, is.na(date_int_state) == FALSE) %>%
    mutate(start_date = date_start,
           start_date = as.Date(start_date, origin='1970-01-01'), 
           from = 1, # Non-intermediate
           to = 2, # intermediate
           start = 0, 
           end = as.numeric(date_int_state - date_start), # days to intermediate
           end = ifelse(end < 0, 0, end)) # if negative assume ventilated on admission
  
  # e) times to non-intermediate (from intermediate)
  if(is.null(int_state_end) == FALSE){ # if end date
    to_non_intermediate = filter(indata, is.na(date_int_state_end) == FALSE) %>%
      mutate(start_date = date_int_state,
             end_date = date_int_state_end,
             from = 2, # intermediate
             to = 1, # non-intermediate
             start = as.numeric(start_date - date_start), # days to intermediate start
             end = as.numeric(end_date - date_start), # days to intermediate end
             start = ifelse(start<0, 0, start)) # assume ventilated on arrival if negative
  }
  
  # combine four/five transitions 
  for_model = bind_rows(to_death, to_discharge, to_censored, to_intermediate)
  if(is.null(int_state_end) == FALSE){ # if end date
    for_model = bind_rows(for_model, to_non_intermediate)
  }
  for_model = mutate(for_model, start = ifelse(start < 0, 0, start)) # assume negative dates mean ventilated on arrival
  
  ## fixes because of same dates
  # if move from transition to non-transition before start (so no transition)
  index = with(for_model, from == 1 & to == 2 & end < 0)
  for_model = for_model[!index,]
  index = with(for_model, from == 2 & to == 1 & end < 0)
  for_model = for_model[!index,]
  # if non-intermediate to intermediate at time zero then assume intermediate on entry (so no transition)
  index = with(for_model, from == 1 & to == 2 & start == 0 & end == 0)
  for_model = for_model[!index,]
  # if intermediate to non-intermediate at time zero then assume intermediate on entry (so no transition)
  index = with(for_model, from == 2 & to == 1 & start == 0 & end == 0)
  for_model = for_model[!index,]
  
  ## add half days
  # add half to discharge if start and end are same
  index = with(for_model, from == 1 & to == 2 & start == end)
  for_model$end[index] = for_model$end[index] + 0.5
  # add half to discharge if start and end are same
  index = with(for_model, from == 2 & to == 1 & start == end)
  for_model$end[index] = for_model$end[index] + 0.5
  # add half to discharge if start and end are same
  index = with(for_model, to == 3 & start == end)
  for_model$end[index] = for_model$end[index] + 0.5
  # add half to death if start and end are same
  index = with(for_model, to == 4 & start == end)
  for_model$end[index] = for_model$end[index] + 0.5
  # add half to censored if start and end are same
  index = with(for_model, to == 0 & start == end)
  for_model$end[index] = for_model$end[index] + 0.5
  
  # remove impossible patients with obvious errors in dates
  exclude = filter(for_model, start > end) 
  for_model = filter(for_model,
                     !pin %in% exclude$pin) %>%
    arrange(pin, start, end)
  
  #remove any remaining NAs
  index = with(for_model, is.na(from)|is.na(to)|is.na(start)|is.na(end))
  for_model = for_model[!index,]
  
  # censor at a specified time (if it exists)
  if(is.na(censor_day) >= FALSE){
    for_model = filter(for_model, start < censor_day) %>% # only transitions before censor day
      mutate(beyond = end > censor_day,
             to = ifelse(beyond == TRUE, 0, to), # change state to censored
             end = ifelse(beyond == TRUE, censor_day, end))
  }
  
  # quick check
  if(check == TRUE){
    t1 = sample_transitions(for_model, in_from = 1, in_to = 2) 
    t2 = sample_transitions(for_model, in_from = 1, in_to = 3) 
    t3 = sample_transitions(for_model, in_from = 1, in_to = 4) 
    t4 = sample_transitions(for_model, in_from = 2, in_to = 3) 
    t5 = sample_transitions(for_model, in_from = 2, in_to = 4) 
    t6 = sample_transitions(for_model, in_from = c(1,2), in_to = 4) 
    t7 = sample_transitions(for_model, in_from = 1, in_to = c(2,3,4)) 
    test = bind_rows(t1, t2, t3, t4, t5, t6, t7)
    ids = unique(test$pin)
    cat('Original data:\n')
    patient_test = filter(indata, pin %in% ids) %>%
      select(pin, date_start, contains('date_int_'), discharge_date, date_death, last_date) %>%
      arrange(pin)
    print(patient_test)
    cat('Transitions:\n')
    transition_test = filter(for_model, pin %in% ids) %>%
      select(pin, from, to, start, end)
    print(transition_test)
    
  }
  
  # rename back
  index = names(for_model) == 'date_start'
  names(for_model)[index] = start_state
  # simple rename for intermediate state
  index = names(for_model) == 'date_int_state'
  names(for_model)[index] = int_state
  if(is.null(int_state_end) == FALSE){
    index = names(for_model) == 'date_int_state_end'
    names(for_model)[index] = int_state_end
  }
  
  return(for_model)
} # end of function

