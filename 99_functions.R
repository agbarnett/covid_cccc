# 99_functions.R
# range of functions used by multiple programs 
# November 2020

# mode
Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

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
    select(pin, stroke_group, date_ecmo, before)
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
           Variable = ifelse(Variable=='days_icu', 'Days in ICU', Variable),
           Variable = ifelse(Variable=='days_admission_mv', 'Days from hospital admission to MV', Variable),
           Variable = ifelse(Variable=='days_admission_ecmo', 'Days from hospital admission to ECMO', Variable),
           Variable = ifelse(Variable=='days_symptoms_hosp', 'Days between symptoms and hospital admission', Variable),
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
  if(is.null(int_state)==FALSE){
    index = names(indata) == int_state
    names(indata)[index] = 'date_int_state'
    if(sum(index)==0){cat('Warning, check name intermediate date\n')}
  }
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
                       end = as.numeric(last_date - date_start)) %>% # end time in days (to discharge)
    filter(!is.na(end)) # remove if there's no end date (no date to base their censoring on)
  
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
  
  #
  #c = select(for_model, pin, date_start, date_int_state, discharge_date, from, to, start, end) %>%
  #  filter(is.na(start))
  
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

### function to make transitions times with no intermediate state
time_transitions_no_intermediate = function(indata, 
                                            start_state = 'date_icu', # starting date/state
                                            discharge_date = 'date_hospital_discharge', # ICU (date_discharge) or hospital discharge (date_hospital_discharge)
                                            check = TRUE, 
                                            censor_day = NA){
  
  # State 0: Censored
  # State 1: starting state (could be ICU admission or ECMO start)
  # State 2: Discharge (ICU or hospital)
  # State 3: Death
  
  # simple rename for starting state
  index = names(indata) == start_state
  names(indata)[index] = 'date_start'
  if(sum(index)==0){cat('Warning, check name starting date\n')}
  # simple rename for discharge
  index = names(indata) == discharge_date
  names(indata)[index] = 'discharge_date'
  if(sum(index)==0){cat('Warning, check name discharge date\n')}
  
  ## fix dates
  ## cannot tell what happened after transfers, so blank dates
  # if transfer date on discharge date then blank transfer
  index = indata$date_transfer == indata$discharge_date & !is.na(indata$discharge_date) & !is.na(indata$date_transfer)
  indata$date_transfer[index] = NA
  # if transfer date on death date then blank transfer
  index = indata$date_transfer == indata$date_death & !is.na(indata$date_death) & !is.na(indata$date_transfer)
  indata$date_transfer[index] = NA
  # if transfer date before start then blank start
  index = indata$date_transfer <= indata$date_start & !is.na(indata$date_start) & !is.na(indata$date_transfer)
  indata$start_state[index] = NA
  # if transfer date before discharge then blank discharge
  index = indata$date_transfer <= indata$discharge_date & !is.na(indata$discharge_date) & !is.na(indata$date_transfer)
  indata$discharge_date[index] = NA
  # if transfer date before death then blank death
  index = indata$date_transfer <= indata$date_death & !is.na(indata$date_death) & !is.na(indata$date_transfer)
  indata$date_death[index] = NA
  # transfer date also trumps last date
  index = indata$date_transfer <= indata$last_date & !is.na(indata$last_date) & !is.na(indata$date_transfer)
  indata$last_date[index] = NA
  # if death and discharge (ICU or hospital) on same day then blank discharge
  index = indata$date_death == indata$discharge_date & !is.na(indata$date_death)
  if(length(index)>0){indata$discharge_date[index] = NA}
  # if last date before known dates then blank last date (assume it is an error)
  for (d in c(start_state, discharge_date, 'date_death')){
    # simple rename
    i = names(indata) == d 
    names(indata)[i] = 'date_check'
    # blank
    index = indata$date_check >= indata$last_date & !is.na(indata$last_date)
    if(length(index)>0){indata$last_date[index] = NA}
    # rename back
    names(indata)[i] = d
  }
  
  # exclude those with no event dates 
  index1 = is.na(indata$date_death) &
    is.na(indata$discharge_date) &
    is.na(indata$last_date)
  index2 = is.na(indata$date_start) # must have start date
  index = index1 | index2
  n_excluded = sum(index)
  indata = indata[!index,]
  
  # a) times to death
  to_death = filter(indata, is.na(date_death) == FALSE) %>% # must have date of death
    mutate(start_date = date_start, # start from admission
           start_date = as.Date(start_date, origin='1970-01-01'), 
           from = 1,
           to = 3, # death
           start = as.numeric(start_date - date_start), # start time in days
           end = as.numeric(date_death - date_start)) # end time in days (to death)
  
  # b) times to discharge
  to_discharge = filter(indata, is.na(date_death) == TRUE, is.na(discharge_date) == FALSE) %>% # must have no death date but discharge date
    mutate( start_date = date_start,
            start_date = as.Date(start_date, origin='1970-01-01'), 
            from = 1,
            to = 2, # discharge
            start = as.numeric(start_date - date_start), # start time in days
            end = as.numeric(discharge_date - date_start)) # end time in days (to discharge)
  
  # c) times to censored
  ## include date_transfer here too???
  to_censored = filter(indata, is.na(date_death) == TRUE, is.na(discharge_date) == TRUE) %>%
    mutate(start_date = date_start, # start from admission
           start_date = as.Date(start_date, origin='1970-01-01'), 
           from = 1,
           to = 0, # censored
           start = as.numeric(start_date - date_start), # start time in days
           end = as.numeric(last_date - date_start)) %>% # end time in days (to discharge)
    filter(!is.na(end)) # remove if there's no end date (no date to base their censoring on)
  
  # combine transitions 
  for_model = bind_rows(to_death, to_discharge, to_censored)
  for_model = mutate(for_model, start = ifelse(start < 0, 0, start)) # assume negative dates mean ventilated on arrival
  
  ## add half days
  # add half to discharge if start and end are same
  index = with(for_model, to == 2 & start == end)
  for_model$end[index] = for_model$end[index] + 0.5
  # add half to death if start and end are same
  index = with(for_model, to == 3 & start == end)
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
    t3 = sample_transitions(for_model, in_from = 1, in_to = 0)  # censored
    test = bind_rows(t1, t2, t3)
    ids = unique(test$pin)
    cat('Original data:\n')
    patient_test = filter(indata, pin %in% ids) %>%
      select(pin, date_start, discharge_date, date_death, date_transfer, last_date) %>%
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
  
  return(for_model)
} # end of function


## function to make multi-state survival analysis given string of dates ##
# data must have a bunch of dates, including `date_transfer` #
time_transitions_string = function(indata, 
                 intermediate_prior_ongoing = FALSE, # are intermediate events prior start date considered to be ongoing on start date (e.g., prone position); keep FALSE for acute events (e.g., stroke)
                 dates,  # vector of dates
                 int_state, # start of intermediate state
                 int_state_end = NULL, # end of intermediate state (optional)
                 start_state, # name of starting state
                 discharge_date, # name of discharge (hospital or ICU)
                 censor_day,
                 check=FALSE){
  
  ## variable fixes
  # simple rename for starting state
  index = names(indata) == start_state
  names(indata)[index] = 'date_start'
  if(sum(index)==0){cat('Warning, check name start date\n')}
  # simple rename for discharge
  index = names(indata) == discharge_date
  names(indata)[index] = 'discharge_date'
  if(sum(index)==0){cat('Warning, check name discharge date\n')}
  # and add old variable back for starting date
  indata$temporary = indata$date_start
  index = names(indata) == 'temporary'
  names(indata)[index] = start_state
  # simple renames for intermediate state start and end
  index = names(indata) == int_state
  names(indata)[index] = 'date_int_state_start'
  if(sum(index)==0){cat('Warning, check name intermediate end date\n')}
  if(is.null(int_state_end) == FALSE){
    index = names(indata) == int_state_end
    names(indata)[index] = 'date_int_state_end'
    if(sum(index)==0){cat('Warning, check name intermediate end date\n')}
  }
  # update dates vector with generic discharge name
  dates[dates==discharge_date] = 'discharge_date'
  dates[dates==int_state] = 'date_int_state_start'
  if(is.null(int_state_end) == FALSE){
    dates[dates==int_state_end] = 'date_int_state_end'
  }
  
  # exclude patients without start date
  indata = filter(indata, !is.na(date_start))

  ## fix dates
  ## cannot tell what happened after transfers, so blank dates
  # if transfer date on discharge date then blank transfer
  index = indata$date_transfer == indata$discharge_date & !is.na(indata$discharge_date) & !is.na(indata$date_transfer)
  indata$date_transfer[index] = NA
  # if transfer date on death date then blank transfer
  index = indata$date_transfer == indata$date_death & !is.na(indata$date_death) & !is.na(indata$date_transfer)
  indata$date_transfer[index] = NA
  # if transfer date before ECMO then blank ECMO
  index = indata$date_transfer <= indata$date_ecmo & !is.na(indata$date_ecmo) & !is.na(indata$date_transfer)
  indata$date_ecmo[index] = NA
  # if transfer date before ECMO end then blank ECMO end
  index = indata$date_transfer <= indata$date_ecmo_discontinued & !is.na(indata$date_ecmo_discontinued) & !is.na(indata$date_transfer)
  indata$date_ecmo_discontinued[index] = NA
  # if transfer date before discharge then blank discharge
  index = indata$date_transfer <= indata$discharge_date & !is.na(indata$discharge_date) & !is.na(indata$date_transfer)
  indata$discharge_date[index] = NA
  # if transfer date before death then blank death
  index = indata$date_transfer <= indata$date_death & !is.na(indata$date_death) & !is.na(indata$date_transfer)
  indata$date_death[index] = NA
  # if transfer date before intermediate start
  index = indata$date_transfer <= indata$date_int_state_start & !is.na(indata$date_int_state_start) & !is.na(indata$date_transfer)
  indata$date_int_state_start[index] = NA
  # if transfer date before intermediate end
  if(is.null(int_state_end) == FALSE){
    index = indata$date_transfer <= indata$date_int_state_end & !is.na(indata$date_int_state_end) & !is.na(indata$date_transfer)
    indata$date_int_state_end[index] = NA
  }
  # transfer date also trumps last date
  index = indata$date_transfer <= indata$last_date & !is.na(indata$last_date) & !is.na(indata$date_transfer)
  indata$last_date[index] = NA
  ##
  # if death and discharge (ICU or hospital) on same day then blank discharge
  index = indata$date_death == indata$discharge_date & !is.na(indata$date_death)
  if(length(index)>0){indata$discharge_date[index] = NA}
  # if last date before known dates then blank
  for (d in dates){
    # simple rename
    i = names(indata) == d 
    names(indata)[i] = 'date_check'
    # blank
    index = indata$date_check >= indata$last_date & !is.na(indata$last_date)
    if(length(index)>0){indata$last_date[index] = NA}
    # rename back
    names(indata)[i] = d
  }
  # if intermediate state stopped is same as death or discharge then delete intermediate end date
  if(is.null(int_state_end) == FALSE){
    indata = mutate(indata,
         date_int_state_end = ifelse(date_int_state_end == discharge_date & !is.na(discharge_date), NA, date_int_state_end),
         date_int_state_end = ifelse(date_int_state_end == date_death & !is.na(date_death), NA, date_int_state_end),
         date_int_state_end = as.Date(date_int_state_end, origin='1970-01-01'))
  }
  # move intermediate start date to overall start date if it is considered an ongoing event
  if(intermediate_prior_ongoing == TRUE & is.null(int_state_end) == FALSE){
    indata = mutate(indata,
        date_int_state_start = ifelse(date_int_state_start < date_start & (date_int_state_end>date_start | is.na(date_int_state_end)) & !is.na(date_int_state_start)
                                      , date_start, date_int_state_start),
        date_int_state_start = as.Date(date_int_state_start, origin='1970-01-01'))
  }

  # arrange dates
  first = select(indata, pin, date_start) # firs
  dates_long = select(indata, pin, all_of(dates), 'last_date') %>%
    gather(key='event', value='date', -`pin`) %>%
    filter(!is.na(date)) %>% # remove any non events
    left_join(first, by='pin') %>% # add start date
    filter(date >= date_start) %>% # exclude any events prior to start date
    arrange(pin, date) %>% # sort by date
    mutate(time = as.numeric(date - date_start)) # get time difference
  # now loop through patients to get one step ahead times
  time_dep_data = NULL
  for (u in unique(dates_long$pin)){
    this_patient = filter(dates_long, pin == u)
    n = nrow(this_patient)
    #if(n == 1){cat(this_patient$pin, '\n')}
    if(n > 1){ # only for patients with at least one transition
      for (r in 1:(n-1)){
        f = data.frame(pin = u, 
                       from = this_patient$event[r], 
                       to = this_patient$event[r+1], 
                       entry = this_patient$time[r], 
                       exit = this_patient$time[r+1],
                       entry_date = this_patient$date[r], # add dates for joint models
                       exit_date = this_patient$date[r+1])
        time_dep_data = bind_rows(time_dep_data, f)
      }
    }
  }
  # edits to transitions:
  time_dep_data = filter(time_dep_data, 
          # remove impossible transitions from absorbing states
          !from %in% c('discharge_date','last_date','date_death'),
          # remove transitions after censor day
          entry <= censor_day) 
  # censor long times
  index = time_dep_data$exit > censor_day
  if(length(index)>0){
    time_dep_data$exit[index] = censor_day # change time
    time_dep_data$to[index] = 'last_date' # change state
  }
  # if non-intermediate to intermediate at time zero then assume intermediate on entry (so no transition)
  index = with(time_dep_data, from %in% c(start_state,'date_int_state_start') & to %in% c('date_int_state_start','date_int_state_end') & entry == 0 & exit == 0)
  if(length(index)>0){time_dep_data = time_dep_data[!index,]}
  # if intermediate to non-intermediate at time zero then assume intermediate on entry (so no transition)
  index = with(time_dep_data, from == 'date_int_state_end' & to == 'date_int_state_start' & entry == 0 & exit == 0)
  if(length(index)>0){time_dep_data = time_dep_data[!index,]}
  # add half day for transitions on the same day
  index = time_dep_data$entry == time_dep_data$exit
  time_dep_data$exit[index] = time_dep_data$exit[index] + 0.5
  with(time_dep_data, table(from, to)) # check
  
  # randomly check 5 patients
  if(check==TRUE){
    s = sample(indata$pin, replace=FALSE, size=5)
    d = dates[dates%in%start_state == FALSE] # drop start date, renamed
    #
    original = filter(indata, pin%in%s) %>% 
      select(pin, 'date_start', all_of(d), 'last_date')
    #
    time_dep = filter(time_dep_data, pin%in%s)
  }
  #
  #filter(patients, pin=='OX_724-0011') %>% select(contains('date_'), 'last_date')
  
  # add back patient information
  time_dep_data = left_join(time_dep_data, indata, by='pin') %>%
    mutate(to = ifelse(to=='date_death', 'death', to), # tidy up names
           to = ifelse(to=='discharge_date', 'discharge', to),
           to = ifelse(to=='last_date', 'censored', to)) # rename for censored
  
  return(time_dep_data)
  
} # end of function


## function to run case-control analysis of platelets/anticoagulants depending on stroke type
case_control = function(s_type, 
                        number_of_controls = 5
                        ){
  set.seed(212121) # for random control selection
  # split into cases and controls by stroke
  cases = filter(patients, 
                 !is.na(complication_stroke_date), # has a stroke ...
                 stroke_group %in% s_type) %>% # ... of this type
    select(pin, age, complication_stroke_date, date_icu) %>%
    mutate(days = as.numeric(complication_stroke_date - date_icu)) %>%
    filter(!is.na(days),
           days > 0) %>% # stroke must be after ICU admission
    mutate(days_c = case_when( # combine days as they get sparse
      days > 20 & days<=30 ~ 25,
      days > 30 & days<=40 ~ 35,
      days > 40 & days<=50 ~ 45,
      days > 50 ~ 55,
      TRUE ~ days
    ))
  controls = filter(patients, 
                    !pin %in% unique(cases$pin), # not in the stroke case data
                    !is.na(date_icu)) %>% # with ICU data
    select(pin, age, date_icu)
  
  ## merge cases with daily data up to day of stroke
  cases = mutate(cases, 
                 eotd_anticoagulants = NA, # start with blanks
                 eotd_anticoagulants_type = NA,
                 platelet_count = NA)
  for (i in 1:nrow(cases)){
    # a) anticoagulants yes/no
    this_daily = filter(daily, # get daily data
                        pin == cases$pin[i], # same patient
                        date_daily <= cases$complication_stroke_date[i], # dates prior to stroke
                        !is.na(eotd_anticoagulants)) %>% # exclude missing
      summarise(eotd_anticoagulants = Modes(eotd_anticoagulants)) # find mode
    # if equal modes take most recent 
    if(nrow(this_daily) > 1){
      this_daily_new = filter(daily, 
                              pin == cases$pin[i], # same patient
                              date_daily <= cases$complication_stroke_date[i],
                              eotd_anticoagulants %in% this_daily$eotd_anticoagulants) %>%
        arrange(date_daily) %>%
        slice(n()) # last date
      this_daily$eotd_anticoagulants = this_daily_new$eotd_anticoagulants # replaces both!
    }
    cases$eotd_anticoagulants[i] = this_daily$eotd_anticoagulants[1]
    # b) anticoagulants type
    this_daily = filter(daily, # get daily data
                        pin == cases$pin[i], # same patient
                        date_daily <= cases$complication_stroke_date[i], # dates prior to stroke
                        !is.na(eotd_anticoagulants_type)) %>% # exclude missing
      summarise(eotd_anticoagulants_type = Modes(eotd_anticoagulants_type)) # find mode
    # if equal modes take most recent 
    if(nrow(this_daily) > 1){
      this_daily_new = filter(daily, 
                              pin == cases$pin[i], # same patient
                              date_daily <= cases$complication_stroke_date[i],
                              eotd_anticoagulants_type %in% this_daily$eotd_anticoagulants_type) %>%
        arrange(date_daily) %>%
        slice(n()) # latest result
      this_daily$eotd_anticoagulants_type = this_daily_new$eotd_anticoagulants_type # replaces both!
    }
    cases$eotd_anticoagulants_type[i] = this_daily$eotd_anticoagulants_type[1]
    # platelets
    this_daily = filter(daily, # get daily data
                        pin == cases$pin[i], # same patient
                        date_daily == cases$complication_stroke_date[i]) # same day
    if(nrow(this_daily) > 0){ cases$platelet_count[i] = this_daily$platelet_count}
  }
  
  ## merge controls with daily data (take all days)
  controls = left_join(controls, daily, by=c('pin'='pin')) %>%
    mutate(days = as.numeric(date_daily - date_icu)) %>%
    filter(!is.na(days),
           days > 0) %>%
    mutate(days_c = case_when( # combine days as they get sparse
      days > 20 & days<=30 ~ 25,
      days > 30 & days<=40 ~ 35,
      days > 40 & days<=50 ~ 45,
      days > 50 ~ 55,
      TRUE ~ days
    ))
  # get eotd_anticoagulants up to day of stroke
  controls = mutate(controls, 
                    eotd_anticoagulants = NA, # start with blanks
                    eotd_anticoagulants_type = NA)
  for (i in 1:nrow(controls)){ # takes a while
    # a) yes/no
    this_daily = filter(daily, # get daily data
                        pin == controls$pin[i], # same patient
                        date_daily <= controls$date_daily[i], # dates prior to stroke
                        !is.na(eotd_anticoagulants)) %>% # exclude missing
      summarise(eotd_anticoagulants = Modes(eotd_anticoagulants)) # find mode
    # if equal modes take most recent 
    if(nrow(this_daily) > 1){
      this_daily_new = filter(daily, 
                              pin == controls$pin[i], # same patient
                              date_daily <= controls$date_daily[i],
                              eotd_anticoagulants %in% this_daily$eotd_anticoagulants) %>%
        arrange(date_daily) %>%
        slice(n()) # latest result
      this_daily$eotd_anticoagulants = this_daily_new$eotd_anticoagulants # replaces both!
    }
    controls$eotd_anticoagulants[i] = this_daily$eotd_anticoagulants[1]
    # b) type
    this_daily = filter(daily, # get daily data
                        pin == controls$pin[i], # same patient
                        date_daily <= controls$date_daily[i], # dates prior to stroke
                        !is.na(eotd_anticoagulants_type)) %>% # exclude missing
      summarise(eotd_anticoagulants_type = Modes(eotd_anticoagulants_type)) # find mode
    # if equal modes take most recent 
    if(nrow(this_daily) > 1){
      this_daily_new = filter(daily, 
                              pin == controls$pin[i], # same patient
                              date_daily <= controls$date_daily[i],
                              eotd_anticoagulants_type %in% this_daily$eotd_anticoagulants_type) %>%
        arrange(date_daily) %>%
        slice(n()) # latest result
      this_daily$eotd_anticoagulants_type = this_daily_new$eotd_anticoagulants_type # replaces both!
    }
    controls$eotd_anticoagulants_type[i] = this_daily$eotd_anticoagulants_type[1]
  }
  
  ## now match controls to case by age and date:
  matched = NULL
  for (k in 1:nrow(cases)){
    this_case = cases[k,] %>%
      mutate(case=1)
    random_controls = filter(controls, 
                             days_c == this_case$days_c, # same (grouped) day in ICU
                             age >= (this_case$age - 2), # age within 2 years
                             age <= (this_case$age + 2)) %>%
      sample_n(number_of_controls) %>% # randomly sample
      mutate(case = 0)
    comb = bind_rows(this_case, random_controls) %>%
      mutate(pair = k)
    matched = bind_rows(matched, comb)
  }
  #
  table = group_by(matched, case, eotd_anticoagulants) %>%
    tally() %>%
    group_by(case) %>%
    mutate(percent = round(prop.table(n)*100),
           cell = paste(n, ' (', percent, ')', sep='')) %>%
    select(-n, -percent) %>%
    ungroup() %>%
    spread(case, cell) %>%
    mutate(eotd_anticoagulants = case_when(
      eotd_anticoagulants == 0 ~ 'No',
      eotd_anticoagulants == 1 ~ 'Yes',
      is.na(eotd_anticoagulants) ~ 'Missing'
    ))
  #
  # table header
  table_header <- data.frame(
    col_keys = c('eotd_anticoagulants','0','1'),
    h1 = c('', 'Stroke', 'Stroke'),
    h2 = c('Anticoagulants','No','Yes'),
    stringsAsFactors = FALSE )
  #
  ftab = flextable(table) %>% 
    set_header_df(mapping = table_header, key = "col_keys" ) %>% # add header
    merge_h(i=1, part = "header") %>% # merge column headers
    theme_box() %>%
    autofit()
  # chi-squared test
  tab = with(matched, table(case, eotd_anticoagulants))
  chisq = chisq.test(tab, simulate.p.value = TRUE)
  # return
  to_return = list()
  to_return$chisq = chisq
  to_return$table = ftab
  to_return$matched = matched
  return(to_return)
} # end of function


## function to run case-control analysis
# uses daily result just before stroke date
case_control_two = function(indata,
                            indaily,
                            scale = 1, 
                            s_type, 
                            rdigits=1, # rounding digits for table
                            log_exp = FALSE, # log exposure (if there is skew)
                            exposure = '', # exposure variable
                        number_of_controls = 5
){
  set.seed(212121) # for random control selection
  # simple rename for exposure variable
  index = names(indaily) == exposure
  names(indaily)[index] = 'exposure'
  # split into cases and controls by stroke
  cases = filter(indata, 
                 !is.na(complication_stroke_date), # has a stroke ...
                 stroke_group %in% s_type) %>% # ... of this type
    select(pin, age, complication_stroke_date, date_icu) %>%
    mutate(days = as.numeric(complication_stroke_date - date_icu)) %>%
    filter(!is.na(days),
           days > 0) %>% # stroke must be after ICU admission
    mutate(days_c = case_when( # combine days as they get sparse
      days > 20 & days<=30 ~ 25,
      days > 30 & days<=40 ~ 35,
      days > 40 & days<=50 ~ 45,
      days > 50 ~ 55,
      TRUE ~ days
    ))
  controls = filter(indata, 
                    !pin %in% unique(cases$pin), # not in the stroke case data
                    !is.na(date_icu)) %>% # with ICU data
    select(pin, age, date_icu)
  
  ## merge cases with daily data; use day before stroke
  cases = mutate(cases, exp = NA) # start with blanks
  for (i in 1:nrow(cases)){
    this_daily = filter(indaily, # get daily data
                        pin == cases$pin[i], # same patient
                        date_daily < cases$complication_stroke_date[i], # dates prior to stroke (less than)
                        !is.na(exposure)) %>% # exclude missing exposure
        arrange(date_daily) %>%
        slice(n()) # latest result
    # update exposure if there's any data
    if(nrow(this_daily) > 0 ) {cases$exp[i] = this_daily$exposure}
  }
  
  ## merge controls with daily data (take all days)
  controls = left_join(controls, indaily, by=c('pin'='pin')) %>%
    mutate(days = as.numeric(date_daily - date_icu)) %>%
    filter(!is.na(days),
           !is.na(exposure), # must have exposure
           days > 0) %>%
    mutate(days_c = case_when( # combine days as they get sparse
      days > 20 & days<=30 ~ 25,
      days > 30 & days<=40 ~ 35,
      days > 40 & days<=50 ~ 45,
      days > 50 ~ 55,
      TRUE ~ days
    )) %>%
    rename('exp' = 'exposure')
  
  ## now match controls to case by age and date:
  matched = NULL
  for (k in 1:nrow(cases)){
    this_case = cases[k,] %>%
      mutate(case=1)
    if(is.na(this_case$exp) == TRUE){next} # skip to next if case has no exposure
    random_controls = filter(controls, 
                             days_c == this_case$days_c, # same (grouped) day in ICU
                             age >= (this_case$age - 2), # age within 2 years
                             age <= (this_case$age + 2)) %>%
      mutate(case = 0)
    if(nrow(random_controls) > number_of_controls){
      random_controls = sample_n(random_controls, number_of_controls) # randomly sample
    }
    comb = bind_rows(this_case, random_controls) %>%
      mutate(pair = k) %>%
      select(pair, pin, age, days, case, exp)
    matched = bind_rows(matched, comb)
  }
  # summary table before log transform
  table = group_by(matched, case) %>%
    summarise(n = n(), median=median(exp), q1=quantile(exp, 0.25), q3=quantile(exp, 0.75))
  # round
  table = mutate(table, median = round(median, rdigits),
         #mean = round(mean, 1),
         q1 = round(q1, rdigits),
         q3 = round(q3, rdigits),
         case = ifelse(case==0, 'Control', 'Case'))
  #
  ftab = flextable(table) %>% 
    theme_box() %>%
    autofit()
  
  ## conditional logistic regression
  # log or not (base 2)
  if(log_exp==TRUE){matched = mutate(matched, exp = log2(exp))}
  # scale, optional
  if(scale !=1){
    matched = mutate(matched, exp = exp / scale) # divide by scale to give per scale increase
  }
  model = clogit(case ~ exp + strata(pair), data=matched)
  ests = tidy(model, conf.int = TRUE, exponentiate = TRUE) %>% # 
    mutate(stype = s_type, exposure=exposure)
  # return
  to_return = list()
  to_return$model = ests
  to_return$table = ftab
  to_return$matched = matched
  return(to_return)
} # end of function

## function to run case-control analysis
# uses maximum exposure from daily data for dates from admission to stroke
case_control_three = function(indata,
                            indaily,
                            s_type, 
                            scale = 1,
                            rdigits=1, # rounding digits for table
                            log_exp = FALSE, # log exposure (if there is skew)
                            exposure = '', # exposure variable
                            number_of_controls = 5
){
  set.seed(212121) # for random control selection
  # simple rename for exposure variable
  index = names(indaily) == exposure
  names(indaily)[index] = 'exposure'
  # split into cases and controls by stroke
  cases = filter(indata, 
                 !is.na(complication_stroke_date), # has a stroke ...
                 stroke_group %in% s_type) %>% # ... of this type
    select(pin, age, complication_stroke_date, date_icu) %>%
    mutate(days = as.numeric(complication_stroke_date - date_icu)) %>%
    filter(!is.na(days),
           days > 0) %>% # stroke must be after ICU admission
    mutate(days_c = case_when( # combine days as they get sparse
      days > 20 & days<=30 ~ 25,
      days > 30 & days<=40 ~ 35,
      days > 40 & days<=50 ~ 45,
      days > 50 ~ 55,
      TRUE ~ days
    ))
  controls = filter(indata, 
                    !pin %in% unique(cases$pin), # not in the stroke case data
                    !is.na(date_icu)) %>% # with ICU data
    select(pin, age, date_icu)
  
  ## merge cases with daily data; use day before stroke
  cases = mutate(cases, exp = NA) # start with blanks
  for (i in 1:nrow(cases)){
    this_daily = filter(indaily, # get daily data
                        pin == cases$pin[i], # same patient
                        date_daily < cases$complication_stroke_date[i], # dates prior to stroke (less than)
                        !is.na(exposure)) # exclude missing exposure
    # update exposure if there's any data
    if(nrow(this_daily) > 0 ) {cases$exp[i] = max(this_daily$exp)} # maximum during this time
  }
  
  ## prepare controls
  controls = left_join(controls, indaily, by=c('pin'='pin')) %>%
    mutate(days = as.numeric(date_daily - date_icu)) %>%
    filter(!is.na(days),
           !is.na(exposure), # must have exposure
           days > 0)
  
  ## go through each control day! controls create multiple exposures
  # takes a long while
  potential_controls = NULL
  for (this_control in unique(controls$pin)) {
    this_control_data = filter(controls, pin == this_control)
    for (d in 1:max(this_control_data$days)){ # loop through all days
      this_data = filter(this_control_data, days<=d) # up to day d
      if(nrow(this_data)>0){
        this_data = group_by(this_data, pin, age) %>% # keep these variables
          summarise(exposure = max(exposure)) %>%
          mutate(days = d) # add day
        potential_controls = bind_rows(potential_controls, this_data)
      }
    }
  }
  potential_controls = ungroup(potential_controls) %>%
    rename('exp' = 'exposure') %>% # for below
    mutate(days = as.numeric(days)) # needed for case_when below
  
  # now truncate days
  potential_controls = mutate(potential_controls, 
                              days_c = case_when( # combine days as they get sparse
    days > 20 & days<=30 ~ 25,
    days > 30 & days<=40 ~ 35,
    days > 40 & days<=50 ~ 45,
    days > 50 ~ 55,
    TRUE ~ days
  )) %>%
    select(-days) %>%
    unique() # remove duplicates
  
  ## now match controls to case by age and date:
  matched = NULL
  for (k in 1:nrow(cases)){
    this_case = cases[k,] %>%
      mutate(case=1)
    if(is.na(this_case$exp) == TRUE){next} # skip to next if case has no exposure
    random_controls = filter(potential_controls, 
                             days_c == this_case$days_c, # same (grouped) day in ICU
                             age >= (this_case$age - 2), # age within 2 years
                             age <= (this_case$age + 2)) %>%
      mutate(case = 0)
    if(nrow(random_controls) > number_of_controls){
      random_controls = sample_n(random_controls, number_of_controls) # randomly sample
    }
    comb = bind_rows(this_case, random_controls) %>%
      mutate(pair = k) %>%
      select(pair, pin, age, days_c, case, exp)
    matched = bind_rows(matched, comb)
  }
  # summary table before log transform
  table = group_by(matched, case) %>%
    summarise(n = n(), 
              #mean=mean(exp), 
              median = median(exp), 
              q1 = quantile(exp, 0.25), 
              q3 = quantile(exp, 0.75)) %>%
    mutate(median = round(median, rdigits),
           q1 = round(q1, rdigits),
           q3 = round(q3, rdigits),
           case = ifelse(case==0, 'Control', 'Case'))
  #
  ftab = flextable(table) %>% 
    theme_box() %>%
    autofit()
  # log or not (base 2)
  if(log_exp==TRUE){
    matched = mutate(matched, exp = log2(exp))
  }
  ## conditional logistic regression
  # scale, optional
  if(scale !=1){
    matched = mutate(matched, exp = exp / scale) # divide by scale to give per scale increase
  }
  model = clogit(case ~ exp + strata(pair), data=matched)
  ests = tidy(model, conf.int = TRUE, exponentiate = TRUE) %>% # 
    mutate(stype = s_type, exposure=exposure) 
  # return
  to_return = list()
  to_return$model = ests
  to_return$table = ftab
  to_return$matched = matched
  return(to_return)
} # end of function


## function to run case-control analysis
# uses maximum exposure from daily data for dates from admission to stroke
# also matches on daily ECMO (yes/no)
case_control_four = function(indata,
                              indaily,
                              s_type, 
                              scale = 1,
                              rdigits=1, # rounding digits for table
                              log_exp = FALSE, # log exposure (if there is skew)
                              exposure = '', # exposure variable
                              number_of_controls = 5
){
  set.seed(212121) # for random control selection
  # simple rename for exposure variable
  index = names(indaily) == exposure
  names(indaily)[index] = 'exposure'
  # split into cases and controls by stroke
  cases = filter(indata, 
                 !is.na(complication_stroke_date), # has a stroke ...
                 stroke_group %in% s_type) %>% # ... of this type
    select(pin, age, complication_stroke_date, date_icu) %>%
    mutate(days = as.numeric(complication_stroke_date - date_icu)) %>%
    filter(!is.na(days),
           days > 0) %>% # stroke must be after ICU admission
    mutate(days_c = case_when( # combine days as they get sparse
      days > 20 & days<=30 ~ 25,
      days > 30 & days<=40 ~ 35,
      days > 40 & days<=50 ~ 45,
      days > 50 ~ 55,
      TRUE ~ days
    ))
  controls = filter(indata, 
                    !pin %in% unique(cases$pin), # not in the stroke case data
                    !is.na(date_icu)) %>% # with ICU data
    select(pin, age, date_icu)
  
  ## merge cases with daily data; use day before stroke
  cases = mutate(cases, exp = NA, # start with blanks
                 last_ecmo = NA)
  for (i in 1:nrow(cases)){
    this_daily = filter(indaily, # get daily data
                        pin == cases$pin[i], # same patient
                        date_daily < cases$complication_stroke_date[i], # dates prior to stroke (less than)
                        !is.na(exposure))  # exclude missing exposure
      # update exposure if there's any data
      if(nrow(this_daily) > 0 ) {
        cases$exp[i] = max(this_daily$exp) # maximum during this time
        # get latest ECMO stats
        last_ecmo = filter(this_daily, !is.na(ecmo)) %>% # exclude any missing ECMO days
          arrange(desc(date_daily)) %>%  # get the latest ECMO status
          slice(1) %>%
          pull(ecmo)
        cases$last_ecmo[i] = last_ecmo
      }
  }
  
  ## prepare controls
  controls = left_join(controls, indaily, by=c('pin'='pin')) %>%
    mutate(days = as.numeric(date_daily - date_icu)) %>%
    filter(!is.na(days), # must have days for matching
           !is.na(exposure), # must have exposure
           days > 0)
  
  ## go through each control day! controls create multiple exposures
  # takes a long while
  potential_controls = NULL
  for (this_control in unique(controls$pin)) {
    this_control_data = filter(controls, pin == this_control)
    for (d in 1:max(this_control_data$days)){ # loop through all days
      this_data = filter(this_control_data, days<=d) # up to day d
      if(nrow(this_data)>0){
        this_data = group_by(this_data, pin, age, ecmo) %>% # keep these variables
          summarise(exposure = max(exposure)) %>%
          mutate(days = d) # add day
        potential_controls = bind_rows(potential_controls, this_data)
      }
    }
  }
  potential_controls = ungroup(potential_controls) %>%
    rename('exp' = 'exposure') %>% # for below
    mutate(days = as.numeric(days)) # needed for case_when below
  
  # now truncate days
  potential_controls = mutate(potential_controls, 
                              days_c = case_when( # combine days as they get sparse
                                days > 20 & days<=30 ~ 25,
                                days > 30 & days<=40 ~ 35,
                                days > 40 & days<=50 ~ 45,
                                days > 50 ~ 55,
                                TRUE ~ days
                              )) %>%
    select(-days) %>%
    unique() # remove duplicates
  
  ## set up missing ECMO as separate group
  cases = mutate(cases, 
                 last_ecmo = ifelse(is.na(last_ecmo)==T, 3, last_ecmo))
  potential_controls = mutate(potential_controls, 
                 last_ecmo = ifelse(is.na(ecmo)==T, 3, ecmo))
  
  ## now match controls to case by age, date and ECMO:
  matched = NULL
  for (k in 1:nrow(cases)){
    this_case = cases[k,] %>%
      mutate(case=1, 
             ecmo=last_ecmo)
    if(is.na(this_case$exp) == TRUE){next} # skip to next if case has no exposure
    random_controls = filter(potential_controls, 
                             ecmo == this_case$last_ecmo, # same ECMO status (yes, no or missing)
                             days_c == this_case$days_c, # same (grouped) day in ICU
                             age >= (this_case$age - 2), # age within 2 years
                             age <= (this_case$age + 2)) %>%
      mutate(case = 0)
    if(nrow(random_controls) > number_of_controls){
      random_controls = sample_n(random_controls, number_of_controls) # randomly sample
    }
    comb = bind_rows(this_case, random_controls) %>%
      mutate(pair = k) %>%
      select(pair, pin, age, ecmo, days_c, case, exp)
    matched = bind_rows(matched, comb)
  }
  # summary table before log or not
  table = group_by(matched, case) %>%
    summarise(n = n(), 
              #mean=mean(exp), 
              median=median(exp), q1=quantile(exp, 0.25), q3=quantile(exp, 0.75)) %>%
    mutate(median = round(median, rdigits),
           q1 = round(q1, rdigits),
           q3 = round(q3, rdigits),
           case = ifelse(case==0, 'Control', 'Case'))
  #
  ftab = flextable(table) %>% 
    theme_box() %>%
    autofit()
  ## conditional logistic regression
  # log or not (base 2)
  if(log_exp==TRUE){
    matched = mutate(matched, exp = log2(exp))
  }
    # scale, optional
  if(scale !=1){
    matched = mutate(matched, exp = exp / scale) # divide by scale to give per scale increase
  }
  model = clogit(case ~ exp + strata(pair), data=matched)
  ests = tidy(model, conf.int = TRUE, exponentiate = TRUE) %>% # 
    mutate(stype = s_type, exposure=exposure) 
  # return
  to_return = list()
  to_return$model = ests
  to_return$table = ftab
  to_return$matched = matched
  return(to_return)
} # end of function


## function to run case-control analysis
# uses daily result just before stroke date
# also matches on daily ECMO (yes/no)
case_control_five = function(indata,
                             indaily,
                             s_type, 
                             scale = 1,
                             rdigits=1, # rounding digits for table
                             log_exp = FALSE, # log exposure (if there is skew)
                             exposure = '', # exposure variable
                             number_of_controls = 5
){
  set.seed(212121) # for random control selection
  # simple rename for exposure variable
  index = names(indaily) == exposure
  names(indaily)[index] = 'exposure'
  # split into cases and controls by stroke
  cases = filter(indata, 
                 !is.na(complication_stroke_date), # has a stroke ...
                 stroke_group %in% s_type) %>% # ... of this type
    select(pin, age, complication_stroke_date, date_icu) %>%
    mutate(days = as.numeric(complication_stroke_date - date_icu)) %>%
    filter(!is.na(days),
           days > 0) %>% # stroke must be after ICU admission
    mutate(days_c = case_when( # combine days as they get sparse
      days > 20 & days<=30 ~ 25,
      days > 30 & days<=40 ~ 35,
      days > 40 & days<=50 ~ 45,
      days > 50 ~ 55,
      TRUE ~ days
    ))
  controls = filter(indata, 
                    !pin %in% unique(cases$pin), # not in the stroke case data
                    !is.na(date_icu)) %>% # with ICU data
    select(pin, age, date_icu)
  
  ## merge cases with daily data; use day before stroke
  cases = mutate(cases, exp = NA, # start with blanks
                 last_ecmo = NA)
  for (i in 1:nrow(cases)){
    this_daily = filter(indaily, # get daily data
                        pin == cases$pin[i], # same patient
                        date_daily < cases$complication_stroke_date[i], # dates prior to stroke (less than)
                        !is.na(exposure))  # exclude missing exposure
    # update exposure if there's any data
    if(nrow(this_daily) > 0 ) {
      # get latest ECMO stats
      just_prior_stroke = filter(this_daily, !is.na(ecmo)) %>% # exclude any missing ECMO days
        arrange(desc(date_daily)) %>%  # get the date closest to the stroke
        slice(1) 
      cases$exp[i] = just_prior_stroke$exp # 
      cases$last_ecmo[i] = just_prior_stroke$ecmo
    }
  }
  
  ## prepare controls
  controls = left_join(controls, indaily, by=c('pin'='pin')) %>%
    mutate(days = as.numeric(date_daily - date_icu)) %>%
    filter(!is.na(days), # must have days for matching
           !is.na(exposure), # must have exposure
           days > 0)
  
  ## go through each control day! controls create multiple exposures
  # takes a long while
  potential_controls = NULL
  for (this_control in unique(controls$pin)) {
    this_control_data = filter(controls, pin == this_control)
    for (d in 1:max(this_control_data$days)){ # loop through all days
      this_data = filter(this_control_data, days<=d) # up to day d
      if(nrow(this_data)>0){
        this_data = group_by(this_data, pin, age, ecmo) %>% # keep these variables
          summarise(exposure = max(exposure)) %>%
          mutate(days = d) # add day
        potential_controls = bind_rows(potential_controls, this_data)
      }
    }
  }
  potential_controls = ungroup(potential_controls) %>%
    rename('exp' = 'exposure') %>% # for below
    mutate(days = as.numeric(days)) # needed for case_when below
  
  # now truncate days
  potential_controls = mutate(potential_controls, 
                              days_c = case_when( # combine days as they get sparse
                                days > 20 & days<=30 ~ 25,
                                days > 30 & days<=40 ~ 35,
                                days > 40 & days<=50 ~ 45,
                                days > 50 ~ 55,
                                TRUE ~ days
                              )) %>%
    select(-days) %>%
    unique() # remove duplicates
  
  ## set up missing ECMO as separate group
  cases = mutate(cases, 
                 last_ecmo = ifelse(is.na(last_ecmo)==T, 3, last_ecmo))
  potential_controls = mutate(potential_controls, 
                              last_ecmo = ifelse(is.na(ecmo)==T, 3, ecmo))
  
  ## now match controls to case by age, date and ECMO:
  matched = NULL
  for (k in 1:nrow(cases)){
    this_case = cases[k,] %>%
      mutate(case=1, 
             ecmo=last_ecmo)
    if(is.na(this_case$exp) == TRUE){next} # skip to next if case has no exposure
    random_controls = filter(potential_controls, 
                             ecmo == this_case$last_ecmo, # same ECMO status (yes, no or missing)
                             days_c == this_case$days_c, # same (grouped) day in ICU
                             age >= (this_case$age - 2), # age within 2 years
                             age <= (this_case$age + 2)) %>%
      mutate(case = 0)
    if(nrow(random_controls) > number_of_controls){
      random_controls = sample_n(random_controls, number_of_controls) # randomly sample
    }
    comb = bind_rows(this_case, random_controls) %>%
      mutate(pair = k) %>%
      select(pair, pin, age, ecmo, days_c, case, exp)
    matched = bind_rows(matched, comb)
  }
  # summary table before log transform
  table = group_by(matched, case) %>%
    summarise(n = n(), 
              #mean=mean(exp), 
              median=median(exp), q1=quantile(exp, 0.25), q3=quantile(exp, 0.75)) %>%
    mutate(median = round(median, rdigits),
           q1 = round(q1, rdigits),
           q3 = round(q3, rdigits),
           case = ifelse(case==0, 'Control', 'Case'))
  #
  ftab = flextable(table) %>% 
    theme_box() %>%
    autofit()
  ## conditional logistic regression
  # log or not (base 2)
  if(log_exp==TRUE){matched = mutate(matched, exp = log2(exp))}
  # scale, optional
  if(scale !=1){
    matched = mutate(matched, exp = exp / scale) # divide by scale to give per scale increase
  }
  model = clogit(case ~ exp + strata(pair), data=matched)
  ests = tidy(model, conf.int = TRUE, exponentiate = TRUE) %>% # 
    mutate(stype = s_type, exposure=exposure) 
  # return
  to_return = list()
  to_return$model = ests
  to_return$table = ftab
  to_return$matched = matched
  return(to_return)
} # end of function



## fractional polynomial functions for lmer/glmer (first and second order)
transform_fp_terms <- function(x.in,p1=1,p2=NULL,shift=0){
  x.in <- x.in+shift  # apply shift to all x values to prevent log(0) (shift=0 by default)
  m = ifelse(is.null(p2),1,2)
  out = array(0,c(length(x.in),m))
  colnames(out) = paste('fp',1:m,sep='.')
  #first-order
  out[,1] = (p1!=0)*(x.in^p1)+(p1==0)*log(x.in) 
  #second-order; multiply by log(x) if p1=p2 i.e. repeated powers
  if (m==2){
      out[,2] = (p2!=0)*(x.in^p2)+(p2==0)*log(x.in)
    if (p1==p2){
      out[,2] = out[,2]*log(x.in)
    }
  }
  out = out 
  return(out)
}





