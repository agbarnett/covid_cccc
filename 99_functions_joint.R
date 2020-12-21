# 99_functions_joint.R
# functions for joint models examining time to stroke/death/discharge
# copied from Nicole's code in joint models folder
# Nov 2000

joint_survival = function(indata.surv, # survival data
                          indata.long, # longitudinal data
                          event_name='stroke', # event, e.g. discharge,death,stroke
                          intercept.survival = FALSE, # random site-specific frailty
                          long_var, # variable to model longitudinally, e.g., pa_co2
                          scale = 1, # scale for longitudinal variable
                          log_long = FALSE, # log-transform (base 2) longitudinal variable
                          censor_day, 
                          outliers = NA) # cut outliers at low end (used for platelet count)
  {

  ## section 0: rename longitudinal and event variables ##
  index = names(indata.long) == long_var
  names(indata.long)[index] = 'long_var'
  # remove missing
  indata.long = filter(indata.long, 
                       !is.na(long_var))
  # cut 'deranged values' (optional)
  if(is.na(outliers)== FALSE){
    indata.long = filter(indata.long, 
                         long_var > outliers)
  }

## section 1: set up the data for the joint models
  
inla.data.surv =  select(indata.surv, pin, date_icu, from, to, entry, exit, site_name, sex, age_c, bmi_c) %>% # 
  mutate(
    clock_reset = exit-entry,
    ecmo_start = as.numeric(from=='ecmo_start'),
    ecmo_end = as.numeric(from=='ecmo_end'),
    sex = as.numeric(sex=='Male'),
    date_start = date_icu+entry, # for merging with longitudinal data
    date_end = date_icu+exit,
    event = as.numeric(to==event_name)) %>% # event is having a stroke/discharge/death
  drop_na() 
# censor long times
index = inla.data.surv$clock_reset > censor_day
inla.data.surv$clock_reset[index] = censor_day
inla.data.surv$event[index] = 0
#
inla.data.long = select(indata.long, pin, date_daily, long_var) %>% 
  mutate(
    long_var_plot = long_var, # keep original scale for plot...
    long_var = long_var / scale, #...but rescale for regression 
    date_daily = as.Date(date_daily, format='%Y-%m-%d')) %>%
  drop_na() #non-missing longitudinal values only

# common pins in longitudinal and survival data
u = intersect(inla.data.surv[['pin']],inla.data.long[['pin']])
inla.data.surv = inla.data.surv %>% filter(pin %in% u)
inla.data.long = inla.data.long %>% filter(pin %in% u)

# now create version of longitudinal data where the dates match those from each entry->exit transition in the survival data
inla.data.surv = mutate(inla.data.surv, row = 1:n())
new_long_data = NULL
for (k in 1:nrow(inla.data.surv)){
  this_long = filter(inla.data.long, 
                     pin == inla.data.surv$pin[k], # same patient
                     date_daily >= inla.data.surv$date_start[k], # within date range  
                     date_daily < inla.data.surv$date_end[k]) %>%
    mutate (row =inla.data.surv$row[k] ) # for later matching
  new_long_data = bind_rows(new_long_data, this_long)
}

# common rows in longitudinal and survival data (need patients to be in both)
u = intersect(inla.data.surv[['row']],new_long_data[['row']])
inla.data.surv = inla.data.surv %>% filter(row %in% u)
new_long_data = new_long_data %>% filter(row %in% u)

# add site_name to long
sites = select(inla.data.surv, pin, site_name) %>%
  unique()
new_long_data = left_join(new_long_data, sites, by='pin')

#reorder pins, sites to tidy up.
#define subject based on unique rows in survival dataset
inla.data.surv = inla.data.surv %>%  
  mutate(site = match(site_name, unique(site_name)),
         subject = match(row, unique(row)))

new_long_data = new_long_data %>% 
  mutate(subject = match(row, unique(inla.data.surv[['row']])),
         site=match(site_name, unique(inla.data.surv[['site_name']])))

# add age, sex to new inla.data.long
small =  select(inla.data.surv, subject, age_c, sex, date_icu)
new_long_data = left_join(new_long_data, small, by='subject') %>%
  mutate(day = as.numeric(date_daily - date_icu)) %>%
  filter(day >= 0) # remove obvious errors

#First, we will create some auxiliary variables that we will use later. 
#Essentially, these variables will help us fill the empty spaces of the variables used in the joint models with NA or zero.
#Number of observations in each dataset
n.l <- nrow(new_long_data)
n.s <- nrow(inla.data.surv)

#Vector of NA's
NAs.l <- rep(NA, n.l)
NAs.s <- rep(NA, n.s)
#Vector of zeros
zeros.l <- rep(0, n.l)
zeros.s <- rep(0, n.s)

#The longitudinal and survival responses, now called Y.long and Y.surv, 
#need to be expanded with NAs so that both variables have the same length and then put together in a list:
#Long. and survival responses
if(log_long == TRUE){Y.long <- c(log2(new_long_data$long_var), NAs.s)} # logged 
if(log_long == FALSE){Y.long <- c(new_long_data$long_var, NAs.s)} # 
Y.surv <- inla.surv(time = c(NAs.l, inla.data.surv$clock_reset),
                    event = c(NAs.l, inla.data.surv$event))

Y.joint <- list(Y.long, Y.surv)

#Expand covariates in the same way 
#Covariates
covariates <- data.frame(
  #Intercepts (as factor with 2 levels)
  inter = as.factor(c(rep("b.l", n.l), rep("b.s", n.s))),
  # Time -linear
  time.l = c(new_long_data$day/10, NAs.s), # scale by 10 days
  # Time - non-linear
  time2.l = c(sqrt(new_long_data$day/10), NAs.s), # scale by 10 days
  #age per 5 years 
  ##longitudinal
  age10.l = c(new_long_data$age_c, NAs.s),
  ##survival
  age10.s = c(NAs.l, inla.data.surv$age_c),
  #sex
  ##longitudinal
  sex.l = c(new_long_data$sex, NAs.s),
  ##survival
  sex.s = c(NAs.l, inla.data.surv$sex),
  # variables for survival only
  ecmo_end.s = c(NAs.l, inla.data.surv$ecmo_end),  
  ecmo_start.s = c(NAs.l, inla.data.surv$ecmo_start)
)

#Indices for random effects
r.effects <- list(
  #Patient id (long.)
  id.l = c(new_long_data$subject, NAs.s),
  #Site id (long.)
  site.l = c(new_long_data$site, NAs.s),  
  #Patient id (surv.)
  id.s = c(NAs.l, inla.data.surv$subject),
  #Site id (surv.)
  site.s = c(NAs.l, inla.data.surv$site)
)

#create a new object, joint.data, with the covariates, the indices of the random effects and the response variables:
joint.data  <- c(covariates, r.effects)
joint.data$Y <- Y.joint

#Unique index for patient from 1 to n.patient
unique.id <- unique(inla.data.surv$subject)
n.id <- length(unique.id)
idx <- 1:n.id

#Unique indices for long. and survival data
idx.l <- match(new_long_data$subject, unique.id)
idx.s <- match(inla.data.surv$subject, unique.id)

## Then, indices are created to indicate the patients in the longitudinal and survival part of the model:
#Longitudinal random effects per subject
joint.data$idx.sh.l <- c(idx.l, NAs.s) #intercept
#joint.data$idx.sh.l2 <- c(idx.l, NAs.s) #slope (time)

#Survival random effects (subject-specific copied from longitudinal component)
joint.data$idx.sh.s <- c(NAs.l, idx.s)
#joint.data$idx.sh.s2 <- c(NAs.l, idx.s)

## Section 2: run model ## 
#fit long and surv models separately using joint model specification
#should give same answer as long, surv models when fitted using separate syntax
if(intercept.survival == FALSE){
  joint.inla <- inla(Y ~ -1 + 
                     # intercepts
                     inter +
                     # Longitudinal part: time, age10, sex as covariates, random intercept
                     time.l + time2.l + age10.l + sex.l + # fixed effects including non-linear time
                     f(idx.sh.l, model = "iid") + # random intercept
                     f(site.l, model = "iid", hyper = list(prec = list(param = c(0.001, 0.001))))+ #site random intercept
                     # Survival model
                     age10.s + sex.s + ecmo_start.s + ecmo_end.s +# fixed effects
                     f(idx.sh.s, copy = "idx.sh.l", fixed = FALSE) 
                   ,
                   data = joint.data, family = c("normal", "loglogistic.surv"),
                   verbose=FALSE)
  }
if(intercept.survival == TRUE){
  joint.inla <- inla(Y ~ -1 + 
                       # intercepts
                       inter +
                       # Longitudinal part: time, age10, sex as covariates, random intercept
                       time.l + time2.l + age10.l + sex.l + # fixed effects including non-linear time
                       f(idx.sh.l, model = "iid") + # random intercept
                       f(site.l, model = "iid", hyper = list(prec = list(param = c(0.001, 0.001))))+ #site random intercept
                       # Survival model
                       age10.s + sex.s + ecmo_start.s + ecmo_end.s +# fixed effects
                       f(site.s, model = "iid", hyper = list(prec = list(param = c(0.001, 0.001))))+ #site specific frailty ### ? - did not converge with random intercept as well
                       f(idx.sh.s, copy = "idx.sh.l", fixed = FALSE) 
                     ,
                     data = joint.data, family = c("normal", "loglogistic.surv"),
                     verbose=FALSE)
}
# shape is "alpha for loglogistic observations"
# summary(joint.inla)

## Section 3: nice table ## 
ests_fixed = select(joint.inla$summary.fixed, mean, '0.025quant',  '0.975quant') %>%
  mutate(term = row.names(joint.inla$summary.fixed))
ests_joint = select(joint.inla$summary.hyperpar, mean, '0.025quant',  '0.975quant') %>%
  mutate(term = row.names(joint.inla$summary.hyperpar)) %>%
  filter(str_detect(term, "^Beta"))
ests_full = bind_rows(ests_fixed, ests_joint)  %>%
  rename('lower' = '0.025quant',
         'upper' = '0.975quant') %>%
  mutate(model = ifelse(str_detect(term, pattern='\\.l$')==TRUE, 'Longitudinal',paste('Time to', event_name)), # extract model type
         long_var = long_var, # store what the longitudinal variable is
         term = str_remove_all(term, pattern='\\.l$|\\.s$'),
         term = ifelse(term=='interb', "Intercept", term),
         term = ifelse(term=='time', "Time (+10 days)", term),
         term = ifelse(term=='time2', "Square-root time", term),
         term = ifelse(term=='age10', "Age (+10 years)", term),
         term = ifelse(term=='sex', "Male", term),
         term = ifelse(term=='ecmo_end', "Post-ECMO", term),
         term = ifelse(term=='ecmo_start', "During ECMO", term),
         is.long = as.numeric(str_detect(term, '^Beta')), # indicator for longitudinal variable
         term = ifelse(str_detect(term, '^Beta'), long_var, term), # join
         mean = ifelse(str_detect(string=model, pattern='^Time to '), exp(mean), mean), # convert to multiplier
         lower = ifelse(str_detect(string=model, pattern='^Time to '), exp(lower), lower), # 
         upper = ifelse(str_detect(string=model, pattern='^Time to '), exp(upper), upper))
ests = mutate(ests_full, 
              mean = roundz(mean, 2), # rounding
              lower = roundz(lower, 2),
              upper = roundz(upper, 2),
              cell = paste(mean, ' (', lower, ' to ', upper, ')', sep='')) %>%
  arrange(model, term) %>%
  select(model, term, cell) %>%
  rename('Mean (95% CI)' = 'cell')
# remove intercept for survival 
index = ests$term=='Intercept'& str_detect(string=ests$model, pattern='^Time to ')
ests = ests[!index,]
# nice table
ftab = flextable(ests) %>%
  theme_box() %>%
  autofit()

## section 4: basic numbers ###
numbers = list()
numbers$N_event = sum(inla.data.surv$event) # number of events
numbers$N_surv = nrow(inla.data.surv)
numbers$N_long = nrow(new_long_data) # number of longitudinal observations
numbers$N_patients = length(unique(inla.data.surv$subject))
tab = table(new_long_data$pin) # results per patient
numbers$long_median = median(tab) # median results per patient
numbers$long_q1 = as.numeric(quantile(tab, 0.25))
numbers$long_q3 = as.numeric(quantile(tab, 0.75))

## Section 5: graph ## 
# add stroke outcome
colours = grey(runif(min=0.1, max=0.7, n = 1000))
small = select(inla.data.surv, subject, event, from)
to_plot = left_join(new_long_data, small, by='subject') %>%
  mutate(event_nice = ifelse(event==0, 'No stroke', 'Stroke'))
# thin out patients to plot
random_p = filter(to_plot, event==0) %>% # controls only
  sample_n(50, replace = FALSE) %>% # randomly sample 50 patients
  pull(subject)
#
to_plot = filter(to_plot, 
                 subject %in% random_p | event==1, # all cases or random controls
                 day <= censor_day) %>%
  mutate(from_nice = case_when(
    from == 'icu' ~ "ICU",
    from == 'ecmo_start' ~ "During ECMO",
    from == 'ecmo_end' ~ "Post-ECMO")
  )
#
lplot = ggplot(data=to_plot, aes(x=day, y=long_var_plot, colour=factor(subject), group=subject)) +
  geom_line()+
  scale_color_manual(NULL, values=colours)+
  stat_smooth(aes(group=1), method='lm', se=FALSE)+ # linear regression
  scale_y_log10()+
  theme(legend.position = 'none')+
  facet_grid(from_nice~event_nice)

## Section 6: return
to_return = list()
to_return$table = ftab # nice table of estimates ...
to_return$ests = ests_full # ... add raw estimates
to_return$graph = lplot
to_return$numbers = numbers
return(to_return)

} # end of function

