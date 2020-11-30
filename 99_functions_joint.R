# 99_functions_joint.R
# functions for joint models
# copied from Nicole's code in joint models folder
# Nov 2000

joint_survival = function(indata.surv, # survival data
                          indata.long, # longitudinal data
                          long_var, # variable to model longitudinally, e.g., pa_co2
                          scale = 1, # scale for longitudinal variable
                          log_long = FALSE, # log-transform (base 2) longitudinal variable
                          censor_day)
  {

  ## section 0: rename longitudinal variable ##
  index = names(indata.long) == long_var
  names(indata.long)[index] = 'long_var'
  
## section 1: set up the data for the joint models
  
inla.data.surv =  select(indata.surv, pin, date_icu, from, to, entry, exit, site_name, sex, age_c, bmi_c) %>%
  mutate(
    clock_reset = exit-entry,
    ecmo_start = as.numeric(from=='ecmo_start'),
    ecmo_end = as.numeric(from=='ecmo_end'),
    sex = as.numeric(sex=='Male'),
    event = as.numeric(to=='stroke')) %>% # event is having a stroke
  drop_na()
# censor long times
index = inla.data.surv$clock_reset > censor_day
inla.data.surv$clock_reset[index] = censor_day
inla.data.surv$event[index] = 0
#
inla.data.long = select(indata.long, pin, date_daily, long_var) %>% 
  mutate(
    long_var_plot = long_var, # keep original scale for plot
    long_var = long_var / scale, # 
    date_daily = as.Date(date_daily, format='%Y-%m-%d')) %>%
  drop_na() #non-missing longitudinal values only

# common pins in longitudinal and survival data
u = intersect(inla.data.surv[['pin']],inla.data.long[['pin']])
inla.data.surv = inla.data.surv %>% filter(pin %in% u)
inla.data.long = inla.data.long %>% filter(pin %in% u)
# add site_name to long
sites = select(inla.data.surv, pin, site_name) %>%
  unique()
inla.data.long = left_join(inla.data.long, sites, by='pin')

#reorder pins, sites to tidy up.
#define subject based on unique pins in survival dataset
inla.data.surv = inla.data.surv %>%  
  mutate(site = match(site_name, unique(site_name)),
         subject = match(pin, unique(pin)))

inla.data.long = inla.data.long %>% 
  mutate(subject = match(pin,unique(inla.data.surv[['pin']])),
         site=match(site_name, unique(inla.data.surv[['site_name']])))

# add age, sex to inla.data.long
small =  select(inla.data.surv, subject, age_c, sex, date_icu)
inla.data.long = left_join(inla.data.long, small, by='subject') %>%
  mutate(day = as.numeric(date_daily - date_icu)) %>%
  filter(day >= 0) # remove obvious errors

#First, we will create some auxiliary variables that we will use later. 
#Essentially, these variables will help us fill the empty spaces of the variables used in the joint models with NA or zero.
#Number of observations in each dataset
n.l <- nrow(inla.data.long)
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
if(log_long == TRUE){Y.long <- c(log2(inla.data.long$long_var) - mean(log2(inla.data.long$long_var)), NAs.s)} # logged and centred
if(log_long == FALSE){Y.long <- c(inla.data.long$long_var - mean(inla.data.long$long_var), NAs.s)} # centred
Y.surv <- inla.surv(time = c(NAs.l, inla.data.surv$clock_reset),
                    event = c(NAs.l, inla.data.surv$event))

Y.joint <- list(Y.long, Y.surv)

#Expand covariates in the same way 
#Covariates
covariates <- data.frame(
  #Intercepts (as factor with 2 levels)
  inter = as.factor(c(rep("b.l", n.l), rep("b.s", n.s))),
  #Time
  time.l = c(inla.data.long$day/10, NAs.s), # scale by 10 days
  #age per 5 years 
  ##longitudinal
  age10.l = c(inla.data.long$age_c, NAs.s),
  ##survival
  age10.s = c(NAs.l, inla.data.surv$age_c),
  #sex
  ##longitudinal
  sex.l = c(inla.data.long$sex, NAs.s),
  ##survival
  sex.s = c(NAs.l, inla.data.surv$sex),
  # variables for survival only
  ecmo_end.s = c(NAs.l, inla.data.surv$ecmo_end),  
  ecmo_start.s = c(NAs.l, inla.data.surv$ecmo_start)
)

#Indices for random effects
r.effects <- list(
  #Patient id (long.)
  id.l = c(inla.data.long$subject, NAs.s),
  #Site id (long.)
  site.l = c(inla.data.long$site, NAs.s),  
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
idx.l <- match(inla.data.long$subject, unique.id)
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
joint.inla <- inla(Y ~ -1 + 
                     # intercepts
                     inter +
                     # Longitudinal part - time, age10, sex as covariates, random intercept
                     time.l + age10.l + sex.l + # fixed effects
                     f(idx.sh.l, model = "iid") + #random intercept
                     f(site.l, model = "iid", hyper = list(prec = list(param = c(0.001, 0.001))))+ #site random intercept
                     # Survival model
                     age10.s + sex.s + ecmo_start.s + ecmo_end.s +# fixed effects
                     #f(site.s, model = "iid", hyper = list(prec = list(param = c(0.001, 0.001))))+ #site specific frailty - did not converge with random intercept as well
                     f(idx.sh.s, copy = "idx.sh.l", fixed = FALSE) 
                   ,
                   data = joint.data, family = c("normal", "loglogistic.surv"),
                   verbose=FALSE
)
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
  mutate(model = ifelse(str_detect(term, pattern='\\.l$')==TRUE, 'Longitudinal','Survival'), # extract model type
         term = str_remove_all(term, pattern='\\.l$|\\.s$'),
         term = ifelse(term=='interb', "Intercept", term),
         term = ifelse(term=='time', "Time (+10 days)", term),
         term = ifelse(term=='age10', "Age (+10 years)", term),
         term = ifelse(term=='sex', "Male", term),
         term = ifelse(term=='ecmo_end', "Post-ECMO", term),
         term = ifelse(term=='ecmo_start', "During ECMO", term),
         term = ifelse(str_detect(term, '^Beta'), long_var, term), # join
         mean = ifelse(model=='Survival', exp(mean), mean), # convert to multiplier
         lower = ifelse(model=='Survival', exp(lower), lower), # 
         upper = ifelse(model=='Survival', exp(upper), upper),
         mean = roundz(mean, 2), # rounding
         lower = roundz(lower, 2),
         upper = roundz(upper, 2),
         cell = paste(mean, ' (', lower, ' to ', upper, ')', sep=''))
ests = select(ests_full, model, term, cell) %>%
  arrange(model, term) %>%
  rename('Mean (95% CI)' = 'cell')
# remove intercept for survival 
index = ests$term=='Intercept'& ests$model=='Survival'
ests = ests[!index,]
# nice table
ftab = flextable(ests) %>%
  theme_box() %>%
  autofit()

## section 4: basic numbers ###
numbers = list()
numbers$N_event = sum(inla.data.surv$event) # number of events
numbers$N_surv = nrow(inla.data.surv)
numbers$N_long = nrow(inla.data.long) # number of longitudinal observations
numbers$N_patients = length(unique(inla.data.surv$subject))

## Section 5: graph ## 
# add stroke outcome
colours = grey(runif(min=0.1, max=0.7, n = 1000))
small = select(inla.data.surv, subject, event)
to_plot = left_join(inla.data.long, small, by='subject') %>%
  mutate(event_nice = ifelse(event==0, 'No stroke', 'Stroke'))
# thin out patients to plot
random_p = filter(to_plot, event==0) %>% # controls only
  sample_n(40, replace = FALSE) %>%
  pull(subject)
#
to_plot = filter(to_plot, 
                 subject %in% random_p | event==1, # all cases or random controls
                 day <= censor_day)
#
lplot = ggplot(data=to_plot, aes(x=day, y=long_var_plot, colour=factor(subject), group=subject)) +
  geom_line()+
  scale_color_manual(NULL, values=colours)+
  stat_smooth(aes(group=1), method='lm', se=FALSE)+ # linear regression
  scale_y_log10()+
  theme(legend.position = 'none')+
  facet_wrap(~event_nice)

## Section 6: return
to_return = list()
to_return$table = ftab # nice table of estimates ...
to_return$ests = ests_full # ... add raw estimates
to_return$graph = lplot
to_return$numbers = numbers
return(to_return)

} # end of function

