# 0_read_data.R
# read the data for: i) the stroke paper, ii) the ecmo papers, iii) obesity paper
# Jan 2021
library(dplyr)
options(dplyr.summarise.inform=FALSE) # suppress annoying warning
library(janitor)
library(stringr)

# date of this data
#source = 'dummy_data' # use dummy data
source = 'data' # use real data
data_date ='_12_21' # _month_day of download from UQ server (stroke)
#data_date ='_01_25' # _month_day of download from UQ server (obesity) - see also save at bottom of file
pfile = paste(source, '/patients', data_date, '.csv', sep='')
dfile = paste(source, '/daily', data_date, '.csv', sep='')

### Section 1 - get the latest data ###

## b) read the patient data ##
patients = read.csv(pfile, stringsAsFactors = FALSE) %>%
  janitor::clean_names() %>%
  unique() %>% # safety net for duplicates
  select(pin, site_name, country, sex, ethnicity, bmi, age, outcome, diagnosis_coronavirus, diagnosis_coronavirus_type,
         starts_with("comorbidity"),
         'stroke_during_treatment','stroke_during_treatment_side',
         starts_with("stroke_during_treatment_type_"),
         'transfer_from_other_facility',
         "diagnosis_coronavirus",
         "healthcare_worker", "lab_worker",
         "pregnant",
         'height','weight',
         # ecmo
         'date_outcome','date_ecmo','date_ecmo_discontinued','any_ecmo','ecmo_type',
         'ecmo_drainage_cannula_insertion_site','ecmo_drainage_cannula_size',
         'ecmo_return_cannula_insertion_site',"ecmo_return_cannula_size",
         "ecmo_cardiac_arrest_before", "ecmo_worst_pa_o2_6hr_before",
         'ecmo_highest_fi_o2_6hr_before',
         "ecmo_worst_pa_co2_6hr_before", 
         'ecmo_worst_pa_o2_6hr_before',
         "ecmo_prone_before",
         "complication_ards", "complication_viral_pheumonitis",
         "complication_stroke", 
         "complication_stroke_date", 
        # "neurology_diag_stroke" , # all missing
        # "neurology_diag_stroke_date" ,
         "complication_heart_failure", # to compare missingness
         "complication_seizure", # to compare missingness
         'ecmo_continuous_rpt_before',
         # clinical:
         'apache_ii','sofa','admission_temperature','admission_heart_rate',
         'admission_respiratory_rate','admission_o2_saturation',
         'admis_worst_pa_o2_6hr',
        # other additions (jan 2021)
        'treatment_trachoestomy','cause_of_death','ecmo_in_elso_registry','ability_to_self_care','mech_vent_cardiac_arrest', 'complication_cardiac_arrest' ,
        'ecmo_neuromuscolar_blockage_before',
        'ecmo_vasoactive_drugs_before', 'mech_vent_vasoactive_drugs_before', 'treatment_inotropes_vasopressors',
         # dates:
         starts_with('date_')) %>% # slim down variables for now
  mutate(healthlab_worker = pmax(healthcare_worker, lab_worker, na.rm = TRUE), # if either NA then use non-NA information
         date_did = NA, # set up empty dates for below
         date_dis = NA,
         date_d = NA,
         date_a = NA,
         date_i = NA,
         date_f = NA,
         date_o = NA,
         date_mv = NA,
         date_mv_end = NA,
         date_e = NA,
         date_s = NA)
# had to loop because these dates would not convert
for(j in 1:nrow(patients)){
  patients$date_did[j] = ifelse(patients$date_hospital_discharge[j] !='', as.Date(patients$date_hospital_discharge[j]), NA)
  patients$date_dis[j] = ifelse(patients$date_discharge[j] !='', as.Date(patients$date_discharge[j]), NA)
  patients$date_d[j] = ifelse(patients$date_death[j] !='', as.Date(patients$date_death[j]), NA)
  patients$date_a[j] = ifelse(patients$date_admission[j] !='', as.Date(patients$date_admission[j]), NA)
  patients$date_i[j] = ifelse(patients$date_icu[j] !='', as.Date(patients$date_icu[j]), NA)
  patients$date_o[j] = ifelse(patients$date_outcome[j] !='', as.Date(patients$date_outcome[j]), NA)
  patients$date_f[j] = ifelse(patients$date_first_symptom[j] !='', as.Date(patients$date_first_symptom[j]), NA)
  patients$date_mv[j] = ifelse(patients$date_mechanical_ventilation[j] !='', as.Date(patients$date_mechanical_ventilation[j]), NA)
  patients$date_mv_end[j] = ifelse(patients$date_mech_vent_discontinued[j] !='', as.Date(patients$date_mech_vent_discontinued[j]), NA)
  patients$date_e[j] = ifelse(patients$date_ecmo[j] !='', as.Date(patients$date_ecmo[j]), NA)
  patients$date_ed[j] = ifelse(patients$date_ecmo_discontinued[j] !='', as.Date(patients$date_ecmo_discontinued[j]), NA)
  patients$date_s[j] = ifelse(patients$complication_stroke_date[j] !='', as.Date(patients$complication_stroke_date[j]), NA)
}

## combine duplicate comorbidity fields: comorbidity_chronic_hematologic_disease, comorbidity_severe_liver_disease 
patients = patients %>% mutate(
  comorbidity_chronic_hematologic_disease_comb = pmax(comorbidity_chronic_hematologic_disease,
              comorbidity_chronic_hematologic_disease_eot,na.rm=T),
  comorbidity_severe_liver_disease_comb = pmax(comorbidity_severe_liver_disease,
              comorbidity_severe_liver_disease_eot,na.rm=T)) %>%
  select(-comorbidity_chronic_hematologic_disease, -comorbidity_severe_liver_disease) %>%
  # rename new variables to old to make life easier:
  rename('comorbidity_chronic_hematologic_disease' = 'comorbidity_chronic_hematologic_disease_comb',
         'comorbidity_severe_liver_disease' = 'comorbidity_severe_liver_disease_comb') 


## make new stroke type variable into list (August 2020) - CRF instructions: "please select up to two (2) options"
yes_responses <- dplyr::select(patients, pin, starts_with('stroke_during_treatment_type_'))  %>%
  mutate(stroke_during_treatment_type_ischemic_stroke = ifelse(stroke_during_treatment_type_ischemic_stroke==1, 'Ischemic stroke', ''), # convert 0,1 numbers to characters
         stroke_during_treatment_type_intraparenchymal_haemorrhage = ifelse(stroke_during_treatment_type_intraparenchymal_haemorrhage==1, 'Intraparenchymal haemorrhage', ''),
         stroke_during_treatment_type_subarachnoid_haemorrhage = ifelse(stroke_during_treatment_type_subarachnoid_haemorrhage==1, 'Subarachnoid haemorrhage', ''),
         stroke_during_treatment_type_hypoxic_ischemic_brain_injury = ifelse(stroke_during_treatment_type_hypoxic_ischemic_brain_injury==1, 'Hypoxic ischemic brain injury', ''),
         stroke_during_treatment_type_cerebral_venous_sinus_thrombosis  = ifelse(stroke_during_treatment_type_cerebral_venous_sinus_thrombosis ==1, 'Cerebral venous sinus thrombosis', ''),
         stroke_during_treatment_type_other  = ifelse(stroke_during_treatment_type_other ==1, 'Other', ''),
         stroke_during_treatment_type_unknown  = ifelse(stroke_during_treatment_type_unknown ==1, 'Unknown', '')
  ) %>%
  as_tibble() %>%
  tidyr::gather(key='type', value='response', -pin) %>%
  filter(!is.na(response),
         !response =='') # remove missing and empty
# now make list from group responses
list_resp = group_by(yes_responses, pin) %>%
  summarise(listed = list(response)) %>%
  ungroup()
names(list_resp)[2] = 'stroke_type'
patients = select(patients, -starts_with('stroke_during_treatment_type_')) # remove from data
patients = left_join(patients, list_resp, by='pin')

# more work on dates
patients = mutate(patients,
                  date_did = as.Date(date_did, origin='1970-01-01'),
                  date_dis = as.Date(date_dis, origin='1970-01-01'),
                  date_a = as.Date(date_a, origin='1970-01-01'),
                  date_i = as.Date(date_i, origin='1970-01-01'),
                  date_o = as.Date(date_o, origin='1970-01-01'),
                  date_d = as.Date(date_d, origin='1970-01-01'),
                  date_f = as.Date(date_f, origin='1970-01-01'),
                  date_mv = as.Date(date_mv, origin='1970-01-01'),
                  date_mv_end = as.Date(date_mv_end, origin='1970-01-01'),
                  date_e = as.Date(date_e, origin='1970-01-01'),
                  date_ed = as.Date(date_ed, origin='1970-01-01'),
                  date_s = as.Date(date_s, origin='1970-01-01')) %>%
  select(-date_discharge, -date_hospital_discharge, -date_death, -date_outcome, -date_first_symptom, -date_icu, -date_admission, -date_mechanical_ventilation, -date_mech_vent_discontinued, -date_ecmo, -date_ecmo_discontinued, -complication_stroke_date) %>% # use previous date names
  rename('date_discharge' = 'date_dis', # ICU discharge
         'date_hospital_discharge' = 'date_did',
         'date_admission' = 'date_a',
         'date_icu' = 'date_i',
         'date_outcome' = 'date_o',
         'date_death' = 'date_d',
         'date_first_symptom' = 'date_f',
         'date_mechanical_ventilation' = 'date_mv',
         'date_mech_vent_discontinued' = 'date_mv_end',
         'date_ecmo' = 'date_e',
         'date_ecmo_discontinued' = 'date_ed',
         'complication_stroke_date' = 'date_s'
         )
# quick check
# select(patients, pin, date_admission, date_icu, date_discharge, date_death) %>% sample_n(size=10)

## tidying: i) centre continuous patient variables, ii) missing categories
patients = mutate(patients,
                  bmi5 = (bmi - 27)/5, # standardise BMI
                  age5 = (age - 60)/5, # standardise age
                  # other tidying:
                  sex = ifelse(sex=='', 'Missing', sex),
                  ethnicity  = ifelse(ethnicity =='', 'Missing', ethnicity ),
                  ethnicity  = ifelse(ethnicity =='notavail', 'Missing', ethnicity ))

## c) read the daily data ##
daily = read.csv(dfile, stringsAsFactors = FALSE) %>%
  janitor::clean_names() %>%
  unique() %>% # safety net for duplicates
  mutate(date_daily = as.Date(date_daily)) %>%
  select(pin, day, date_daily, in_icu, mean_arterial_pressure, compliance_respiratory_system,
         wbc, pa_co2, serum_creatinine, sodium, potassium, crp, fibrinogen, troponin_i, haemoglobin, il_6,
         glucose, neutrophil_count, lymphocyte_count,
         # 
         'glasgow_coma_score','avpu',
         'prone_positioning',
         'ecmo', platelet_count, 'p_h', "aptt", "aptr", "aptt_aptr", "inr", alt_sgpt, "ast_sgot", 'bilirubin','d_dimer',
         'blood_urea_nitrogen', 'hco3', 'pa_co2','pa_o2', 'pa_o2_fi_o2', 'fi_o2',
         'eotd_anticoagulants', "eotd_anticoagulants_type",
         "eotd_haemorrhagic_complication", # 4.52 Haemorrhagic Complication 1
         'eotd_haemorrhagic_complication_source', # 4.53 Source Of Haemorrhagic Complication 1
         'eotd_peep','eotd_tidal_volume','eotd_ventilatory_mode','eotd_fi_o2','eotd_respiratory_rate',
         'eotd_tidal_volume_ideal','eotd_peep','eotd_airway_plateau_pressure',
         'tracheostomy_inserted',
         #
         'neuromuscular_blocking_agents', 'eotd_neuromuscolar_blockage',
         'inotropic_support', 'eotd_vasoactive_drugs',
         contains('infection')) %>% # slim down variables
  filter(day >= 0) %>% # remove clear errors in negative days
  arrange(pin, day)

## last date
admin_censor_date = as.Date(paste('2020', data_date, sep=''), format='%Y_%m_%d') # using date of data update
# add maximum date for censoring
max_date = group_by(daily, pin) %>%
  summarise(last_date = max(date_daily)) %>% 
  # if no last date then use admin censoring date
  mutate(last_date = ifelse(is.na(last_date) == TRUE, admin_censor_date, last_date), # fix a few impossible dates
         last_date = as.Date(last_date, origin='1970-01-01')) %>%
  ungroup()
# add date for administrative censoring (people still in hospital)

# add last date back to patients:
patients = left_join(patients, max_date, by='pin') %>% # add last date to patients
  # make a date of transfer if the outcome is transfer, censor patients at this time (later)
  mutate(date_transfer = ifelse(str_detect(string=outcome, pattern='^Transfer'), date_outcome, NA), # not for Hospitalization
         date_transfer = as.Date(date_transfer, origin='1970-01-01'))

## add dates of prone positioning and ecmo end
# a) dates ecmo 
dates_ecmo = filter(daily, ecmo==1) %>%
  select(pin, date_daily) %>%
  group_by(pin) %>%
  summarise(ecmo_start = min(date_daily),
            ecmo_end = max(date_daily)) %>%
  ungroup() 
# b) prone positioning
dates_prone = filter(daily, prone_positioning==1) %>%
  select(pin, date_daily) %>%
  group_by(pin) %>%
  summarise(prone_start = min(date_daily),
            prone_end = max(date_daily)) %>%
  ungroup() 
# add prone/ecmo dates to patients
patients = left_join(patients, dates_ecmo, by='pin')
patients = left_join(patients, dates_prone, by='pin')

# add stroke groups for those with a stroke, need to use loop because of lists
patients$stroke_group = ''
for (k in 1:nrow(patients)){
  types = unlist(patients$stroke_type[k])
  if(is.null(types)==FALSE){
    i = max(str_detect(pattern="Ischemic stroke|Hypoxic ischemic brain injury", string=types))
    h = max(str_detect(pattern="Intraparenchymal haemorrhage|Subarachnoid haemorrhage", string=types))
    o = max(str_detect(pattern="Other|Unknown", string=types))
    if(o == 1){patients$stroke_group[k] = 'Other/Unknown'}
    if(h == 1){patients$stroke_group[k] = 'Hemorrhage'}
    if(i == 1){patients$stroke_group[k] = 'Ischemic'} # if more than one pick this one
#    cat(k, i, h, o, patients$stroke_group[k], '\n', sep=' ')
  }
}
# other fixes
patients = mutate(patients,
          # change missing to `no stroke`
          stroke_group = ifelse(stroke_group=='', 'None', stroke_group),
          # if there is a date then make sure they are not in "None" group
          stroke_group = ifelse(!is.na(complication_stroke_date) & stroke_group=='None', 'Other/Unknown', stroke_group))


### Section 2: save ###

# dates for data
month = month.name[as.numeric(stringr::str_split(data_date, '_')[[1]][2])]
day = as.numeric(stringr::str_split(data_date, '_')[[1]][3])
data_date = list(month=month, day=day)

# save
#filename = 'data/COVID_Data_obese.RData' # obesity
filename = 'data/COVID_Data_stroke.RData' # stroke
if(source == 'dummy_data'){filename = 'data/COVID_Data_stroke_dummy.RData'}
save(daily, patients, data_date, file=filename)

# for finding variables:
# see also https://github.com/Samreay/COVID19/blob/master/definitions/oxford.yml
#names(patients)[grep('ecmo', names(patients))] 
#names(daily)[grep('_x', names(daily))] 


# check dates
select(patients, date_admission, date_icu, date_discharge, date_hospital_discharge, outcome, last_date) %>%
  filter(!is.na(date_discharge)) %>%
  sample_n(10)
