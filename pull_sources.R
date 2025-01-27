#remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"), INSTALL_opts = "--no-multiarch")

# install.packages("BiocManager")
# Load Packages -----------------------------------------------------------
pacman::p_load(tidyverse, pdftools, jsonlite, httr, lubridate, tabulapdf, MMWRweek)

# Functions ---------------------------------------------------------------
pull_api = function(endpoint, limit = "5000000") {
  return(read_csv(paste0(endpoint, "?$limit=", limit)))
}



# Pull NRVSS labs ---------------------------------------------------------


server = "https://nominatim.openstreetmap.org"

# Response 
participating_labs = extract_tables('https://www.cdc.gov/surveillance/nrevss/labs/list.pdf', 
                                    guess = F, 
                                    area = list(c(40, 18, 800, 800)), 
                                    output = "tibble") |> 
  map(.f = function(x){ 
    x %>% 
      janitor::remove_empty('cols') %>%
      set_names(c('state', 'institution'))}) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(url = paste0(server, "/search?q=", 
                      gsub(" ", "+", enc2utf8(institution), fixed = TRUE), 
                      "&format=jsonv2"), 
         try = list(fromJSON(url))) |> 
  mutate(out = list(tmaptools::geocode_OSM(institution)))

# Pull Hospitalization Data -----------------------------------------------
# Using old health data endpoint:
# hospitalizations_age <<- read_csv("https://healthdata.gov/resource/g62h-syeh.csv?$limit=500000") %>%
#   select(state, date,
#          starts_with("previous_day_admission_adult_covid_confirmed_"),
#          starts_with("previous_day_admission_pediatric_covid_confirmed_")) %>%
#   pivot_longer(starts_with("previous_day_admission")) %>%
#   filter(!str_detect(name, "_coverage$")) %>%
#   filter(!str_detect(name, "unknown$")) %>%
#   mutate(name = if_else(name == "previous_day_admission_adult_covid_confirmed_80+",
#                         "previous_day_admission_adult_covid_confirmed_80-99",
#                         name)) %>%
#   mutate(name = str_replace(name, "^.*confirmed_(.*)$", "\\1"), 
#          mmwr = MMWRweek::MMWRweek(as.Date(date)), 
#          mmwr_date = paste0(mmwr$MMWRyear, '-', mmwr$MMWRweek)) #%>%
  # gcmdata::splitAgeGroups("name", delim = "[^[:alnum:]]+") %>%
  # drop_na() %>%
  # mutate(value = as.integer(value),
  #        date = as.Date(date)) %>%
  # gcmdata::rescaleAgeGroupsValue() %>%
  # gcmdata::expandAgeGroups() %>%
  # gcmdata::tidyAgeGroupColumns() %>%
  # gcmdata::aggregateByAgeGroup(columnsToSummarize = "value",
  #                              ageGroupCrosswalk = gcmdata::fetchAgeGroupsCrosswalk("StandardCovidVaccineAgeGroups")) %>%
  # arrange(date) %>%
  # group_by(state, AgeGroupNames) %>%
  # mutate(value = zoo::rollmean(c(rep(NA, 3), value, rep(NA, 3)), 7)) %>%
  # ungroup %>%
  # rename(age_group = AgeGroupNames, data = value)

# Using new data.cdc.gov endpoint:
hospitalizations_age <<- read.csv("https://data.cdc.gov/resource/aemt-mg7g.csv?$limit=500000") %>% 
  select(date = week_end_date, state = jurisdiction, weekly_actual_days_reporting_any_data, weekly_percent_days_reporting_any_data, 
         total_admissions_all_covid_confirmed, total_admissions_adult_covid_confirmed, total_admissions_pediatric_covid_confirmed, 
         avg_admissions_all_covid_confirmed, percent_adult_covid_admissions,
         num_hospitals_admissions_all_covid_confirmed) %>%
  mutate(date = as.Date(date), 
         mmwr = MMWRweek::MMWRweek(date), 
         mmwr_date = paste0(mmwr$MMWRyear, '-', if_else(nchar(mmwr$MMWRweek) == 1, 
                                                        paste0("0", mmwr$MMWRweek), 
                                                        as.character(mmwr$MMWRweek))))



# National level
national_hosps <- hospitalizations_age %>%
  group_by(date, state) %>%
  summarise(hospitalizations = sum(total_admissions_all_covid_confirmed, na.rm = T)) %>%
  ungroup() %>%
  mutate(date = lubridate::ymd(date))


# Pull data from CDC API --------------------------------------------------
scrape_list = list(NWSS_Wastewater_Metric = 'https://data.cdc.gov/resource/2ew6-ywp6.csv', # 
                   NWSS_Concentration = 'https://data.cdc.gov/resource/g653-rqe2.csv',
                   NSSP_ED_Visit_Trajectory = 'https://data.cdc.gov/resource/rdmq-nq56.csv', # 
                   #MakeMyTestCount = 'https://data.cdc.gov/resource/i2a4-xk9k.csv', # 
                   NRVESS_Test_Positivity = 'https://data.cdc.gov/resource/gvsb-yw6g.csv') # 
                   #NVSN_child_Test_Positivity = 'https://data.cdc.gov/Public-Health-Surveillance/Percent-Positivity-of-Viral-Detections-Among-Enrol/kipu-qxy8.csv', 
                   #NRVESS_Var_Props = 'https://data.cdc.gov/resource/jr58-6ysp.csv')

full_scrape = lapply(scrape_list, function(x) try(pull_api(x)))



# Make plots --------------------------------------------------------------
full_scrape$MakeMyTestCount %>%
  filter(test_result == "Positive") %>%
  mutate(date = lubridate::ymd(date)) %>%
  group_by(date) %>%
  filter(date >= '2022-08-01') %>%
  summarise(total_tests = sum(total_tests)) %>%
  ungroup() %>%
  left_join(national_hosps) %>%
  ggplot(aes(x = date, y = total_tests)) +
  geom_line() +
  geom_point(aes(y = hospitalizations)) +
  scale_x_date(date_breaks = 'months')


#   
# tmaptools::geocode_OSM(participating_labs$institution, keep.unfound = T)
# addr <- paste0(server, "/search?q=", q2, "&format=xml&polygon=0&addressdetails=0")

  # separate_wider_delim(cols = df, delim = "\n", too_few = "align_start", names_sep = "_")


# crosswalk <- load_variables("acs1", year = 2022)
# 
# table <- c('B01001', str_c('B01001', LETTERS[1:9]))
# acs_data <- lapply(table, function(x){
#                   get_acs(geography = "county", table = x, year = 2022, survey = "acs1", 
#                         show_call = T)
#             }) %>%
#             bind_rows() %>% 
#             left_join(crosswalk, by = join_by("variable" == "name"))
# 
# write_rds(acs_data, 'acs_dat.RDS')
# acs_data <- read_rds('acs_dat.RDS') %>%
#               mutate(race = gsub())