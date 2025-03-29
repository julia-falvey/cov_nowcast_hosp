# Load Packages -----------------------------------------------------------
pacman::p_load(tidyverse, pdftools, jsonlite, httr, lubridate, tabulapdf, MMWRweek)

# Functions ---------------------------------------------------------------
pull_api = function(endpoint, limit = "5000000") {
  return(read_csv(paste0(endpoint, "?$limit=", limit)))
}

# Pull Hospitalization Data -----------------------------------------------
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
scrape_list = list(NWSS_Wastewater_Metric = 'https://data.cdc.gov/resource/2ew6-ywp6.csv?$limit=5000000', # 
                   NSSP_ED_Visit_Trajectory = 'https://data.cdc.gov/resource/rdmq-nq56.csv?$limit=5000000', # 
                   NRVESS_Test_Positivity = 'https://data.cdc.gov/resource/gvsb-yw6g.csv?$limit=5000000')

full_scrape = lapply(scrape_list, function(x) try(pull_api(x)))
