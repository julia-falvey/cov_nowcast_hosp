# Load Packages -----------------------------------------------------------
pacman::p_load(tidyverse, pdftools, jsonlite, httr, lubridate, tabulapdf, MMWRweek)

# Function ---------------------------------------------------------------
pull_api = function(endpoint, limit = "5000000") {
  return(read_csv(paste0(endpoint, "?$limit=", limit)))
}

# Create State/Abb Crosswalk DF -------------------------------------------

states <- state.abb %>% bind_cols(state.name) %>%
  bind_rows(list(`...1` = 'DC', `...2` = 'District of Columbia')) %>%
  bind_rows(list(`...1` = 'USA', `...2` = 'United States'))
colnames(states) <- c('state', 'geography')


# Pull Hospitalization Data -----------------------------------------------
# Using new data.cdc.gov endpoint:
hospitalizations_age <- read.csv("https://data.cdc.gov/resource/aemt-mg7g.csv?$limit=500000") %>% 
  select(date = week_end_date, state = jurisdiction, weekly_actual_days_reporting_any_data, weekly_percent_days_reporting_any_data, 
         total_admissions_all_covid_confirmed, total_admissions_adult_covid_confirmed, total_admissions_pediatric_covid_confirmed, 
         avg_admissions_all_covid_confirmed, percent_adult_covid_admissions,
         num_hospitals_admissions_all_covid_confirmed) %>%
  mutate(date = as.Date(date), 
         mmwr = MMWRweek::MMWRweek(date), 
         mmwr_date = paste0(mmwr$MMWRyear, '-', if_else(nchar(mmwr$MMWRweek) == 1, 
                                                        paste0("0", mmwr$MMWRweek), 
                                                        as.character(mmwr$MMWRweek))),
           #add HHS region to hospitalizations to map NVSS TP
           level = case_when(grepl('CT|ME|MA|NH|RI|VT', state) ~ 'Region 1', 
                                    grepl('NJ|NY|PR|VI', state) ~ 'Region 2', 
                                    grepl('DE|DC|MD|PA|VA|WV', state) ~ 'Region 3',
                                    grepl('AL|FL|GA|KY|MS|NC|SC|TN', state) ~ 'Region 4',
                                    grepl('IL|IN|MI|MN|OH|WI', state) ~ 'Region 5', 
                                    grepl('AR|LA|NM|OK|TX', state) ~ 'Region 6', 
                                    grepl('IA|KS|MO|NE', state) ~ 'Region 7', 
                                    grepl('CO|MT|ND|SD|UT|WY', state) ~ 'Region 8', 
                                    grepl('AZ|CA|HI|NV|AS|GU|MP', state) ~ 'Region 9', 
                                    grepl('AK|ID|OR|WA', state) ~ 'Region 10', 
                                    state == 'USA' ~ 'National',
                                    T ~ 'UNMAPPED'))

# Pull data from CDC API --------------------------------------------------
scrape_list = list(NWSS_Wastewater_Metric = 'https://data.cdc.gov/resource/2ew6-ywp6.csv?$limit=5000000', # 
                   NSSP_ED_Visit_Trajectory = 'https://data.cdc.gov/resource/rdmq-nq56.csv?$limit=5000000', # 
                   NRVESS_Test_Positivity = 'https://data.cdc.gov/resource/gvsb-yw6g.csv?$limit=5000000')

full_scrape = lapply(scrape_list, function(x) try(pull_api(x)))

nvss_tp <- full_scrape$NRVESS_Test_Positivity %>% 
  mutate(date = as.Date(mmwrweek_end), 
         mmwr = MMWRweek::MMWRweek(date), 
         mmwr_date = paste0(mmwr$MMWRyear, '-', 
                            if_else(nchar(mmwr$MMWRweek) == 1, paste0("0", 
                                                                      mmwr$MMWRweek),  as.character(mmwr$MMWRweek)))) %>%
  select(date, level, percent_pos, number_tested, mmwr_date, posted) %>%
  mutate(count_pos = percent_pos * number_tested) %>%
  group_by(date, level) %>%
  filter(posted== max(posted), 
         date <= "2024-05-01") %>%
  ungroup() %>%
  distinct()

nssp <- full_scrape$NSSP_ED_Visit_Trajectory.csv  %>% 
  filter(county == "All", 
         week_end <= "2024-05-01") %>% 
  mutate(date = as.Date(week_end), 
         mmwr = MMWRweek::MMWRweek(date), 
         mmwr_date = paste0(mmwr$MMWRyear, '-', 
                            if_else(nchar(mmwr$MMWRweek) == 1, paste0("0", 
                                                                      mmwr$MMWRweek),  as.character(mmwr$MMWRweek)))) 

# Ensure full set of data for each state, date possible from the dataset
dates <- unique(nssp  %>% select(date, mmwr_date) %>% distinct() %>% 
                  filter(date <= '2024-05-01'))
nssp <- nssp %>% 
  full_join(states %>% 
              cross_join(dates)) %>%
  # Replace NA values with a 0
  mutate(percent_visits_covid = replace_na(percent_visits_covid, 0)) %>%
  select(ed_visit_date = date, state, geography, percent_visits_covid, mmwr_date)

nwss_metric <-  full_scrape$NWSS_Wastewater_Metric.csv %>% 
  filter(date_start <= '2024-05-01', 
         date_start >= '2022-08-01') %>%
  select(wwtp_jurisdiction:key_plot_id, population_served, date = date_end,
         detect_prop_15d) %>%
  distinct() %>%
  mutate(date = as.Date(date), 
         mmwr = MMWRweek::MMWRweek(date), 
         mmwr_date = paste0(mmwr$MMWRyear, '-', 
                            if_else(nchar(mmwr$MMWRweek) == 1, paste0("0", 
                             mmwr$MMWRweek),  as.character(mmwr$MMWRweek))))  %>%
  filter(!is.nan(detect_prop_15d)) %>%
  group_by(reporting_jurisdiction, key_plot_id, sample_location, mmwr_date) %>%
  summarise(population_served = mean(population_served, na.rm = T),
            detect_prop_15d = mean(na.omit(detect_prop_15d), na.rm = T)) %>%
  group_by(reporting_jurisdiction, sample_location, key_plot_id) %>%
  arrange(reporting_jurisdiction, key_plot_id, mmwr_date) %>%
  tidyr::fill(., detect_prop_15d)

nwss <- nwss_metric %>%
  filter(sample_location == 'Treatment plant')  %>%
  mutate(location = gsub("90_|89_", "",
                         str_split_i(key_plot_id, "Treatment plant_", 2)))

nwss_dates <- cross_join(states, 
               cross_join(nwss  %>% ungroup() %>% select(location) %>% distinct(), 
                        nwss  %>% ungroup() %>% select(mmwr_date) %>% distinct()))

nwss <- nwss %>%
  ungroup() %>%
  select(state = reporting_jurisdiction, key_plot_id, mmwr_date, location, 
         population_served, detect_prop_15d) # pcr_conc_lin, normalization,

nwss_metric_wide <- nwss %>%
  select(state, key_plot_id, mmwr_date, location, population_served, 
         detect_prop_15d) %>%
  filter(!is.nan(detect_prop_15d), !is.na(population_served)) %>%
  group_by(state, mmwr_date, location) %>%
  summarise(total_pop = sum(population_served, na.rm = T), 
            avg_detect_prop_weighted = weighted.mean(detect_prop_15d, population_served, 
                                                     na.rm = T), 
            avg_detect_prop_unweighted = mean(detect_prop_15d,  na.rm = T)) %>%
  select(!total_pop)

nwss_wide <- nwss_metric_wide %>%
  ungroup() %>%
  rename(geography = state) %>%
  full_join(nwss_dates) %>%
  group_by(state, location) %>%
  fill(avg_detect_prop_weighted:avg_detect_prop_unweighted, .direction = "down") %>%
  ungroup() %>%
  select(!geography) %>%
  filter(!is.na(state)) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)/100)) |>
  pivot_wider(names_from = location, 
              values_from = avg_detect_prop_weighted:avg_detect_prop_unweighted)

# Combine Data ------------------------------------------------------------
joined_df <- hospitalizations_age %>%
  left_join(nvss_tp, by = c('level', 'mmwr_date', 'date')) %>%
  left_join(nssp, by = c('state', 'mmwr_date')) %>%
  left_join(nwss_wide, by = c('state', 'mmwr_date')) %>%
  filter(date <= "2024-05-01")

write_csv('data/joined_df.csv')