pacman::p_load(plotly, tidyverse, gridExtra, RColorBrewer)
# CCF of data and outcome

## Combine datasets
hospitalizations_age <- read_csv('data/hospitalizations.csv') %>%
  #add HHS region to hospitalizations to map NVSS TP
  mutate(level = case_when(grepl('CT|ME|MA|NH|RI|VT', state) ~ 'Region 1', 
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


states <- state.abb %>% bind_cols(state.name) %>%
  bind_rows(list(`...1` = 'DC', `...2` = 'District of Columbia')) %>%
  bind_rows(list(`...1` = 'USA', `...2` = 'United States'))
colnames(states) <- c('state', 'geography')


nvss_tp <- read_csv('data/NRVESS_Test_Positivity.csv') %>% 
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

nssp <- read_csv('data/NSSP_ED_Visit_Trajectory.csv')  %>% 
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

nwss_metric <-  read_csv('data/NWSS_Wastewater_Metric.csv') %>% 
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

nwss_concentration <- read_csv('data/NWSS_Concentration.csv') %>% 
  filter(date <= '2024-05-01', 
         date >= '2022-08-01') %>%
  mutate(date = as.Date(date), 
         mmwr = MMWRweek::MMWRweek(date), 
         mmwr_date = paste0(mmwr$MMWRyear, '-', 
                            if_else(nchar(mmwr$MMWRweek) == 1, paste0("0", 
                        mmwr$MMWRweek),  as.character(mmwr$MMWRweek)))) %>%
  group_by(key_plot_id, normalization, mmwr_date) %>%
  summarise(pcr_conc_lin = sum(pcr_conc_lin, na.rm = T))
  
nwss <- full_join(nwss_metric, nwss_concentration) %>%
  filter(sample_location == 'Treatment plant')  %>%
  mutate(location = gsub("90_|89_", "",
                         str_split_i(key_plot_id, "Treatment plant_", 2)))


# Ensure full set of data for each state, date, location possible from the dataset
nwss_dates <- cross_join(states, 
              cross_join(nwss  %>% ungroup() %>% select(location) %>% distinct(), 
                         nwss  %>% ungroup() %>% select(mmwr_date) %>% distinct()))

nwss <- nwss %>%
  ungroup() %>%
  select(state = reporting_jurisdiction, key_plot_id, mmwr_date, location, 
         pcr_conc_lin, population_served, normalization, detect_prop_15d) 

nwss_concentration_wide = nwss %>% 
  select(state, key_plot_id, mmwr_date, location, normalization, pcr_conc_lin) %>%
  group_by(state, key_plot_id, location) %>%
  fill(normalization) %>%
  filter(!is.nan(pcr_conc_lin), !is.na(normalization)) %>%
  group_by(state, mmwr_date, location, normalization) %>%
  summarise(avg_detect = mean(pcr_conc_lin, na.rm = T))   %>%
  pivot_wider(names_from = normalization, values_from = avg_detect, 
              values_fill = 0) %>%
  # Update to handle outliers in concentration reporting that cause issues - 
  # i.e. one day there's a value that is orders of magnitude higher than all
  # other days - might be week new plant is coming online or data error, 
  # but keeping in causes major issues
  group_by(state, location) %>%
  mutate(microbial = if_else(microbial >= quantile(microbial, 0.95, na.rm = T) * 10, 
                             0, microbial), 
         `flow-population` = if_else(`flow-population` >= quantile(`flow-population`, 0.95, na.rm = T) * 10, 
                             0, `flow-population`))

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

nwss_wide <- nwss_concentration_wide %>%
  full_join(nwss_metric_wide) %>% 
  ungroup() %>%
  rename(geography = state) %>%
  full_join(nwss_dates) %>%
  group_by(state, location) %>%
  fill(microbial:avg_detect_prop_unweighted, .direction = "down") %>%
  ungroup() %>%
  select(!geography) %>%
  filter(!is.na(state)) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) |>
  pivot_wider(names_from = location, 
              values_from = microbial:avg_detect_prop_unweighted)

joined_df <- hospitalizations_age %>%
  left_join(nvss_tp, by = c('level', 'mmwr_date', 'date')) %>%
  left_join(nssp, by = c('state', 'mmwr_date')) %>%
  left_join(nwss_wide, by = c('state', 'mmwr_date')) %>%
  filter(date <= "2024-05-01")

write_csv(joined_df, "data/full_df_nwss_fixed.csv")

get_ccf <- function(df, state_nm){
  print(state_nm)
  filtered_df <- df %>% 
    filter(state == state_nm) 
  
  nvss = filtered_df %>%
    select(state, date, percent_pos, count_pos, 
           total_admissions_all_covid_confirmed) %>%
    filter(date >= '2022-08-01') %>%
    mutate(percent_pos = replace_na(percent_pos, 0), 
           count_pos = replace_na(count_pos, 0))
  
  ed = filtered_df %>%
    select(state, date, percent_visits_covid, 
           total_admissions_all_covid_confirmed)  %>%  
    filter(date >= '2022-10-01') %>%
    mutate(total_admissions_all_covid_confirmed = 
             replace_na(total_admissions_all_covid_confirmed, 0))
  
  
  nwss = filtered_df %>%
    select(state, date, `microbial_post grit removal`:`avg_detect_prop_unweighted_primary effluent`,
           total_admissions_all_covid_confirmed) %>%
    drop_na()
  
  pct_pos = ccf(nvss$percent_pos,
                nvss$total_admissions_all_covid_confirmed,
                lag.max = 5,
                plot = F)
  count_pos = ccf(nvss$count_pos,
                  nvss$total_admissions_all_covid_confirmed, 
                  lag.max = 5,
                  na.action = na.pass, plot = F)
  
   ed_visit = ccf(ed$percent_visits_covid,
                   ed$total_admissions_all_covid_confirmed, lag.max = 5,
                   na.action = na.pass, plot = F)


    microbial_post_grit_removal = ccf(nwss$`microbial_post grit removal`,
                             nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                             na.action = na.pass, plot = F)
    microbial_raw_wastewater = ccf(nwss$`microbial_raw wastewater`,
                                      nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                      na.action = na.pass, plot = F)
    microbial_primary_sludge = ccf(nwss$`microbial_primary sludge`,
                                   nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                   na.action = na.pass, plot = F)
    microbial_primary_effluent = ccf(nwss$`microbial_primary effluent`,
                                   nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                   na.action = na.pass, plot = F)
    
    flow_population_post_grit_removal = ccf(nwss$`flow-population_post grit removal`,
                                      nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                      na.action = na.pass, plot = F)
    flow_population_raw_wastewater = ccf(nwss$`flow-population_raw wastewater`,
                                   nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                   na.action = na.pass, plot = F)
    flow_population_primary_sludge = ccf(nwss$`flow-population_primary sludge`,
                                   nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                   na.action = na.pass, plot = F)
    flow_population_primary_effluent = ccf(nwss$`flow-population_primary effluent`,
                                     nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                     na.action = na.pass, plot = F)
    
    avg_detect_prop_weighted_post_grit_removal = ccf(nwss$`avg_detect_prop_weighted_post grit removal`,
                                      nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                      na.action = na.pass, plot = F)
    avg_detect_prop_weighted_raw_wastewater = ccf(nwss$`avg_detect_prop_weighted_raw wastewater`,
                                   nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                   na.action = na.pass, plot = F)
    avg_detect_prop_weighted_primary_sludge = ccf(nwss$`avg_detect_prop_weighted_primary sludge`,
                                   nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                   na.action = na.pass, plot = F)
    avg_detect_prop_weighted_primary_effluent = ccf(nwss$`avg_detect_prop_weighted_primary effluent`,
                                     nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                     na.action = na.pass, plot = F)
    
    
    avg_detect_prop_unweighted_post_grit_removal = ccf(nwss$`avg_detect_prop_unweighted_post grit removal`,
                                                     nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                                     na.action = na.pass, plot = F)
    avg_detect_prop_unweighted_raw_wastewater = ccf(nwss$`avg_detect_prop_unweighted_raw wastewater`,
                                                  nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                                  na.action = na.pass, plot = F)
    avg_detect_prop_unweighted_primary_sludge = ccf(nwss$`avg_detect_prop_unweighted_primary sludge`,
                                                  nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                                  na.action = na.pass, plot = F)
    avg_detect_prop_unweighted_primary_effluent = ccf(nwss$`avg_detect_prop_unweighted_primary effluent`,
                                                    nwss$total_admissions_all_covid_confirmed, lag.max = 5,
                                                    na.action = na.pass, plot = F)


  tibble("acf" = pct_pos$acf, "lag" = pct_pos$lag, "type" = '% Positive') %>%
    bind_rows(tibble("acf" = count_pos$acf, "lag" = count_pos$lag,
                     "type" = 'Count Positive')) %>%
    bind_rows(tibble("acf" = ed_visit$acf, "lag" = ed_visit$lag,
                     "type" = 'ED Visit')) %>%
    bind_rows(tibble("acf" = microbial_post_grit_removal$acf,
                     "lag" = microbial_post_grit_removal$lag,
                     "type" = 'Microbial Post-Grit Removal WW Detect')) %>%
    # bind_rows(tibble("acf" = microbial_raw_wastewater$acf,
    #                  "lag" = microbial_raw_wastewater$lag,
    #                  "type" = 'Microbial Raw WW Detect')) %>%
    bind_rows(tibble("acf" = microbial_primary_sludge$acf,
                     "lag" = microbial_primary_sludge$lag,
                     "type" = 'Microbial Primary Sludge WW Detect')) %>%
    # bind_rows(tibble("acf" = microbial_primary_effluent$acf,
    #                  "lag" = microbial_primary_effluent$lag,
    #                  "type" = 'Microbial Primary Effluent WW Detect')) %>%
    bind_rows(tibble("acf" = flow_population_post_grit_removal$acf,
                     "lag" = flow_population_post_grit_removal$lag,
                     "type" = 'Flow-Pop Post-Grit Removal WW Detect')) %>%
    bind_rows(tibble("acf" = flow_population_raw_wastewater$acf,
                     "lag" = flow_population_raw_wastewater$lag,
                     "type" = 'Flow-Pop Raw WW Detect')) %>%
    # bind_rows(tibble("acf" = flow_population_primary_sludge$acf,
    #                  "lag" = flow_population_primary_sludge$lag,
    #                  "type" = 'Flow-Pop Primary Sludge WW Detect')) %>%
    # bind_rows(tibble("acf" = flow_population_primary_effluent$acf,
    #                  "lag" = flow_population_primary_effluent$lag,
    #                  "type" = 'Flow-Pop Primary Effluent WW Detect')) %>%
    bind_rows(tibble("acf" = avg_detect_prop_weighted_post_grit_removal$acf,
                     "lag" = avg_detect_prop_weighted_post_grit_removal$lag,
                     "type" = 'Weighted WW Post-Grit Removal % Detect')) %>%
    bind_rows(tibble("acf" = avg_detect_prop_unweighted_raw_wastewater$acf,
                     "lag" = avg_detect_prop_unweighted_raw_wastewater$lag,
                     "type" = 'Weighted WW Raw % Detect')) %>%
    # bind_rows(tibble("acf" = avg_detect_prop_weighted_primary_sludge$acf,
    #                  "lag" = avg_detect_prop_weighted_primary_sludge$lag,
    #                  "type" = 'Weighted WW Primary Sludge % Detect')) %>%
    # bind_rows(tibble("acf" = avg_detect_prop_weighted_primary_effluent$acf,
    #                  "lag" = avg_detect_prop_weighted_primary_effluent$lag,
    #                  "type" = 'Weighted WW Primary Effluent % Detect')) %>%
    
    # bind_rows(tibble("acf" = avg_detect_prop_unweighted_post_grit_removal$acf,
    #                  "lag" = avg_detect_prop_unweighted_post_grit_removal$lag,
    #                  "type" = 'Unweighted WW Post-Grit Removal % Detect')) %>%
    # bind_rows(tibble("acf" = avg_detect_prop_unweighted_raw_wastewater$acf,
    #                  "lag" = avg_detect_prop_unweighted_raw_wastewater$lag,
    #                  "type" = 'Unweighted WW Raw % Detect')) %>%
    # bind_rows(tibble("acf" = avg_detect_prop_unweighted_primary_sludge$acf,
    #                  "lag" = avg_detect_prop_unweighted_primary_sludge$lag,
    #                  "type" = 'Unweighted WW Primary Sludge % Detect')) %>%
    # bind_rows(tibble("acf" = avg_detect_prop_unweighted_primary_effluent$acf,
    #                  "lag" = avg_detect_prop_unweighted_primary_effluent$lag,
    #                  "type" = 'Unweighted WW Primary Effluent % Detect')) %>%
    mutate(acf = replace_na(acf, 0), 
           abs_acf = abs(acf)) %>%
    rowwise() %>%
    filter(any(lag <= 0)) %>%
    ungroup() %>%
    arrange(desc(abs_acf)) %>%
    slice_head(n = 1, by = "type") %>%
    mutate(state = state_nm, 
           lag = if_else(acf == 0, 0, lag))
}

state_ccf= lapply(states$state, function(st){
  if(st == "USA"){
    return()
  } else {
    get_ccf(joined_df, st)
  }
}) %>%
  bind_rows()


ccf_df <- state_ccf %>% 
  select(state, type, acf) %>%
  ungroup()  # %>%
#   pivot_wider(names_from = "type", values_from = "acf") %>%
#   as.data.frame()
# 
# rownames(ccf_df) = ccf_df$state
# ccf_df <- as.matrix(select(ccf_df, !state))
# 
# heatmap(ccf_df, Colv = NA, Rowv = NA, scale="column", 
#         col= colorRampPalette(brewer.pal(8, "YlGnBu"))(25))
#         


write.csv(state_ccf, 'results/ccf.csv')
write_csv(joined_df, 'data/updated_joined_df.csv')

ccf_df <- read_csv('results/ccf.csv') %>% 
  select(state, type, acf) %>%
  ungroup()
p <- ggplot(ccf_df, aes(x = state, y = type, fill = acf)) + 
  geom_tile() + 
  geom_text(aes(label = if_else(acf == 0, '-', as.character(round(acf, 2)))), size=2.25) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 12)) +
  scale_fill_distiller("Correlation", palette = "RdBu", limits = c(-1,1)) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 0.45))
ggsave(paste0("results/CCF Heatmap.pdf"),
       plot = p, width = 16, height = 10, units = "in", bg = "white", 
       dpi = "retina")

