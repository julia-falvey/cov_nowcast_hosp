# Load Packages -----------------------------------------------------------
pacman::p_load(plotly, tidyverse, gridExtra, RColorBrewer)

# Function ---------------------------------------------------------------
get_ccf <- function(df, state_nm){
  print(state_nm)
  #Subset full data for state
  filtered_df <- df %>% 
    filter(state == state_nm) 

  #Subset dataframe for NREVSS columsn
  nvss = filtered_df %>%
    select(state, date, percent_pos, count_pos, 
           total_admissions_all_covid_confirmed) %>%
    filter(date >= '2022-08-01') %>%
    mutate(percent_pos = replace_na(percent_pos, 0), 
           count_pos = replace_na(count_pos, 0))

  #Subset dataframe for ED Visit columsn
  ed = filtered_df %>%
    select(state, date, percent_visits_covid, 
           total_admissions_all_covid_confirmed)  %>%  
    filter(date >= '2022-10-01') %>%
    mutate(total_admissions_all_covid_confirmed = 
             replace_na(total_admissions_all_covid_confirmed, 0))

  #Subset dataframe for NWSS columsn
  nwss = filtered_df %>%
    select(state, date, 
           `avg_detect_prop_weighted_post grit removal`:`avg_detect_prop_unweighted_primary effluent`,
           total_admissions_all_covid_confirmed) %>%
    drop_na()

  # Run CCF
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

  #Join CCF Reults
  tibble("acf" = pct_pos$acf, "lag" = pct_pos$lag, "type" = '% Positive') %>%
    bind_rows(tibble("acf" = count_pos$acf, "lag" = count_pos$lag,
                     "type" = 'Count Positive')) %>%
    bind_rows(tibble("acf" = ed_visit$acf, "lag" = ed_visit$lag,
                     "type" = 'ED Visit')) %>%
    bind_rows(tibble("acf" = avg_detect_prop_weighted_post_grit_removal$acf,
                     "lag" = avg_detect_prop_weighted_post_grit_removal$lag,
                     "type" = 'Weighted WW Post-Grit Removal % Detect')) %>%
    bind_rows(tibble("acf" = avg_detect_prop_unweighted_raw_wastewater$acf,
                     "lag" = avg_detect_prop_unweighted_raw_wastewater$lag,
                     "type" = 'Weighted WW Raw % Detect')) %>%
    bind_rows(tibble("acf" = avg_detect_prop_weighted_primary_sludge$acf,
                     "lag" = avg_detect_prop_weighted_primary_sludge$lag,
                     "type" = 'Weighted WW Primary Sludge % Detect')) %>%
    # Replace NA and find absolute value
    mutate(acf = replace_na(acf, 0), 
           abs_acf = abs(acf)) %>%
    rowwise() %>%
    # Only look at lag, not lead
    filter(any(lag <= 0)) %>%
    ungroup() %>%
    arrange(desc(abs_acf)) %>%
    # Grab only the first row once sorted
    slice_head(n = 1, by = "type") %>%
    mutate(state = state_nm, 
           # Replace lag with 0 if ACF was NA (now 0)
           lag = if_else(acf == 0, 0, lag))
}


# Iterate through states --------------------------------------------------
state_ccf= lapply(states$state, function(st){
  if(st == "USA"){
    return()
  } else {
    get_ccf(joined_df, st)
  }
}) %>%
  bind_rows()

write.csv(state_ccf, 'results/ccf.csv')


# Plot Heatmap ------------------------------------------------------------
ccf_df <- read_csv('results/ccf.csv') %>% 
  select(state, type, acf) %>%
  ungroup()

p <- ggplot(ccf_df, aes(x = state, y = type, fill = acf)) + 
  geom_tile() + 
  geom_text(aes(label = if_else(acf == 0, '-', as.character(round(acf, 2)))), 
            size=2.25) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 12)) +
  scale_fill_distiller("Correlation", palette = "RdBu", limits = c(-1,1)) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 0.45))

ggsave(paste0("results/CCF Heatmap.pdf"),
       plot = p, width = 16, height = 10, units = "in", bg = "white", 
       dpi = "retina")
