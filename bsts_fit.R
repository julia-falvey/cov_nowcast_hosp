pacman::p_load(bsts, stats, tidyverse, scales, xtable)
setwd("~/Desktop/git/cov_nowcast_hosp")


hospitalizations_age <<- read.csv("https://data.cdc.gov/resource/aemt-mg7g.csv?$limit=500000") %>% 
  dplyr::select(date = week_end_date, state = jurisdiction, weekly_actual_days_reporting_any_data, weekly_percent_days_reporting_any_data, 
                total_admissions_all_covid_confirmed, total_admissions_adult_covid_confirmed, total_admissions_pediatric_covid_confirmed, 
                avg_admissions_all_covid_confirmed, percent_adult_covid_admissions,
                num_hospitals_admissions_all_covid_confirmed) %>%
  mutate(date = as.Date(date), 
         mmwr = MMWRweek::MMWRweek(date), 
         mmwr_date = paste0(mmwr$MMWRyear, '-', if_else(nchar(mmwr$MMWRweek) == 1, 
                                                        paste0("0", mmwr$MMWRweek), 
                                                        as.character(mmwr$MMWRweek)))) %>%
  filter(date <= "2024-05-01") # Filter out optional reporting period 

max_date <- max(hospitalizations_age$date)
prediction_horizons <- c(max_date %m-% months(3), max_date %m-% months(6), 
                         max_date %m-% months(9), max_date %m-% months(12)) %>%
  as.character()


get_yhats <- function(fit){
  burn <- SuggestBurn(0.2, fit)
  X     <- fit$state.contributions
  niter <- dim(X)[1]
  ncomp <- dim(X)[2]
  nobs  <- dim(X)[3]
  
  # initialize final fit/residuals matrices with zeros
  predictions <- matrix(data = 0, nrow = niter - burn, ncol = nobs)
  p0 <- predictions
  
  comps <- seq_len(ncomp)
  for (comp in comps) {
    # pull out the state contributions for this component and transpose to
    # a niter x (nobs - burn) array
    compX <- X[-seq_len(burn), comp, ]
    
    # accumulate the predictions across each component
    predictions <- predictions + compX
  }
  
  return(predictions)
}

fit_bsts <- function(cut_date, data, state) {
  train_data <- data %>% filter(date < cut_date)
  n = n_distinct(data %>% filter(date >= cut_date) %>% pull(date))
  ss <- list()
  ss <- AddLocalLevel(list(),train_data$hosp_transformed)
  ss <- AddStudentLocalLinearTrend(list(), train_data$hosp_transformed)
  ss <- AddSeasonal(ss, train_data$hosp_transformed, nseasons = 52, season.duration = 1)
  #ss <- AddTrig(ss, train_data$hosp, period = 26, 1:3)
  model1 <- bsts(train_data$hosp_transformed,
                 state.specification = ss,
                # family = 'poisson',
                 niter = 1000)
  fitted <- get_yhats(model1) %>% 
    as_tibble(.name_repair = "unique") %>%
    mutate(sim_number = row_number()) %>%
    pivot_longer(!sim_number) %>%
    mutate(timepoint = as.numeric(gsub('...', '', name, fixed = T))) %>%
    group_by(timepoint) %>%
    summarise(fitted = exp(mean(value)))

  forecasted <- data %>% filter(date >= cut_date) %>% 
    select(date) %>%
    mutate(forecasted = predict(model1, horizon = n, 
                                burn = SuggestBurn(0.2, model1))$mean %>% 
             exp() - 1) %>% 
    inner_join(data %>% filter(date >= cut_date))

  train_data <- train_data %>% 
    mutate(timepoint = row_number()) %>%
    left_join(fitted) %>%
    select(!timepoint) %>%
    bind_rows(forecasted) %>%
    dplyr::select(date, hosp, fitted, `Point Forecast` = forecasted)

  test_metric <- forecasted %>%
    ungroup() %>% 
    mutate(mse = (hosp-forecasted)^2) %>%
    summarise(rmse = sqrt(sum(mse))/n()) %>%
    pull(rmse)
  
  train_data <- train_data %>%
    mutate(state = state, cut_date = cut_date) %>%
    dplyr::select(state, cut_date, date, hosp, everything())
  return(list(data = train_data, rmse = test_metric))
}


out_tibble <- tribble(~state, ~cut_date, ~rmse_val)
df_list <- tibble()
for(state_nm in c(state.abb, "DC")){
  print(state_nm)
  data <- hospitalizations_age %>%
    filter(state == state_nm) %>%
    dplyr::select(date, hosp = total_admissions_all_covid_confirmed) %>%
    # # use box-cox transform on hospital data to avoid negatives in predictions
    # # we assume box-cox order 0 (so log() transform, but add 1 to avoid fitting
    # # log(0) in our model)
    mutate(hosp_transformed = log(hosp+1))
  for(predict_cut in prediction_horizons){
    rmse_spec = fit_bsts(predict_cut, data, state_nm)
    df_list <- bind_rows(df_list, rmse_spec$data)
    out_tibble <- bind_rows(out_tibble, 
                            tibble_row(state = state_nm, 
                                       cut_date = predict_cut, 
                                       rmse_val = rmse_spec$rmse))
  }
}


write_rds(df_list, "bsts_data.rds")
df_list <- read_rds("bsts_data.rds")
for(predict_cut in prediction_horizons){
  p <- df_list %>%
    filter(cut_date == predict_cut, 
           date >= "2022-07-01") %>%
    rename(`Fitted Period` = fitted, `Forecast Period` = `Point Forecast`) %>%
    pivot_longer(`Fitted Period`:`Forecast Period`) %>%
    drop_na() %>%
    filter(cut_date == as.Date(predict_cut), 
           date >= "2022-07-01") %>%
    ggplot(aes(x = date)) + 
    geom_line(aes(y = hosp, color = "Reported \nHospitalizations")) +
    geom_point(aes(y = value, color = name), size = 0.5) + 
    theme_minimal() + 
    scale_color_manual(values = c("#1B9E77", "#7570B3", "#000000")) + 
    labs(y = "Hospitalizations", x = "Date", color = "") +
    facet_wrap(state ~ ., ncol =  4, scales = "free_y")
  ggsave(paste0("BSTS Results for ", as.Date(predict_cut), ".pdf"), 
         plot = p, width = 16, height = 10, units = "in", bg = "white", 
         dpi = "retina")
}


change_table <- out_tibble %>%
  group_by(state) %>%
  mutate(rmse_farthest_timepoint = round(rmse_val[row_number()==1], 1), 
         pct_change = scales::percent((rmse_val - rmse_val[row_number()==1])/rmse_val[row_number()==1], 
                                      accuracy = 2)) %>%
  dplyr::select(-rmse_val) %>%
  ungroup() %>%
  pivot_wider(names_from = "cut_date", values_from = "pct_change") %>%
  dplyr::select(-`2024-01-27`) %>%
  rename(`2024-01-27` = rmse_farthest_timepoint)


options(xtable.floating = FALSE)
options(xtable.timestamp = "change_table")
print(xtable(change_table), include.rownames=FALSE) 


# With Regressors ---------------------------------------------------------

fit_bsts_regressor <- function(cut_date, data, state, n_reg) {
  train_data <- data %>% filter(date < cut_date)
  new_data <- data %>% filter(date >= cut_date) %>% select(!date)
  n = n_distinct(data %>% filter(date >= cut_date) %>% pull(date))
  ss <- list()
  ss <- AddLocalLevel(list(),train_data$hosp_transformed)
  ss <- AddLocalLinearTrend(list(), train_data$hosp_transformed)
  ss <- AddSeasonal(ss, train_data$hosp_transformed, nseasons = 52, season.duration = 1)
  #ss <- AddDynamicRegression(ss, hosp_transformed ~ ., data = train_data %>% select(-date, -hosp))
  #ss <- AddTrig(ss, train_data$hosp, period = 26, 1:3)
  model1 <- bsts(hosp_transformed ~ .,
                 state.specification = ss,
                 data = train_data %>% select(-date, -hosp),
                 # family = 'poisson',
                 niter = 3000,
                 expected.model.size = n_reg)
  fitted <- get_yhats(model1) %>% 
    as_tibble(.name_repair = "unique") %>%
    mutate(sim_number = row_number()) %>%
    pivot_longer(!sim_number) %>%
    mutate(timepoint = as.numeric(gsub('...', '', name, fixed = T))) %>%
    group_by(timepoint) %>%
    summarise(fitted = exp(mean(value)))
  
  forecasted <- data %>% filter(date >= cut_date) %>% 
    mutate(forecasted = predict(model1, horizon = n, 
                                burn = SuggestBurn(0.2, model1), 
                                newdata = new_data %>% 
                                  select(-hosp, -hosp_transformed) %>% as.matrix())$mean %>% 
             exp() - 1) %>% 
    inner_join(data %>% filter(date >= cut_date))
  
  train_data <- train_data %>% 
    mutate(timepoint = row_number()) %>%
    left_join(fitted) %>%
    select(!timepoint) %>%
    bind_rows(forecasted) %>%
    dplyr::select(date, hosp, fitted, `Point Forecast` = forecasted)
  
  train_data <- train_data %>%
    mutate(state = state, cut_date = cut_date) %>%
    dplyr::select(state, cut_date, date, hosp, everything())
  return(train_data)
}


df_list <- tibble()
full_data <- read_csv('data/updated_joined_df.csv')
ccf = read_csv('results/ccf.csv') 
for(state_nm in c(state.abb, "DC")){
  ccf_state = ccf %>%
    filter(state == state_nm)
  print(state_nm)
  data <- full_data %>%
    filter(state == state_nm, 
           date <= '2024-05-01') %>%
    janitor::clean_names() %>%
    dplyr::select(date, hosp = total_admissions_all_covid_confirmed, count_pos, 
                  percent_pos, percent_visits_covid, 
                  weighted_raw_ww = avg_detect_prop_weighted_raw_wastewater,
                  weighted_postgrit_ww = avg_detect_prop_weighted_post_grit_removal,
                  microbial_postgrit = microbial_post_grit_removal, 
                  microbial_sludge = microbial_primary_sludge,
                  flow_postgrit = flow_population_post_grit_removal, 
                  flow_raw = flow_population_raw_wastewater) %>%
    # use box-cox transform on hospital data to avoid negatives in predictions
    # we assume box-cox order 0 (so log() transform, but add 1 to avoid fitting
    # log(0) in our model)
    # also, add random perturbation to all 0 predictors to match approach taken 
    # with ARIMA to avoid fitting failures 
      mutate(count_pos = lag(replace_na(count_pos, 0), 
                           n = ccf_state %>% 
                             filter(type == 'Count Positive') %>%
                             pull(lag) %>% abs()),
           count_transformed = log(count_pos+1),
           percent_pos = lag(replace_na(percent_pos, 0), 
                             n = ccf_state %>% 
                               filter(type == '% Positive') %>%
                               pull(lag) %>% abs()),
           percent_transformed = log(percent_pos+1),
           percent_visits_covid = lag(replace_na(percent_visits_covid, 0), 
                                      n = ccf_state %>% 
                                        filter(type == 'ED Visit') %>%
                                        pull(lag) %>% abs()),
           ed_transformed = log(percent_visits_covid + 1),
           weighted_raw_ww = lag(replace_na(weighted_raw_ww, 0), 
                                 n = ccf_state %>% 
                                   filter(type == 'Weighted WW Raw % Detect') %>%
                                   pull(lag) %>% abs()),
           raw_pct_transformed = log(weighted_raw_ww + 1),
           weighted_postgrit_ww = lag(replace_na(weighted_postgrit_ww, 0), 
                                      n = ccf_state %>% 
                                        filter(type == 'Weighted WW Post-Grit Removal % Detect') %>%
                                        pull(lag) %>% abs()), 
           postgrit_pct_transformed = log(weighted_postgrit_ww + 1),
           microbial_postgrit = lag(replace_na(microbial_postgrit, 0), 
                                    n = ccf_state %>% 
                                      filter(type == 'Microbial Post-Grit Removal WW Detect') %>%
                                      pull(lag) %>% abs()), 
           microbial_postgrit_transformed = log(microbial_postgrit + 1),
           microbial_sludge = lag(replace_na(microbial_sludge, 0), 
                                  n = ccf_state %>% 
                                    filter(type == 'Microbial Primary Sludge WW Detect') %>%
                                    pull(lag) %>% abs()),
           microbial_sludge_transformed = log(microbial_sludge + 1),
           flow_raw = lag(replace_na(flow_raw, 0), 
                          n = ccf_state %>% 
                            filter(type == 'Flow-Pop Raw WW Detect') %>%
                            pull(lag) %>% abs()),
           flow_raw_transformed = log(flow_raw + 1),
           flow_postgrit = lag(replace_na(flow_postgrit, 0), 
                               n = ccf_state %>% 
                                 filter(type == 'Flow-Pop Post-Grit Removal WW Detect') %>%
                                 pull(lag) %>% abs()),
           flow_postgrit_transformed = log(flow_postgrit + 1),
           hosp_transformed = log(hosp+1))  %>%
    filter(           date >= '2022-10-01') %>%
    drop_na()

  for(predict_cut in prediction_horizons){
    print(paste0(state_nm, "- ", predict_cut))
    ed_reg = fit_bsts_regressor(predict_cut, data %>%
                              select(date, hosp, hosp_transformed, ed_transformed), 
                            state_nm, n_reg = 1)  %>%
      mutate(model_type = "% ED Visits")
    ed_count_reg = fit_bsts_regressor(predict_cut, data %>%
                                          select(date, hosp, hosp_transformed, 
                                                 ed_transformed, 
                                                 count_transformed, 
                                                 percent_transformed), 
                                        state_nm, n_reg = 3) %>%
      mutate(model_type = "% ED Visits and Test Positivity")
    ww_reg = fit_bsts_regressor(predict_cut, data %>%
                              select(date, hosp, hosp_transformed,
                                     raw_pct_transformed),
                            state_nm, n_reg = 1) %>%
   mutate(model_type = "Wastewater % Detection")
    ww_all_reg = fit_bsts_regressor(predict_cut, data %>%
                                      select(date, hosp, hosp_transformed,
                                             raw_pct_transformed,
                                             postgrit_pct_transformed,
                                             microbial_postgrit_transformed,
                                             microbial_sludge_transformed,
                                             flow_raw_transformed,
                                             flow_postgrit_transformed),
                                    state_nm, n_reg = 6) %>%
      mutate(model_type = "Wastewater % Detection and Concentrations")
    all_reg = fit_bsts_regressor(predict_cut, data %>%
                                       select(date, hosp, hosp_transformed, 
                                              ed_transformed, 
                                              count_transformed,
                                              percent_transformed,
                                              raw_pct_transformed, 
                                              postgrit_pct_transformed, 
                                              microbial_postgrit_transformed,
                                              microbial_sludge_transformed, 
                                              flow_raw_transformed, 
                                              flow_postgrit_transformed), state_nm,
                                 n_reg = 9)  %>%
      mutate(model_type = "All Regressors")
    df_list <- bind_rows(df_list, ed_reg) %>%
      bind_rows(ed_count_reg) %>%
      bind_rows(ww_reg) %>%
      bind_rows(ww_all_reg) %>%
      bind_rows(all_reg)
  }
}


write_rds(df_list, "results/Multivariate/bsts_lagged_regressors_no_ww_outlier_model.rds")

for(predict_cut in prediction_horizons){
  for(type in unique(df_list$model_type)){
    p <- df_list %>%
      filter(cut_date == predict_cut, 
             model_type == type) %>%    
      rename(`Fitted Period` = fitted, `Forecast Period` = `Point Forecast`) %>%
      pivot_longer(`Fitted Period`:`Forecast Period`) %>%
      drop_na() %>%
      filter(cut_date == as.Date(predict_cut), 
             date >= "2022-07-01") %>%
      ggplot(aes(x = date)) + 
      geom_line(aes(y = hosp, color = "Reported \nHospitalizations")) +
      geom_point(aes(y = value, color = name), size = 0.5) + 
      theme_minimal() + 
      scale_color_manual(values = c("#1B9E77", "#7570B3", "#000000")) + 
      labs(y = "Hospitalizations", x = "Date", color = "") +
      facet_wrap(state ~ ., ncol =  4, scales = "free_y")
        
    ggsave(paste0("results/Multivariate/Multivariate Plots/BSTS Regressor Results for ", as.Date(predict_cut), "_", 
                  gsub("%", "PCT", type, fixed = T), ".pdf"), 
           plot = p, width = 16, height = 10, units = "in", bg = "white", 
           dpi = "retina")
  }
}
