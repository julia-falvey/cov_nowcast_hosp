# Load Packages -----------------------------------------------------------
pacman::p_load(aTSA, stats, forecast, prophet, urca, tidyverse, xtable, xgboost, 
               uroot)

# Functions ---------------------------------------------------------------
# Function to add random perturbation to predictors to avoid fit issues
add_perturbation <- function(x){
  if(sum(x, na.rm = T) == 0){
    return(log(runif(length(x))) + 1)
  } else {
    return(x)
  }
} 

## Faster than auto.arima way to get optimal model, first CSS then ML
extract_loglik <- function(p, d, q, P, D, Q){
  # print(paste("Order is", p, ",",d, ",", q, " | Season is", P, ",",D, ",", Q))
  arima(ts_dat, order = c(p, d, q), method = "CSS",
        seasonal = c(P, D, Q))$loglik
}
extract_aic <- function(p, d, q, P, D, Q){
  print(paste("Order is", p, ",",d, ",", q, " | Season is", P, ",",D, ",", Q))
  arima(ts_dat, order = c(p, d, q), method = "ML", transform.pars=FALSE,
        seasonal = c(P, D, Q))$aic
}

# Fit Arima for each state/univariate model
fit_arima <- function(cut_date, data, fit, seasonal, state, freq) {
  train_data <- data %>% filter(date < cut_date)
  ts_dat <- ts(train_data$hosp, frequency = freq)

  model <- tryCatch({Arima(ts_dat, order = fit,
                           method = "ML", transform.pars=FALSE,
                           seasonal = seasonal)}, 
                      error =  function(err){
                        mod <- Arima(ts_dat, order = fit, 
                                     transform.pars=FALSE,
                                     seasonal = seasonal)
                        return(mod)
                      })

  train_data$fitted <- model$fitted
  n = n_distinct(data %>% filter(date >= cut_date) %>% pull(date))
  forecasted <- as_tibble(forecast(model, h = n)) %>%
    bind_cols(data %>% filter(date >= cut_date)) %>%
    select(date, hosp, everything())

  train_data <- train_data %>%
    bind_rows(forecasted) %>%
    mutate(state = state, cut_date = cut_date, 
           fit = list(fit), 
           seasonal = list(seasonal), 
           freq = freq) %>%
    select(state, cut_date, fit, date, hosp, everything())

  return(train_data)
}

# Fit Arima for each state/multivariate model
fit_arima_regressors <- function(cut_date, data, state, freq, xreg = NULL) {
  train_data <- data %>% filter(date < cut_date)
  ts_dat <- ts(train_data$hosp_transformed, frequency = freq)
  
  training_xreg = xreg %>% filter(date < cut_date) %>%
    select(!date) %>%
    as.matrix()
  
  test_xreg = xreg %>% filter(date >= cut_date) %>%
    select(!date) %>%
    as.matrix()
  
  tryCatch({
    model <- auto.arima(ts_dat, max.order = 30, 
                        # Max values come from arima fits of TS alone
                        max.p = 7, max.q = 8, 
                        max.P = 3, max.Q = 3,
                        max.d = 2, max.D = 2,
                        stepwise = F, approximation = F, 
                        trace = T, ic = "aic", seasonal.test = "ch", 
                        xreg = training_xreg)
    
    train_data$fitted <- model$fitted
    n = n_distinct(data %>% filter(date >= cut_date) %>% pull(date))
    forecasted <- as_tibble(forecast(model, h = n, xreg = test_xreg)) %>%
      bind_cols(data %>% filter(date >= cut_date)) %>%
      select(date, hosp, everything())
    
    train_data <- train_data %>%
      bind_rows(forecasted) %>%
      mutate(state = state, cut_date = cut_date, 
             freq = freq, 
             fitted = exp(fitted) - 1,
             `Point Forecast` = exp(`Point Forecast`) - 1) %>%
      select(state, cut_date, date, hosp, everything())
    return(train_data)}, 
    error =  function(err){
      print(err)
      return(tibble())
    })
}

# Pull in Hospitalization Data --------------------------------------------
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



# Set Prediction Horizons -------------------------------------------------
max_date <- max(hospitalizations_age$date)
prediction_horizons <- c(max_date %m-% months(3), max_date %m-% months(6), 
                         max_date %m-% months(9), max_date %m-% months(12)) %>%
  as.character()

# ARIMA -------------------------------------------------------------------
## Only maintain the 26 week frequency period
freq_period = c(26)
# Run ARIMA
for(freq in freq_period){
  df_list <- tibble()
  for(state_nm in c(state.abb, "DC")){
    print(state_nm)
    # Subset data for just the given state
    data <- hospitalizations_age %>%
      filter(state == state_nm) %>%
      mutate(hosp = log(total_admissions_all_covid_confirmed+1)) %>%
      dplyr::select(date, hosp) 
    # Make timeseries object
    ts_dat <- ts(data$hosp, frequency = freq)
    # Get `d`
    d = ndiffs(ts_dat)
    p = 1:12
    q = 1:12
    # Create DF of orders
    full_orders = expand_grid(p, d, q)
    # force seasonal difference D >= 1
    D = max(1, nsdiffs(ts_dat, test = "ch", max.D = 2))
    P = 0:3
    Q = 0:3
    # Create DF of seasonal orders
    full_seasonal = bind_rows(expand_grid(P, D, Q))
    # Find optimal model
    model_fit = expand_grid(full_orders, full_seasonal) %>%
      rowwise() %>%
      mutate(log_lik = possibly(extract_loglik, otherwise = NA_real_)(p, d, q, P, D, Q)) %>%
      ungroup() %>%
      filter(!is.na(log_lik)) %>%
      slice_min(order_by=log_lik, n = 25) %>%
      rowwise() %>%
      mutate(aic = possibly(extract_aic, otherwise = NA_real_)(p, d, q, P, D, Q))%>%
      ungroup() %>%
      filter(aic == min(aic, na.rm = T)) %>% 
      ungroup() %>%
      # if there is a tie somehow, pick one at random
      slice_sample(n = 1)
    # Get fit data
    fit = c(model_fit$p, model_fit$d, model_fit$q)
    # If no model, fit with order(1,1,1)
    if(length(fit) == 0){
      fit = c(1,1,1)
    }
    # Get seasonal data
    seasonal = c(model_fit$P, model_fit$D, model_fit$Q)
    # If no seasonal model, fit with seasonal(1,1,1)
    if(length(seasonal) == 0){
      seasonal = c(1, 1, 1)
    }
    rm(model_fit)
    gc()
    # Run arima for the prediction horizons
    for(predict_cut in prediction_horizons){
      print(predict_cut)
      rmse_spec = tryCatch({fit_arima(predict_cut, data, fit, seasonal, state_nm, freq)}, 
                           # If fails, fit ARIMA with seasonal(1,0,1) model
                           error =  function(err){
                             fit_arima(predict_cut, data, fit, c(1,0,1), state_nm, freq)
                           })
      df_list <- bind_rows(df_list, rmse_spec %>%
                             mutate(hosp = hosp - 1, fitted = fitted - 1, 
                             `Point Forecast` = `Point Forecast` - 1))
    }
  }
  # Write out results
  write_rds(df_list, 
            paste0("results/arima_", freq, "_seasonal_fit.rds"))
  # Plot results
  for(predict_cut in prediction_horizons){
    p <- df_list %>%
      select(state:`Point Forecast`) %>%
      arrange(state) %>%
      rename(`Fitted Period` = fitted, `Forecast Period` = `Point Forecast`) %>%
      pivot_longer(`Fitted Period`:`Forecast Period`) %>%
      mutate(value = exp(value) - 1, 
             hosp = exp(hosp) - 1) %>%
      drop_na() %>%
      filter(cut_date == predict_cut, 
             date >= "2022-07-01") %>%
      ggplot(aes(x = date)) + 
      geom_line(aes(y = hosp, color = "Reported \nHospitalizations")) +
      geom_point(aes(y = value, color = name), size = 0.5) + 
      theme_minimal() + 
      scale_color_manual(values = c("#1B9E77", "#7570B3", "#000000")) + 
      labs(y = "Hospitalizations", x = "Date", color = "") +
      facet_wrap(state ~ ., ncol =  4, scales = "free_y")
    ggsave(paste0("results/Arima seasonal results for ", freq, 
                  " freqency time series at ", predict_cut, " horizon", ".pdf"),
           plot = p, width = 16, height = 10, units = "in", bg = "white", 
           dpi = "retina")
  }
}

# ARIMAX ------------------------------------------------------------------
full_data <- read_csv('data/updated_joined_df.csv')
ccf <- read_csv('results/ccf.csv') 

freq = 26
df_list <- tibble()

# For each state
for(state_nm in df_list$state){
  print(state_nm)
  # Pull CCF
  ccf_state = ccf %>%
    filter(state == state_nm)
  data <- full_data %>%
    filter(state == state_nm, 
           date <= '2024-05-01') %>%
    janitor::clean_names() %>%
    dplyr::select(date, hosp = total_admissions_all_covid_confirmed, count_pos, 
                  percent_pos, percent_visits_covid, 
                  weighted_raw_ww = avg_detect_prop_weighted_raw_wastewater,
                  weighted_postgrit_ww = avg_detect_prop_weighted_post_grit_removal,
                  weighted_sludge_ww = avg_detect_prop_weighted_primary_sludge) %>%
    # use box-cox transform on hospital data to avoid negatives in predictions
    # we assume box-cox order 0 (so log() transform, but add 1 to avoid fitting
    # log(0) in our model)
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
           weighted_sludge_ww = lag(replace_na(weighted_sludge_ww, 0), 
                                    n = ccf_state %>% 
                                      filter(type == 'Weighted WW Primary Sludge % Detect') %>%
                                      pull(lag) %>% abs()), 
           sludge_pct_transformed = log(weighted_sludge_ww + 1),
           hosp_transformed = log(hosp+1)) %>%
    filter(           date >= '2022-10-01') %>%
    drop_na()  %>%
    mutate(across(count_transformed:sludge_pct_transformed, 
                  ~ add_perturbation(.x)))
  for(predict_cut in prediction_horizons){
    # Run models for each prediction horizon and set of covariates
    print(paste0(state_nm, predict_cut))
    ed_reg = fit_arima_regressors(predict_cut, data, state_nm, freq,
                                xreg = data %>%
                                  select(date, ed_transformed)) %>%
      mutate(model_type = "% ED Visits")
    ed_count_reg = fit_arima_regressors(predict_cut, data, state_nm, freq,
                             xreg = data %>%
                               select(date, ed_transformed,
                                      count_transformed,
                                      percent_transformed))%>%
      mutate(model_type = "% ED Visits and Test Positivity")
    ww_reg = fit_arima_regressors(predict_cut, data, state_nm, freq,
                       xreg = data %>%
                         select(date, raw_pct_transformed)) %>%
      mutate(model_type = "Wastewater % Detection (Raw)")
    ww_all_reg = fit_arima_regressors(predict_cut, data, state_nm, freq, 
                           xreg = data %>%
                             select(date, raw_pct_transformed, 
                                    postgrit_pct_transformed, 
                                    sludge_pct_transformed)) %>% 
      mutate(model_type = "Wastewater % Detection (Raw, Post-Grit, Sludge)")
    all_reg = fit_arima_regressors(predict_cut, data, state_nm, freq, 
                        xreg = data %>%
                          select(date, ed_transformed, 
                                 count_transformed,
                                 percent_transformed,
                                 raw_pct_transformed, 
                                 postgrit_pct_transformed, 
                                 sludge_pct_transformed)) %>% 
      mutate(model_type = "All Regressors")
    df_list <- bind_rows(df_list, ed_reg) %>%
      bind_rows(ed_count_reg) %>%
      bind_rows(ww_reg) %>%
      bind_rows(ww_all_reg) %>%
      bind_rows(all_reg)
  }
}

# Write out results
write_rds(df_list, paste0("results/Multivariate/arimaX_", freq, "final.rds"))

# Plot
for(predict_cut in prediction_horizons){
  for(type in unique(df_list$model_type)) {
    p <- df_list %>%
      select(state, cut_date, date, hosp, hosp_transformed, 
             fitted, `Point Forecast`, model_type) %>%
      arrange(state) %>%
      rename(`Fitted Period` = fitted, `Forecast Period` = `Point Forecast`) %>%
#      filter(is.na(`Forecast Period`) | `Forecast Period` <= max(hosp) * 100 ) %>%
      pivot_longer(`Fitted Period`:`Forecast Period`) %>%
      # mutate(value = exp(value) - 1,
      #        hosp = exp(hosp) - 1) %>%
      drop_na() %>%
      filter(cut_date == predict_cut,
             date >= "2022-07-01", 
             model_type == type) %>%
      ggplot(aes(x = date)) +
      geom_line(aes(y = hosp, color = "Reported \nHospitalizations")) +
      geom_point(aes(y = value, color = name), size = 0.5) +
      theme_minimal() +
      scale_color_manual(values = c("#1B9E77", "#7570B3", "#000000")) +
      labs(y = "Hospitalizations", x = "Date", color = "") +
      facet_wrap(state ~ ., ncol =  4, scales = "free_y")
    ggsave(paste0("results/Multivariate/Multivariate Plots/ArimaX seasonal results for ", gsub("%", "PCT", type, fixed = T),
                  " freqency time series at ", predict_cut, " horizon", ".pdf"),
           plot = p, width = 16, height = 10, units = "in", bg = "white",
           dpi = "retina")
  }

}
