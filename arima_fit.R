pacman::p_load(aTSA, stats, forecast, prophet, urca, tidyverse, xtable, xgboost, 
               uroot)


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

fit_arima <- function(cut_date, data, fit, seasonal, state) {
  train_data <- data %>% filter(date < cut_date)
  model <- forecast::Arima(train_data$hosp, order=fit,
                           seasonal = seasonal, method = "ML",
                 include.drift = FALSE )
  train_data$fitted <- model$fitted
  n = n_distinct(data %>% filter(date >= cut_date) %>% pull(date))
  forecasted <- as_tibble(forecast(model, h = n)) %>%
    bind_cols(data %>% filter(date >= cut_date)) %>%
    select(date, hosp, everything())
  
  test_metric <- forecasted %>%
    ungroup() %>% 
    mutate(sse = (hosp-as_tibble(forecast(model, h = n))$`Point Forecast`)^2) %>%
    summarise(rmse = sqrt(sum(sse)/n())) %>%
    pull(rmse)
  
  train_data <- train_data %>%
    bind_rows(forecasted) %>%
    mutate(state = state, cut_date = cut_date, 
           fit = list(fit), 
           seasonal = list(seasonal)) %>%
    select(state, cut_date, fit, date, hosp, everything())
  return(list(data = train_data, rmse = test_metric))
}

# start_point <- read_rds("results/arima_stepwise_non_seasonal_fit.rds") %>% 
#   dplyr::select(state, fit) %>%
#   unnest_wider(fit, names_sep = "_") %>% distinct() 

extract_loglik <- function(p, d, q, P, D, Q){
  print(paste("Order is", p, ",",d, ",", q, " | Season is", P, ",",D, ",", Q))
  arima(ts_dat, order = c(p, d, q), method = "CSS",
        seasonal = c(P, D, Q))$loglik
}

extract_aic <- function(p, d, q, P, D, Q){
  print(paste("Order is", p, ",",d, ",", q, " | Season is", P, ",",D, ",", Q))
  arima(ts_dat, order = c(p, d, q), method = "ML", transform.pars=FALSE,
        seasonal = c(P, D, Q))$aic
}
  
out_tibble <- tribble(~state, ~cut_date, ~rmse_val)
df_list <- tibble()
for(state_nm in lubridate::setdiff(c(state.abb, "DC"), c(unique(df_list$state), "UT"))){
  print(state_nm)
  data <- hospitalizations_age %>%
    filter(state == state_nm) %>%
    select(date, hosp = total_admissions_all_covid_confirmed) 
  ts_dat <- ts(data$hosp, frequency = 52)
  d = ndiffs(ts_dat)
  D = nsdiffs(ts_dat, test = "hegy", max.D = 2)
  p = 0:12
  q = 0:12
  P = 0:3
  Q = 0:3
  full_orders = expand_grid(p, d, q)
  full_seasonal = bind_rows(expand_grid(P, D, Q), tribble(~P, ~D, ~Q, 0, 0, 0))
  model_fit = expand_grid(full_orders, full_seasonal) %>%
    rowwise() %>%
    mutate(log_lik = possibly(extract_loglik, otherwise = NA_real_)(p, d, q, P, D, Q)) %>%
    ungroup() %>%
    filter(!is.na(log_lik)) %>%
    slice_min(order_by=log_lik, prop = 0.05) %>%
    rowwise() %>%
    mutate(aic = possibly(extract_aic, otherwise = NA_real_)(p, d, q, P, D, Q))%>%
    ungroup() %>%
    filter(aic == min(aic, na.rm = T)) %>% 
    mutate(order_sum = p+q, 
           p = if_else(order_sum == 0, 1, p),
           q = if_else(order_sum == 0, 1, q), 
           season_sum = P+Q, 
           P = if_else(season_sum == 0, 1, P), 
           Q = if_else(season_sum == 0, 1, Q)) %>%
    ungroup() %>%
    # if there is a tie somehow, pick one at random
    slice_sample(n = 1)
  fit = c(model_fit$p, model_fit$d, model_fit$q)
  seasonal = c(model_fit$P, model_fit$D, model_fit$Q)
  rm(model_fit)
  gc()
  for(predict_cut in prediction_horizons){
    print(predict_cut)
    rmse_spec = fit_arima(predict_cut, data, fit, seasonal, state_nm)
    df_list <- bind_rows(df_list, rmse_spec$data)
    out_tibble <- bind_rows(out_tibble, 
                            tibble_row(state = state_nm, 
                                       cut_date = predict_cut, 
                                       rmse_val = rmse_spec$rmse))
  }
}

write_rds(df_list, "results/arima_stepwise_seasonal_fit_hegy_nsdiffs.rds")

for(predict_cut in prediction_horizons){
  p <- df_list %>%
    select(state:`Point Forecast`) %>%
    arrange(state) %>%
    rename(`Fitted Period` = fitted, `Forecast Period` = `Point Forecast`) %>%
    pivot_longer(`Fitted Period`:`Forecast Period`) %>%
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
  ggsave(paste0("Arima seasonal results for ch ", predict_cut, ".pdf"), 
         plot = p, width = 16, height = 10, units = "in", bg = "white", 
         dpi = "retina")
}


change_table <- out_tibble %>%
  group_by(state) %>%
  arrange(cut_date) %>%
  mutate(rmse_farthest_timepoint = round(rmse_val[row_number()==1], 1), 
         pct_change = scales::percent((rmse_val - rmse_val[row_number()==1])/rmse_val[row_number()==1], 
                                      accuracy = 2)) %>%
  dplyr::select(-rmse_val) %>%
  ungroup() %>%
  pivot_wider(names_from = "cut_date", values_from = "pct_change") %>%
  dplyr::select(-`2023-04-27`) %>%
  rename(`2023-04-27` = rmse_farthest_timepoint)


options(xtable.floating = FALSE)
options(xtable.timestamp = "change_table")
print(xtable(change_table), include.rownames=FALSE) 
