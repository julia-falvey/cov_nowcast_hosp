pacman::p_load(aTSA, stats, forecast, prophet, urca, tidyverse, xtable, xgboost)


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

start_point <- read_rds("results/arima_stepwise_non_seasonal_fit.rds") %>% 
  dplyr::select(state, fit) %>%
  unnest_wider(fit, names_sep = "_") %>% distinct() 

out_tibble <- tribble(~state, ~cut_date, ~rmse_val)
df_list <- tibble()
for(state_nm in lubridate::setdiff(c(state.abb, "DC"), c(unique(df_list$state), "UT"))){
  print(state_nm)
  data <- hospitalizations_age %>%
    filter(state == state_nm) %>%
    select(date, hosp = total_admissions_all_covid_confirmed) 
  ts_dat <- ts(data$hosp, frequency = 52)
  d = ndiffs(ts_dat)
  D = nsdiffs(ts_dat)
  model_fit = auto.arima(ts_dat, d = d, allowdrift = F,
                         max.order = 34, seasonal = T,
                         max.p = 6, max.q = 6, max.P = 3, max.Q = 3,
                         start.p = start_point$fit_1, start.P = 1,
                         start.q = start_point$fit_3,
                         stepwise = F, approximation = T, trace = T) %>%
    broom::tidy()
  
  if(nrow(model_fit) == 0){
    ar = 1
    ma = 1
    sar = 1
    sam = 1
  } else{
    seasonal <- model_fit %>%
      filter(grepl('s', term)) %>%
      mutate(order = as.integer(str_sub(term, 4)),
             term = str_sub(term, 2, 3))
      sar = seasonal %>% filter(term == 'ar') %>%
        filter(order == max(order, 0)) %>% pull(order)
      sam = seasonal %>% filter(term == 'ma') %>%
        filter(order == max(order, 0)) %>% pull(order)
    model_fit <- model_fit %>% 
      filter(!grepl('s', term)) %>% 
      mutate(order = as.integer(str_sub(term, 3)), 
           term = str_sub(term, 1, 2))
    ar = model_fit %>% filter(term == 'ar') %>% 
      filter(order == max(order, 0)) %>% pull(order)
    ma = model_fit %>% filter(term == 'ma')  %>% 
      filter(order == max(order, 0)) %>% pull(order)
  }
  fit = c(max(ar, 0, na.rm = T), d, max(ma, 0, na.rm = T))
  seasonal = c(max(sar, 0, na.rm = T), D, max(sam, 0, na.rm = T))
  print(fit)
  print(seasonal)
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

write_rds(df_list, "results/arima_stepwise_seasonal_fit.rds")

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
  ggsave(paste0("Arima seasonal results for ", predict_cut, ".pdf"), 
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
