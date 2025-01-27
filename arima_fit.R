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

fit_arima <- function(cut_date, data, fit, seasonal, state, freq) {
  train_data <- data %>% filter(date < cut_date)
  ts_dat <- ts(train_data$hosp, frequency = freq)

  model <- tryCatch({Arima(ts_dat, order = fit,
                           method = "ML", transform.pars=FALSE, seasonal = seasonal)}, 
                    error =  function(err){
                      mod <- Arima(ts_dat, order = fit,
                            method = "CSS", seasonal = seasonal)
                      return(mod)
                    })
  
  train_data$fitted <- model$fitted
  n = n_distinct(data %>% filter(date >= cut_date) %>% pull(date))
  forecasted <- as_tibble(forecast(model, h = n)) %>%
    bind_cols(data %>% filter(date >= cut_date)) %>%
    select(date, hosp, everything())
  
  test_metric <- forecasted %>%
    ungroup() %>% 
    mutate(se = (exp(hosp)-exp(as_tibble(forecast(model, h = n))$`Point Forecast`))^2) %>%
    summarise(sse = sum(se), n = n()) %>%
    mutate(rmse = sqrt(sse/n)) %>%
    pull(rmse)
  
  train_data <- train_data %>%
    bind_rows(forecasted) %>%
    mutate(state = state, cut_date = cut_date, 
           fit = list(fit), 
           seasonal = list(seasonal), 
           freq = freq,
           hosp = exp(hosp) - 1, 
           fitted = exp(fitted) - 1, 
           `Point Forecast` = exp(`Point Forecast`) - 1) %>%
    select(state, cut_date, fit, date, hosp, everything())
  return(list(data = train_data, rmse = test_metric))
}

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


freq_period = c(26)
for(freq in freq_period){
  out_tibble <- tribble(~state, ~cut_date, ~rmse_val)
  df_list <- tibble()
  for(state_nm in c(state.abb, "DC")){
    print(state_nm)
    data <- hospitalizations_age %>%
      filter(state == state_nm) %>%
      mutate(hosp = log(total_admissions_all_covid_confirmed+1)) %>%
      dplyr::select(date, hosp) 
    ts_dat <- ts(data$hosp, frequency = freq)
    d = ndiffs(ts_dat)
    # force seasonal difference
    D = max(1, nsdiffs(ts_dat, test = "ch", max.D = 2))
    p = 1:12
    q = 1:12
    P = 0:3
    Q = 0:3
    full_orders = expand_grid(p, d, q)
    full_seasonal = bind_rows(expand_grid(P, D, Q))
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
    fit = c(model_fit$p, model_fit$d, model_fit$q)
    if(length(fit) == 0){
      fit = c(1,1,1)
    }
    seasonal = c(model_fit$P, model_fit$D, model_fit$Q)
    if(length(seasonal) == 0){
      seasonal = c(1, 1, 1)
    }
    rm(model_fit)
    gc()
    for(predict_cut in prediction_horizons){
      print(predict_cut)
      rmse_spec = tryCatch({fit_arima(predict_cut, data, fit, seasonal, state_nm, freq)}, 
                           error =  function(err){
                             fit_arima(predict_cut, data, fit, c(1,0,1), state_nm, freq)
                           })
      df_list <- bind_rows(df_list, rmse_spec$data %>%
                             mutate(hosp = hosp - 1, fitted = fitted - 1, 
                             `Point Forecast` = `Point Forecast` - 1))
      out_tibble <- bind_rows(out_tibble, 
                              tibble_row(state = state_nm, 
                                         cut_date = predict_cut, 
                                         rmse_val = rmse_spec$rmse))
    }
  }
  write_rds(df_list, 
            paste0("results/arima_", freq, "_seasonal_fit.rds"))
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
  
    write_tsv(change_table, 
              paste0("results/change_table_arima_", freq, "wk_seasonablity.tsv"))
}
# main_list <- read_rds('results/arima_12_seasonal_fit.rds') %>%
#   bind_rows(read_rds('results/arima_20_seasonal_fit.rds')) %>%
#   bind_rows(read_rds('results/arima_22_seasonal_fit.rds')) %>%
#   bind_rows(read_rds('results/arima_24_seasonal_fit.rds')) %>%
#   bind_rows(read_rds('results/arima_26_seasonal_fit.rds')) %>%
#   bind_rows(df_list)

for(st_nm in unique(main_list$state)){
  p <- main_list %>%
    select(state:`Point Forecast`, seasonal, freq) %>%
    arrange(cut_date) %>%
    rename(`Fitted Period` = fitted, `Forecast Period` = `Point Forecast`) %>%
    pivot_longer(`Fitted Period`:`Forecast Period`) %>%
    drop_na() %>%
    filter(state == st_nm) %>%
    rowwise() %>%
    mutate(freq = paste0(freq, " - order: ", paste(fit, collapse = ","), 
                         " | season: ", paste(seasonal, collapse = ","))) %>%
    ungroup() %>%
    ggplot(aes(x = date)) + 
    geom_line(aes(y = hosp, color = "Reported \nHospitalizations")) +
    geom_point(aes(y = value, color = name), size = 0.5) + 
    theme_minimal() + 
    scale_color_manual(values = c("#1B9E77", "#7570B3", "#000000")) + 
    labs(y = "Hospitalizations", x = "Date", color = "") +
    facet_grid(cut_date ~ freq, scales = "free_y") + 
    ggtitle(label = paste(st_nm, "plots by period and cut date"))
  ggsave(paste0("results/Arima 52 results for ", st_nm, 
                ".pdf"),
         plot = p, width = 16, height = 10, units = "in", bg = "white", 
         dpi = "retina")
}
# 
# options(xtable.floating = FALSE)
# options(xtable.timestamp = "change_table")
# print(xtable(change_table), include.rownames=FALSE) 


# 
# fits <- main_list %>%
#   filter(freq == 52) %>%
#   select(state, fit, seasonal) %>%
#   distinct()
# refined_tibble <- tribble(~state, ~cut_date, ~rmse_val)
# refined_list <- tibble()
# for(state_nm in c(state.abb, "DC")){
#   print(state_nm)
#   data <- hospitalizations_age %>%
#     filter(state == state_nm) %>%
#     select(date, hosp = total_admissions_all_covid_confirmed) 
#   ts_dat <- ts(data$hosp, frequency = 52)
#   st_fits <- fits %>% filter(state == state_nm)
#   fit = c(st_fits$fit[[1]][1], st_fits$fit[[1]][2], st_fits$fit[[1]][3])
#   seasonal = c(st_fits$seasonal[[1]][1], 
#                max(st_fits$seasonal[[1]][2], 1), #try to force D=1 if D=0 picked
#                st_fits$seasonal[[1]][3])
#   for(predict_cut in prediction_horizons){
#     print(predict_cut)
#     rmse_spec = fit_arima(predict_cut, data, fit, seasonal, state_nm, freq)
#     refined_list <- bind_rows(refined_list, rmse_spec$data)
#     refined_tibble <- bind_rows(refined_tibble, 
#                             tibble_row(state = state_nm, 
#                                        cut_date = predict_cut, 
#                                        rmse_val = rmse_spec$rmse))
#   }
# }

arima_fits <- read_rds("arima_26_seasonal_fit.rds")
