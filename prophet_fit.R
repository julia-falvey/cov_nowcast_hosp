pacman::p_load(aTSA, stats, forecast, prophet, urca, tidyverse, ldsr, scales, xtable)
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


fit_prophet <- function(cut_date, data, state) {
  train_data <- data %>% filter(date < cut_date) %>%
    # mutate(cap = max(hosp), floor = 0) %>%
    rename(ds = date, y = hosp_transformed)
  model <- prophet(train_data, growth = "linear", 
                   changepoint.prior.scale = 0.01,
                   yearly.seasonality = "auto", 
                   mcmc.samples = 1000)
  n = n_distinct(data %>% filter(date >= cut_date) %>% pull(date))
  future <- make_future_dataframe(model, periods = n, freq = "week") #%>%
    # mutate(floor = 0, cap = max(train_data$cap))
  forecasted <- predict(model, future) %>%
    dplyr::select(date = ds, fitted = yhat)
  
  train_data <- train_data %>% 
    dplyr::select(date = ds, hosp_transformed = y) %>%
    left_join(forecasted)

  forecasted <- forecasted %>% 
    inner_join(data %>% filter(date >= cut_date)) %>%
    dplyr::select(date, hosp_transformed, fitted) %>%
    mutate(hosp = exp(hosp_transformed) - 1, 
           `Point Forecast` = exp(fitted) - 1)
  
  test_metric <- forecasted %>%
    ungroup() %>% 
    mutate(mse = (hosp-`Point Forecast`)^2) %>%
    summarise(rmse = sqrt(sum(mse))/n()) %>%
    pull(rmse)
  
  train_data <- train_data %>%
    bind_rows(forecasted) %>%
    mutate(state = state, cut_date = cut_date, 
           hosp = exp(hosp_transformed - 1)) %>%
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
    # use box-cox transform on hospital data to avoid negatives in predictions
    # we assume box-cox order 0 (so log() transform, but add 1 to avoid fitting
    # log(0) in our model)
    mutate(hosp_transformed = log(hosp+1))
  for(predict_cut in prediction_horizons){
    rmse_spec = fit_prophet(predict_cut, data, state_nm)
    df_list <- bind_rows(df_list, rmse_spec$data)
    out_tibble <- bind_rows(out_tibble, 
                            tibble_row(state = state_nm, 
                                       cut_date = predict_cut, 
                                       rmse_val = rmse_spec$rmse))
  }
}


write_rds(df_list, "prophet_new_horizons_boxcox.rds")

# Do boxcox and inverse boxcox and move MSE into supplement 
for(predict_cut in prediction_horizons){
  p <- df_list %>%
    filter(cut_date == predict_cut, 
           date >= "2022-07-01") %>%
    mutate(fitted = exp(fitted) - 1,
           hosp = exp(hosp_transformed) - 1) %>%
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
  ggsave(paste0("Prophet Results for ", as.Date(predict_cut), ".pdf"), 
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
