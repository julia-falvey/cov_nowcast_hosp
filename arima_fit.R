pacman::p_load(aTSA, stats, forecast, prophet)

prediction_horizons <- c('2023-01-01', '2023-06-24', '2024-01-01')

us_data <- hospitalizations_age %>%
  filter(state == "USA") %>%
  select(date, hosp = total_admissions_all_covid_confirmed) %>%
  filter(date <= '2024-05-01')

us_fit <- auto.arima(us_data$hosp, trace = T, max.p = 20, max.q = 20, max.D = 3, max.d = 3, max.order = 45, 
                     stepwise = F, approximation = T, allowdrift = F)


ggplot(data = bind_cols(us_fit$fitted, us_fit$x) %>%
         mutate(row = row_number()), aes(x = row)) + 
  geom_line(aes(y = `...2`)) + geom_point(aes(y = `...1`)) + 
  labs(title = "Full US fitted (line) vs actual (dots) Auto-Arima")

rmse_func <- function(x, y){
  sqrt((as.vector(x)-as.vector(y))^2)
}
fit_arima <- function(cut_date, data) {
  train_data <- data %>% filter(date <= cut_date)
  model <- Arima(as.vector(train_data$hosp), order=c(2,1,4), include.drift = FALSE )
  test_data <- data %>% 
    filter(date >= cut_date) %>% 
    mutate(n = n(), 
           forecast = as_tibble(forecast(model, h = n))$`Point Forecast`, 
           mse = (hosp-forecast)^2) %>%
    ungroup() %>% 
    summarise(rmse = sqrt(sum(mse))) %>%
    pull(rmse)

  return(test_data)  
}

out_tibble <- tribble(~state, ~cut_date, ~rmse_val)

for(state_nm in unique(hospitalizations_age$state)){
  print(state_nm)
  for(predict_cut in prediction_horizons){
    data <- hospitalizations_age %>%
      filter(state == state_nm) %>%
      select(date, hosp = total_admissions_all_covid_confirmed) %>%
      filter(date <= '2024-05-01')
    rmse_spec = fit_arima(predict_cut, data)
    out_tibble <- bind_rows(out_tibble, 
                            tibble_row(state = state_nm, 
                                       cut_date = predict_cut, 
                                       rmse_val = rmse_spec))
  }
}

View(out_tibble %>%
       filter(state != 'USA') %>%
       pivot_wider(names_from = "cut_date", values_from = "rmse_val"))  

