pacman::p_load(plotly, tidyverse, gridExtra, RColorBrewer)

setwd("~/Desktop/git/cov_nowcast_hosp")

get_baseline_vals <- function(df, algorithm){
  df %>%
    group_by(state, cut_date) %>%
    filter(date >= cut_date) %>%
    mutate(error = `Point Forecast` - hosp, 
           square_error = error^2) %>%
    summarise(mean_square_error = mean(square_error, na.rm = T), 
              rmse = sqrt(mean_square_error), 
              .groups = "keep") %>%
    mutate(algorithm_type = algorithm, 
           model_type = "baseline")
}

get_regressor_vals <- function(df, algorithm){
  df %>%
    group_by(state, cut_date, model_type) %>%
    filter(date >= cut_date) %>%
    mutate(error = `Point Forecast` - hosp, 
           square_error = error^2) %>%
    summarise(mean_square_error = mean(square_error, na.rm = T), 
              rmse = sqrt(mean_square_error), 
              .groups = "keep")  %>%
    group_by(cut_date) %>%
    filter(rmse <= quantile(rmse, 0.95)*100) %>%
    ungroup()  %>%
    mutate(algorithm_type = algorithm)
}

get_median_rmse <- function(df){
  df %>%
    group_by(cut_date, algorithm_type) %>%
    summarise(median_rmse = median(rmse, na.rm = T))
}

get_rmse_change <- function(df){
  df %>%
    group_by(state, algorithm_type) %>%
    mutate(baseline_rmse = rmse[cut_date=="2023-04-27"], 
           pct_change = (rmse - baseline_rmse)/baseline_rmse)
}

heatmap_rmse_change <- function(df){
  df %>%
    group_by(state, algorithm_type, cut_date) %>%
    mutate(baseline_rmse = rmse[model_type=="baseline"], 
           pct_change = (rmse - baseline_rmse)/baseline_rmse, 
           pct_change_color = if_else(pct_change > 1, 1, pct_change))
}
# ARIMA -------------------------------------------------------------------
baseline_arima <- read_rds('results/Univariate/arima_simple_df_list.rds') %>%
  get_baseline_vals(., "ARIMA")

models_arima <- read_rds("results/Multivariate/arimaX_26_no_ww_outliers_seasonal_fit.rds") %>%
  get_regressor_vals(., "ARIMA")

# Prophet -----------------------------------------------------------------
baseline_prophet <- read_rds('results/Univariate/prophet_simple_df_list.rds') %>%
  get_baseline_vals(., "Prophet")

models_prophet <- read_rds("results/Multivariate/prophet_regressors_no_ww_outlier_lagged.rds") %>%
  get_regressor_vals(., "Prophet")

# BSTS --------------------------------------------------------------------
baseline_bsts <- read_rds('results/Univariate/bsts_simple_df_list.rds') %>%
  get_baseline_vals(., "BSTS")

models_bsts <- read_rds("results/Multivariate/bsts_lagged_regressors_no_ww_outlier_model.rds") %>%
  get_regressor_vals(., "BSTS")



# Elbow Plot for Baseline -------------------------------------------------

p <- bind_rows(get_median_rmse(baseline_arima), 
               get_median_rmse(baseline_prophet)) %>%
  bind_rows(get_median_rmse(baseline_bsts)) %>%
  ggplot(aes(x = cut_date, y = median_rmse, group = algorithm_type, 
             color = algorithm_type)) +
  geom_line() +
  theme_minimal() + 
  labs(y = "Median RMSE", x = "Prediction Window", color = "") 

ggsave("results/Median RMSE.pdf", 
       plot = p, width = 16, height = 10, units = "in", bg = "white", 
       dpi = "retina")

p <- bind_rows(baseline_arima, baseline_prophet) %>%
  bind_rows(baseline_bsts) %>%
  ggplot(aes(x = cut_date, y = rmse, group = algorithm_type, 
             color = algorithm_type)) +
  geom_line() + 
  theme_minimal() + 
  labs(y = "RMSE", x = "Prediction Window", 
       color = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(state ~ ., nrow =  6, scales = "free_y")

ggsave("results/RMSE.pdf", 
       plot = p, width = 16, height = 10, units = "in", bg = "white", 
       dpi = "retina")

p <- bind_rows(get_rmse_change(baseline_arima), 
               get_rmse_change(baseline_prophet)) %>%
  bind_rows(get_rmse_change(baseline_bsts)) %>%
  ggplot(aes(x = cut_date, y = pct_change, group = algorithm_type, 
             color = algorithm_type)) +
  geom_line() + 
  theme_minimal() + 
  labs(y = "RMSE Change from Farthest Prediction Window", x = "Prediction Window", 
       color = "") +
  ylim(c(-1,NA)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(state ~ ., nrow =  6)

ggsave("results/RMSE Rel Change from Baseline.pdf", 
       plot = p, width = 16, height = 10, units = "in", bg = "white", 
       dpi = "retina")


# Tables ------------------------------------------------------------------
arima_table <- get_rmse_change(baseline_arima) %>%
  select(state, cut_date, rmse) %>%
  group_by(state) %>%
  mutate(first_val = rmse[cut_date == '2023-04-27'], 
         val_pct_change = if_else(cut_date == '2023-04-27', rmse, (rmse-first_val)/first_val)) %>%
  select(state, cut_date, val_pct_change) %>%
  pivot_wider(names_from = "cut_date", values_from = "val_pct_change")

arima_table %>% ungroup() %>% 
  summarise(mean(`2024-01-27`))

prophet_table <- get_rmse_change(baseline_prophet) %>%
  select(state, cut_date, rmse) %>%
  group_by(state) %>%
  mutate(first_val = rmse[cut_date == '2023-04-27'], 
         val_pct_change = if_else(cut_date == '2023-04-27', rmse, (rmse-first_val)/first_val)) %>%
  select(state, cut_date, val_pct_change) %>%
  pivot_wider(names_from = "cut_date", values_from = "val_pct_change")

prophet_table %>% ungroup() %>% 
  summarise(mean(`2024-01-27`))

bsts_table <- get_rmse_change(baseline_bsts) %>%
  select(state, cut_date, rmse) %>%
  group_by(state) %>%
  mutate(first_val = rmse[cut_date == '2023-04-27'], 
         val_pct_change = if_else(cut_date == '2023-04-27', rmse, (rmse-first_val)/first_val)) %>%
  select(state, cut_date, val_pct_change) %>%
  pivot_wider(names_from = "cut_date", values_from = "val_pct_change")

bsts_table %>% ungroup() %>% 
  summarise(mean(`2024-01-27`))


options(xtable.floating = FALSE)
options(xtable.timestamp = "bsts_table")
print(xtable(bsts_table), include.rownames=FALSE) 

ccf = read_csv('results/ccf.csv') 
print(xtable(ccf %>% 
               group_by(type) %>% 
               summarise(mean = mean(abs_acf[abs_acf != 0]), 
                         n_states = sum(abs_acf!=0)) %>%
               ungroup() %>%
               arrange(desc(mean)) %>%
               select(Covariate = type, `Avg. Cross Correlation` = mean, 
                      `N. States with Data` = n_states), 
      caption = "Avg. Cross Correlation With Hospitalization Time Series Across States", 
      label = "tab:ccf_avg"), include.rownames = F)
# Heatmap -----------------------------------------------------------------

heatmap_df <- bind_rows(heatmap_rmse_change(bind_rows(baseline_arima, models_arima)), 
                        heatmap_rmse_change(bind_rows(baseline_prophet, models_prophet))) %>%
  bind_rows(heatmap_rmse_change(bind_rows(baseline_bsts, models_bsts))) %>%
  filter(model_type != "baseline") %>%
  #mutate(type = paste(model_type, algorithm_type, sep = ": "))
  mutate(type = model_type)

for(algorithm in unique(heatmap_df$algorithm_type)){
  p <- ggplot(aes(x = state, y = type, fill = pct_change), 
              data = heatmap_df %>% filter(algorithm_type == algorithm)) + 
    geom_tile() + 
    scale_y_discrete(labels = function(x) str_wrap(x, width = 12)) +
    scale_fill_distiller("Rel RMSE", palette = "RdBu", limits = c(-1, 1), 
                         na.value = "gray10") +
    theme_minimal() +
    theme(axis.title.y = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 0.45)) +
    facet_wrap(~ cut_date)
  ggsave(paste0("results/Predcitors RMSE Change from Baseline at ", cut, ".pdf"),
         plot = p, width = 16, height = 10, units = "in", bg = "white", 
         dpi = "retina") 
}


# Updates -----------------------------------------------------------------

p_df <- bind_rows(get_median_rmse(models_arima  %>%
                                 mutate(algorithm_type = model_type)) %>% 
                 mutate(algorithm = "ARIMA"), 
               get_median_rmse(models_prophet  %>%
                                 mutate(algorithm_type = model_type)) %>% 
                 mutate(algorithm = "Prophet")) %>%
  bind_rows(get_median_rmse(models_bsts  %>%
              mutate(algorithm_type = model_type)) %>% 
                mutate(algorithm = "BSTS"))

p <- p_df %>%
  ggplot(aes(x = cut_date, y = median_rmse, group = algorithm, 
             color = algorithm)) +
  geom_line() +
  theme_minimal() + 
  facet_wrap(. ~ algorithm_type) +
  labs(y = "Median RMSE", x = "Prediction Window", color = "") 



p <- p_df %>% 
  bind_rows(bind_rows(get_median_rmse(baseline_arima), 
          get_median_rmse(baseline_prophet)) %>%
  bind_rows(get_median_rmse(baseline_bsts)) %>%
  mutate(algorithm = algorithm_type, 
         algorithm_type = "Baseline")) %>% 
  mutate(algorithm_type = factor(algorithm_type, levels = c('Baseline', 'All Regressors', 
                                                            '% ED Visits', '% ED Visits and Test Positivity', 
                                                            'Wastewater % Detection', 'Wastewater % Detection and Concentrations'))) %>%
  ggplot(aes(x = cut_date, y = median_rmse, group = algorithm_type, color = algorithm_type)) + 
  geom_line(linewidth = 1.15) + 
  scale_color_manual(values = c('#9e0142', '#2d004b', '#4393c3', 
                                '#2166ac', '#66bd63', '#006837')) +
  labs(color = "Model",  x = "Prediction Window", y = "Median RMSE") +
  theme_minimal() + 
  theme(text = element_text(size = 11.5), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text = element_text(size = 13)) +
  facet_wrap(.~algorithm)

ggsave("results/Median_Multivariate RMSE.pdf", 
       plot = p, width = 16, height = 10, units = "in", bg = "white", 
       dpi = "retina")
