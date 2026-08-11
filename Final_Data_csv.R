library(tidyverse)
library(lubridate)

macro_predictors_zhang <- read_csv("/cloud/project/data/macro_predictors_zhang.csv") %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

technical_indicators <- read_csv("/cloud/project/data/technical_indicators_yahoo_full18.csv") %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

forecast_data_full32 <- macro_predictors_zhang %>%
  left_join(technical_indicators, by = "date") %>%
  select(
    date,
    r,
    
    # 14 macro predictors
    tbill,
    long_yield,
    inflation,
    svar,
    epu,
    kilian,
    oil_prod_growth,
    oil_stock_growth,
    oil_import_growth,
    m2_growth,
    ip_crude_growth,
    unrate_change,
    cfnai,
    capacity_util,
    
    # 18 technical indicators
    MA_1_9, MA_1_12, MA_2_9, MA_2_12, MA_3_9, MA_3_12,
    MOM_1, MOM_2, MOM_3, MOM_6, MOM_9, MOM_12,
    OBV_1_9, OBV_1_12, OBV_2_9, OBV_2_12, OBV_3_9, OBV_3_12
  ) %>%
  drop_na() %>%
  arrange(date)

forecast_data_full32 %>%
  summarise(
    n = n(),
    first_date = min(date),
    last_date = max(date),
    across(everything(), ~ sum(is.na(.x)))
  ) %>%
  print(width = Inf)

ncol(forecast_data_full32)

write_csv(
  forecast_data_full32,
  "/cloud/project/data/forecast_data_full32_yahoo_volume_period.csv"
)
