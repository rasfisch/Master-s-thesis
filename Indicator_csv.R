library(quantmod)
library(tidyverse)
library(lubridate)
library(slider)

macro_predictors_zhang <- read_csv("/cloud/project/data/macro_predictors_zhang.csv") %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

nymex_monthly <- read_csv("/cloud/project/data/nymex_front_month_yahoo_monthly.csv") %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)


#Helper functions for technical indicators


make_ma <- function(x, window) {
  slide_dbl(
    x,
    mean,
    .before = window - 1,
    .complete = TRUE
  )
}

make_ma_signal <- function(price, s, l) {
  ma_s <- make_ma(price, s)
  ma_l <- make_ma(price, l)
  as.numeric(ma_s >= ma_l)
}

make_mom_signal <- function(price, m) {
  as.numeric(price >= lag(price, m))
}

make_obv <- function(price, volume) {
  d <- case_when(
    price - lag(price) >= 0 ~ 1,
    price - lag(price) < 0 ~ -1,
    TRUE ~ 0
  )
  
  cumsum(volume * d)
}

make_obv_signal <- function(obv, s, l) {
  ma_s <- make_ma(obv, s)
  ma_l <- make_ma(obv, l)
  as.numeric(ma_s >= ma_l)
}

#Construct all 18 technical indicators

technical_indicators <- nymex_monthly %>%
  arrange(date) %>%
  mutate(
    obv = make_obv(futures_price, volume),
    
    #Moving-average indicators
    MA_1_9  = make_ma_signal(futures_price, 1, 9),
    MA_1_12 = make_ma_signal(futures_price, 1, 12),
    MA_2_9  = make_ma_signal(futures_price, 2, 9),
    MA_2_12 = make_ma_signal(futures_price, 2, 12),
    MA_3_9  = make_ma_signal(futures_price, 3, 9),
    MA_3_12 = make_ma_signal(futures_price, 3, 12),
    
    #Momentum indicators
    MOM_1  = make_mom_signal(futures_price, 1),
    MOM_2  = make_mom_signal(futures_price, 2),
    MOM_3  = make_mom_signal(futures_price, 3),
    MOM_6  = make_mom_signal(futures_price, 6),
    MOM_9  = make_mom_signal(futures_price, 9),
    MOM_12 = make_mom_signal(futures_price, 12),
    
    #On-balance volume indicators
    OBV_1_9  = make_obv_signal(obv, 1, 9),
    OBV_1_12 = make_obv_signal(obv, 1, 12),
    OBV_2_9  = make_obv_signal(obv, 2, 9),
    OBV_2_12 = make_obv_signal(obv, 2, 12),
    OBV_3_9  = make_obv_signal(obv, 3, 9),
    OBV_3_12 = make_obv_signal(obv, 3, 12)
  ) %>%
  select(
    date,
    MA_1_9, MA_1_12, MA_2_9, MA_2_12, MA_3_9, MA_3_12,
    MOM_1, MOM_2, MOM_3, MOM_6, MOM_9, MOM_12,
    OBV_1_9, OBV_1_12, OBV_2_9, OBV_2_12, OBV_3_9, OBV_3_12
  ) %>%
  arrange(date)

technical_indicators %>%
  summarise(
    n = n(),
    first_date = min(date),
    last_date = max(date),
    across(everything(), ~ sum(is.na(.x)))
  ) %>%
  print(width = Inf)

write_csv(technical_indicators, "data/technical_indicators_yahoo_full18.csv")



