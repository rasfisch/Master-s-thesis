library(quantmod)
library(tidyverse)
library(lubridate)
library(slider)

#Load clean macro dataset
macro_predictors_zhang <- read_csv("/cloud/project/data/macro_predictors_zhang.csv") %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)


#Download free front-month crude oil futures data


getSymbols(
  Symbols = "CL=F",
  src = "yahoo",
  from = "2000-01-01",
  to = "2017-01-31",
  auto.assign = TRUE
)

cl_daily <- tibble(
  date = as.Date(index(`CL=F`)),
  futures_price = as.numeric(Cl(`CL=F`)),
  volume = as.numeric(Vo(`CL=F`))
) %>%
  filter(!is.na(futures_price)) %>%
  arrange(date)

#Check daily range
range(cl_daily$date)
head(cl_daily)
tail(cl_daily)


#Convert daily futures data to monthly


nymex_monthly <- cl_daily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(date = month) %>%
  summarise(
    futures_price = last(futures_price[!is.na(futures_price)]),
    volume = sum(volume, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(date)

glimpse(nymex_monthly)

nymex_monthly <- cl_daily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(date = month) %>%
  summarise(
    futures_price = last(futures_price[!is.na(futures_price)]),
    volume = sum(volume, na.rm = TRUE),
    n_daily_obs = n(),
    .groups = "drop"
  ) %>%
  arrange(date) %>%
  filter(date >= as.Date("2000-09-01"))

write_csv(nymex_monthly, "data/nymex_front_month_yahoo_monthly.csv")



