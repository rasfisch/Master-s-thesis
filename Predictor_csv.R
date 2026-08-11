library(quantmod)
library(tidyverse)
library(lubridate)
library(readxl)
library(janitor)

fred_codes <- c(
  "TB3MS",      #3-month Treasury bill
  "GS10",       #10-year Treasury yield
  "CPIAUCSL",   #CPI
  "USEPUINDXM", #Economic policy uncertainty
  "IGREA",      #Kilian/Dallas Fed global real economic activity
  "M2SL",       #M2 money stock
  "IPG21112N",  #Industrial production: crude oil
  "UNRATE",     #Unemployment rate
  "CFNAI",      #Chicago Fed National Activity Index
  "MCUMFN"      #Capacity utilization: manufacturing
)

getSymbols(fred_codes, src = "FRED")

eia_folder <- "data/eia_raw"

read_eia_monthly <- function(file, value_name) {
  
  path <- file.path("data/eia_raw", file)
  
  raw <- readxl::read_excel(path, sheet = "Data 1", skip = 2)
  
  raw <- raw[, 1:2]
  names(raw) <- c("date", value_name)
  
  raw %>%
    mutate(
      date = as.Date(date),
      date = lubridate::floor_date(date, unit = "month"),
      "{value_name}" := as.numeric(.data[[value_name]])
    ) %>%
    filter(!is.na(date)) %>%
    arrange(date)
}

wti <- read_eia_monthly("RWTCm.xls", "wti")

crude_prod <- read_eia_monthly("MCRFPUS2m.xls", "crude_prod")

crude_stocks_old <- read_eia_monthly(
  "M_EPC0_SAXL_NUS_MBBLm.xls",
  "crude_stocks_old"
)

crude_stocks_new <- read_eia_monthly(
  "MCESTUS1m.xls",
  "crude_stocks_new"
)

crude_stocks <- crude_stocks_old %>%
  full_join(crude_stocks_new, by = "date") %>%
  mutate(
    crude_stocks = if_else(
      date <= as.Date("2015-12-01"),
      crude_stocks_old,
      crude_stocks_new
    )
  ) %>%
  select(date, crude_stocks) %>%
  arrange(date)

crude_imports <- read_eia_monthly("MCRIMUS2m.xls", "crude_imports")

rac <- read_eia_monthly("R1300____3m.xls", "rac")


#Convert FRED xts objects to clean monthly data frames
xts_to_tbl <- function(x, name) {
  tibble(
    date = as.Date(index(x)),
    "{name}" := as.numeric(x[, 1])
  ) %>%
    mutate(date = floor_date(date, "month")) %>%
    arrange(date)
}

fred_data <- list(
  xts_to_tbl(TB3MS, "tbill"),
  xts_to_tbl(GS10, "long_yield"),
  xts_to_tbl(CPIAUCSL, "cpi"),
  xts_to_tbl(USEPUINDXM, "epu"),
  xts_to_tbl(IGREA, "kilian"),
  xts_to_tbl(M2SL, "m2"),
  xts_to_tbl(IPG21112N, "ip_crude"),
  xts_to_tbl(UNRATE, "unrate"),
  xts_to_tbl(CFNAI, "cfnai"),
  xts_to_tbl(MCUMFN, "capacity_util")
) %>%
  reduce(full_join, by = "date") %>%
  arrange(date)

glimpse(fred_data)


#Merge EIA series
eia_data <- list(
  wti,
  crude_prod,
  crude_stocks,
  crude_imports,
  rac
) %>%
  reduce(full_join, by = "date") %>%
  arrange(date)

glimpse(eia_data)


#Merge FRED and EIA data
macro_raw <- fred_data %>%
  full_join(eia_data, by = "date") %>%
  arrange(date)

glimpse(macro_raw)

#Construct Zhang-style macro variables
macro_data <- macro_raw %>%
  arrange(date) %>%
  mutate(
    #Real oil prices
    real_wti = wti / cpi * 100,
    real_rac = rac / cpi * 100,
    
    #Oil returns
    wti_return = log(real_wti) - lag(log(real_wti)),
    rac_return = log(real_rac) - lag(log(real_rac)),
    
    #Macroeconomic predictors
    inflation = log(cpi) - lag(log(cpi)),
    oil_prod_growth = log(crude_prod) - lag(log(crude_prod)),
    oil_stock_growth = log(crude_stocks) - lag(log(crude_stocks)),
    oil_import_growth = log(crude_imports) - lag(log(crude_imports)),
    m2_growth = log(m2) - lag(log(m2)),
    ip_crude_growth = log(ip_crude) - lag(log(ip_crude)),
    unrate_change = unrate - lag(unrate)
  )


#Add publication lag for macro variables
#Zhang et al. lag some macro variables by an additional month
#to account for publication delays.
#
#Financial variables such as interest rates are generally available
#contemporaneously, while macro quantity/activity variables are lagged.

macro_data <- macro_data %>%
  arrange(date) %>%
  mutate(
    tbill_lagged = tbill,
    long_yield_lagged = long_yield,
    
    inflation_lagged = lag(inflation, 1),
    svar_lagged = NA_real_,       # filled later
    epu_lagged = lag(epu, 1),
    kilian_lagged = lag(kilian, 1),
    oil_prod_growth_lagged = lag(oil_prod_growth, 1),
    oil_stock_growth_lagged = lag(oil_stock_growth, 1),
    oil_import_growth_lagged = lag(oil_import_growth, 1),
    m2_growth_lagged = lag(m2_growth, 1),
    ip_crude_growth_lagged = lag(ip_crude_growth, 1),
    unrate_change_lagged = lag(unrate_change, 1),
    cfnai_lagged = lag(cfnai, 1),
    capacity_util_lagged = lag(capacity_util, 1)
  )


#Construct stock return variance
#Zhang uses the sum of squared daily S&P 500 returns within each month.
#This downloads daily S&P 500 from FRED and constructs monthly variance.

getSymbols("^GSPC", src = "yahoo", from = "1985-01-01", to = "2017-01-31")

sp500_daily <- tibble(
  date = as.Date(index(GSPC)),
  sp500 = as.numeric(Ad(GSPC))
) %>%
  filter(!is.na(sp500)) %>%
  arrange(date) %>%
  mutate(
    daily_return = log(sp500) - lag(log(sp500)),
    month = floor_date(date, "month")
  )

svar_monthly <- sp500_daily %>%
  group_by(date = month) %>%
  summarise(
    svar = sum(daily_return^2, na.rm = TRUE),
    .groups = "drop"
  )

macro_data <- macro_data %>%
  select(-any_of(c("svar", "svar_lagged"))) %>%
  left_join(svar_monthly, by = "date") %>%
  arrange(date) %>%
  mutate(
    svar_lagged = svar
  )


#Create final 14 macro predictors
macro_predictors <- macro_data %>%
  transmute(
    date,
    
    # Dependent variables
    wti_return,
    rac_return,
    
    # 14 macro predictors
    tbill = tbill_lagged,
    long_yield = long_yield_lagged,
    inflation = inflation_lagged,
    svar = svar_lagged,
    epu = epu_lagged,
    kilian = kilian_lagged,
    oil_prod_growth = oil_prod_growth_lagged,
    oil_stock_growth = oil_stock_growth_lagged,
    oil_import_growth = oil_import_growth_lagged,
    m2_growth = m2_growth_lagged,
    ip_crude_growth = ip_crude_growth_lagged,
    unrate_change = unrate_change_lagged,
    cfnai = cfnai_lagged,
    capacity_util = capacity_util_lagged
  ) %>%
  arrange(date)

glimpse(macro_predictors)


#Construct dependent variable for forecasting
#Row t contains predictors known at t.
#r is the return from t to t+1.

macro_predictors <- macro_predictors %>%
  mutate(
    r = lead(wti_return),
    r_rac = lead(rac_return)
  )


#Keep Zhang sample
macro_predictors_zhang <- macro_predictors %>%
  filter(
    date >= as.Date("1986-02-01"),
    date <= as.Date("2016-11-01")
  )

# Check missing values
macro_predictors_zhang %>%
  summarise(
    n = n(),
    across(everything(), ~ sum(is.na(.x)))
  ) %>%
  print(width = Inf)

write_csv(macro_predictors_zhang, "data/macro_predictors_zhang.csv")

glimpse(crude_stocks)

range(crude_stocks$date)

crude_stocks %>%
  filter(date >= as.Date("1986-02-01"),
         date <= as.Date("2016-11-01")) %>%
  summarise(
    n = n(),
    missing = sum(is.na(crude_stocks)),
    first = min(date),
    last = max(date)
  )



