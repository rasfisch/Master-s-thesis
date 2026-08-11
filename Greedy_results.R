
#Greedy algorithms extension (Screening + componentwise boosting)


rm(list = ls())

#Packages
packages <- c(
  "tidyverse",
  "lubridate",
  "sandwich",
  "lmtest",
  "knitr"
)

new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

lapply(packages, library, character.only = TRUE)

set.seed(123)

#User settings
data_path <- "/cloud/project/data/forecast_data_full32_yahoo_volume_period.csv"
baseline_performance_path <- "/cloud/project/output_corrected/baseline_performance_raw_corrected.csv"

output_dir <- "/cloud/project/output_greedy"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

initial_window <- 60
validation_share <- 0.40

#Screening settings
d_grid <- 1:15

#Boosting settings
boost_B_max <- 300
boost_nu <- 0.10
boost_B_grid <- c(1, 5, 10, 25, 50, 100, 200, 300)

#Clark-West HAC lags
cw_lags <- 0

#Load final dataset
df <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

y <- df$r
date_vec <- df$date

x_cols <- setdiff(names(df), c("date", "r"))
X <- as.matrix(df[, x_cols])

n <- nrow(df)
p <- ncol(X)

cat("Observations:", n, "\n")
cat("Predictors:", p, "\n")
cat("Sample:", as.character(min(date_vec)), "to", as.character(max(date_vec)), "\n")
cat("Missing y:", sum(is.na(y)), "\n")
cat("Missing X:", sum(is.na(X)), "\n")

stopifnot(!anyNA(y))
stopifnot(!anyNA(X))
stopifnot(p == 32)


#OOS split
oos_start <- initial_window + 1
oos_index <- oos_start:n

actual_oos <- y[oos_index]
dates_oos <- date_vec[oos_index]

cat("Initial estimation window:", initial_window, "months\n")
cat("First OOS forecast date:", as.character(dates_oos[1]), "\n")
cat("Last OOS forecast date:", as.character(dates_oos[length(dates_oos)]), "\n")
cat("Number of OOS forecasts:", length(oos_index), "\n")


#Helper functions
standardize_train_test <- function(X_train_raw, X_test_raw) {
  
  mu <- apply(X_train_raw, 2, mean, na.rm = TRUE)
  sig <- apply(X_train_raw, 2, sd, na.rm = TRUE)
  
  sig[is.na(sig) | sig == 0] <- 1
  
  X_train <- scale(X_train_raw, center = mu, scale = sig)
  X_test <- scale(matrix(X_test_raw, nrow = 1), center = mu, scale = sig)
  
  list(
    X_train = as.matrix(X_train),
    X_test = as.numeric(X_test),
    mu = mu,
    sig = sig
  )
}

mspe <- function(actual, forecast) {
  mean((actual - forecast)^2, na.rm = TRUE)
}

r2_os <- function(actual, forecast, benchmark) {
  1 - mspe(actual, forecast) / mspe(actual, benchmark)
}

success_ratio <- function(actual, forecast) {
  mean((actual > 0) == (forecast > 0), na.rm = TRUE)
}

clark_west_test <- function(actual, forecast_model, forecast_bench, lags = 0) {
  
  e_b <- actual - forecast_bench
  e_m <- actual - forecast_model
  
  f <- e_b^2 - (e_m^2 - (forecast_bench - forecast_model)^2)
  
  fit <- lm(f ~ 1)
  
  if (lags == 0) {
    se <- sqrt(vcovHC(fit, type = "HC0")[1, 1])
  } else {
    se <- sqrt(NeweyWest(fit, lag = lags, prewhite = FALSE)[1, 1])
  }
  
  stat <- coef(fit)[1] / se
  pval <- 1 - pnorm(stat)
  
  list(
    stat = as.numeric(stat),
    p_value = as.numeric(pval)
  )
}

directional_accuracy_test <- function(actual, forecast) {
  
  a <- as.integer(actual > 0)
  f <- as.integer(forecast > 0)
  
  valid <- complete.cases(a, f)
  a <- a[valid]
  f <- f[valid]
  
  n <- length(a)
  
  p_hat <- mean(a == f)
  p_y <- mean(a == 1)
  p_f <- mean(f == 1)
  
  p_star <- p_y * p_f + (1 - p_y) * (1 - p_f)
  var_p <- p_star * (1 - p_star) / n
  
  if (var_p <= 0 || is.na(var_p)) {
    return(list(stat = NA_real_, p_value = NA_real_))
  }
  
  stat <- (p_hat - p_star) / sqrt(var_p)
  pval <- 1 - pnorm(stat)
  
  list(
    stat = as.numeric(stat),
    p_value = as.numeric(pval)
  )
}

add_stars <- function(value, pval, digits = 3) {
  
  stars <- case_when(
    is.na(pval) ~ "",
    pval <= 0.01 ~ "***",
    pval <= 0.05 ~ "**",
    pval <= 0.10 ~ "*",
    TRUE ~ ""
  )
  
  paste0(sprintf(paste0("%.", digits, "f"), value), stars)
}

safe_lm_selected_predict <- function(y_train, X_train, X_test, selected_idx) {
  
  if (length(selected_idx) == 0) {
    return(mean(y_train, na.rm = TRUE))
  }
  
  Xs_train <- X_train[, selected_idx, drop = FALSE]
  Xs_test <- X_test[selected_idx]
  
  train_df <- data.frame(y = y_train, Xs_train)
  test_df <- data.frame(t(Xs_test))
  names(test_df) <- names(train_df)[-1]
  
  fit <- tryCatch(
    lm(y ~ ., data = train_df),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(mean(y_train, na.rm = TRUE))
  }
  
  pred <- tryCatch(
    predict(fit, newdata = test_df),
    error = function(e) mean(y_train, na.rm = TRUE)
  )
  
  as.numeric(pred)
}


#Screening functions
rank_predictors_by_marginal_fit <- function(y_train, X_train) {
  
  #Since X is standardised, absolute correlation and univariate R2 rankings
  #are equivalent up to monotonic transformation.
  cors <- apply(X_train, 2, function(x) {
    suppressWarnings(cor(y_train, x, use = "complete.obs"))
  })
  
  cors[is.na(cors)] <- 0
  
  order(abs(cors), decreasing = TRUE)
}

fit_screening_fixed_d <- function(y_train, X_train_raw, X_test_raw, d) {
  
  std <- standardize_train_test(X_train_raw, X_test_raw)
  X_train <- std$X_train
  X_test <- std$X_test
  
  ranking <- rank_predictors_by_marginal_fit(y_train, X_train)
  d <- min(d, ncol(X_train), length(ranking))
  
  selected_idx <- ranking[1:d]
  
  forecast <- safe_lm_selected_predict(
    y_train = y_train,
    X_train = X_train,
    X_test = X_test,
    selected_idx = selected_idx
  )
  
  list(
    forecast = forecast,
    selected_idx = selected_idx,
    n_selected = length(selected_idx)
  )
}

fit_screening_validation <- function(y_train, X_train_raw, X_test_raw,
                                     d_grid = 1:15,
                                     validation_share = 0.40) {
  
  n_train <- length(y_train)
  
  val_n <- floor(validation_share * n_train)
  ins_n <- n_train - val_n
  
  if (ins_n < 20 || val_n < 10) {
    stop("Training window too small for validation split.")
  }
  
  X_ins_raw <- X_train_raw[1:ins_n, , drop = FALSE]
  y_ins <- y_train[1:ins_n]
  
  X_val_raw <- X_train_raw[(ins_n + 1):n_train, , drop = FALSE]
  y_val <- y_train[(ins_n + 1):n_train]
  
  #Scale inner train and validation using inner training statistics
  mu_ins <- apply(X_ins_raw, 2, mean, na.rm = TRUE)
  sig_ins <- apply(X_ins_raw, 2, sd, na.rm = TRUE)
  sig_ins[is.na(sig_ins) | sig_ins == 0] <- 1
  
  X_ins <- scale(X_ins_raw, center = mu_ins, scale = sig_ins)
  X_val <- scale(X_val_raw, center = mu_ins, scale = sig_ins)
  
  ranking_inner <- rank_predictors_by_marginal_fit(y_ins, X_ins)
  
  val_mspes <- rep(NA_real_, length(d_grid))
  
  for (g in seq_along(d_grid)) {
    
    d <- min(d_grid[g], ncol(X_ins), length(ranking_inner))
    selected_idx <- ranking_inner[1:d]
    
    val_forecasts <- rep(NA_real_, length(y_val))
    
    fit <- tryCatch(
      lm(y_ins ~ X_ins[, selected_idx, drop = FALSE]),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      val_forecasts <- rep(mean(y_ins), length(y_val))
    } else {
      coefs <- coef(fit)
      coefs[is.na(coefs)] <- 0
      
      val_forecasts <- as.numeric(
        cbind(1, X_val[, selected_idx, drop = FALSE]) %*% coefs
      )
    }
    
    val_mspes[g] <- mean((y_val - val_forecasts)^2)
  }
  
  best_d <- d_grid[which.min(val_mspes)]
  
  #Refit on full recursive training window using best_d
  std_full <- standardize_train_test(X_train_raw, X_test_raw)
  X_train <- std_full$X_train
  X_test <- std_full$X_test
  
  ranking_full <- rank_predictors_by_marginal_fit(y_train, X_train)
  selected_idx <- ranking_full[1:min(best_d, length(ranking_full))]
  
  forecast <- safe_lm_selected_predict(
    y_train = y_train,
    X_train = X_train,
    X_test = X_test,
    selected_idx = selected_idx
  )
  
  list(
    forecast = forecast,
    selected_idx = selected_idx,
    n_selected = length(selected_idx),
    best_d = best_d,
    val_mspes = val_mspes
  )
}

#Componentwise boosting functions
boosting_validation_select_B <- function(y_ins, X_ins, y_val, X_val,
                                         B_max = 300,
                                         nu = 0.10) {
  
  p <- ncol(X_ins)
  xss <- colSums(X_ins^2)
  xss[xss == 0 | is.na(xss)] <- 1
  
  intercept <- mean(y_ins)
  
  fitted_ins <- rep(intercept, length(y_ins))
  fitted_val <- rep(intercept, length(y_val))
  
  beta <- rep(0, p)
  
  val_mspe <- rep(NA_real_, B_max + 1)
  val_mspe[1] <- mean((y_val - fitted_val)^2)
  
  selected_path <- integer(B_max)
  
  for (b in 1:B_max) {
    
    residual <- y_ins - fitted_ins
    
    slopes <- colSums(X_ins * residual) / xss
    slopes[is.na(slopes)] <- 0
    
    improvement <- slopes^2 * xss
    j <- which.max(improvement)
    
    update <- nu * slopes[j]
    
    beta[j] <- beta[j] + update
    
    fitted_ins <- fitted_ins + update * X_ins[, j]
    fitted_val <- fitted_val + update * X_val[, j]
    
    val_mspe[b + 1] <- mean((y_val - fitted_val)^2)
    selected_path[b] <- j
  }
  
  B_star <- which.min(val_mspe) - 1
  
  list(
    B_star = B_star,
    val_mspe = val_mspe,
    selected_path = selected_path
  )
}

boosting_refit_predict <- function(y_train, X_train, X_test,
                                   B,
                                   nu = 0.10) {
  
  p <- ncol(X_train)
  xss <- colSums(X_train^2)
  xss[xss == 0 | is.na(xss)] <- 1
  
  intercept <- mean(y_train)
  
  fitted_train <- rep(intercept, length(y_train))
  forecast <- intercept
  
  beta <- rep(0, p)
  selected_path <- integer(0)
  
  if (B == 0) {
    return(list(
      forecast = forecast,
      beta = beta,
      selected_idx = integer(0),
      n_selected = 0,
      selected_path = selected_path
    ))
  }
  
  for (b in 1:B) {
    
    residual <- y_train - fitted_train
    
    slopes <- colSums(X_train * residual) / xss
    slopes[is.na(slopes)] <- 0
    
    improvement <- slopes^2 * xss
    j <- which.max(improvement)
    
    update <- nu * slopes[j]
    
    beta[j] <- beta[j] + update
    
    fitted_train <- fitted_train + update * X_train[, j]
    forecast <- forecast + update * X_test[j]
    
    selected_path <- c(selected_path, j)
  }
  
  selected_idx <- which(beta != 0)
  
  list(
    forecast = as.numeric(forecast),
    beta = beta,
    selected_idx = selected_idx,
    n_selected = length(selected_idx),
    selected_path = selected_path
  )
}

fit_boosting_validation <- function(y_train, X_train_raw, X_test_raw,
                                    B_max = 300,
                                    nu = 0.10,
                                    validation_share = 0.40) {
  
  n_train <- length(y_train)
  
  val_n <- floor(validation_share * n_train)
  ins_n <- n_train - val_n
  
  if (ins_n < 20 || val_n < 10) {
    stop("Training window too small for validation split.")
  }
  
  X_ins_raw <- X_train_raw[1:ins_n, , drop = FALSE]
  y_ins <- y_train[1:ins_n]
  
  X_val_raw <- X_train_raw[(ins_n + 1):n_train, , drop = FALSE]
  y_val <- y_train[(ins_n + 1):n_train]
  
  #Inner scaling
  mu_ins <- apply(X_ins_raw, 2, mean, na.rm = TRUE)
  sig_ins <- apply(X_ins_raw, 2, sd, na.rm = TRUE)
  sig_ins[is.na(sig_ins) | sig_ins == 0] <- 1
  
  X_ins <- scale(X_ins_raw, center = mu_ins, scale = sig_ins)
  X_val <- scale(X_val_raw, center = mu_ins, scale = sig_ins)
  
  B_selection <- boosting_validation_select_B(
    y_ins = y_ins,
    X_ins = as.matrix(X_ins),
    y_val = y_val,
    X_val = as.matrix(X_val),
    B_max = B_max,
    nu = nu
  )
  
  B_star <- B_selection$B_star
  
  #Full recursive training scaling
  std_full <- standardize_train_test(X_train_raw, X_test_raw)
  X_train <- std_full$X_train
  X_test <- std_full$X_test
  
  final_fit <- boosting_refit_predict(
    y_train = y_train,
    X_train = X_train,
    X_test = X_test,
    B = B_star,
    nu = nu
  )
  
  list(
    forecast = final_fit$forecast,
    selected_idx = final_fit$selected_idx,
    n_selected = final_fit$n_selected,
    B_star = B_star,
    beta = final_fit$beta,
    selected_path = final_fit$selected_path,
    val_mspe = B_selection$val_mspe
  )
}


#Storage
model_names <- c(
  "Historical average",
  paste0("Screening ", 1:15),
  "Screening CV",
  paste0("Boosting ", boost_B_grid),
  "Boosting CV"
)

q <- length(oos_index)

forecasts <- matrix(NA_real_, nrow = q, ncol = length(model_names))
colnames(forecasts) <- model_names

screen_cv_selected <- vector("list", q)
screen_5_selected <- vector("list", q)
screen_10_selected <- vector("list", q)
boost_selected <- vector("list", q)
# Store selected variables for all fixed screening models
screen_fixed_selected <- setNames(
  vector("list", length(d_grid)),
  paste0("Screening_", d_grid)
)

for (d in d_grid) {
  screen_fixed_selected[[paste0("Screening_", d)]] <- vector("list", q)
}

#Store selected variables for all fixed boosting models
boost_fixed_selected <- setNames(
  vector("list", length(boost_B_grid)),
  paste0("Boosting_", boost_B_grid)
)

for (B_fixed in boost_B_grid) {
  boost_fixed_selected[[paste0("Boosting_", B_fixed)]] <- vector("list", q)
}

screen_cv_d <- rep(NA_integer_, q)
boost_B <- rep(NA_integer_, q)

n_selected_screen_cv <- rep(NA_integer_, q)
n_selected_screen_5 <- rep(NA_integer_, q)
n_selected_screen_10 <- rep(NA_integer_, q)
n_selected_boost <- rep(NA_integer_, q)


#Recursive expanding-window forecasting loop
for (ii in seq_along(oos_index)) {
  
  t_idx <- oos_index[ii]
  
  train_idx <- 1:(t_idx - 1)
  
  y_train <- y[train_idx]
  X_train_raw <- X[train_idx, , drop = FALSE]
  X_test_raw <- X[t_idx, , drop = FALSE]
  
  #Historical average
  forecasts[ii, "Historical average"] <- mean(y_train, na.rm = TRUE)
  
  #Screening fixed d = 1,..., 15
  for (d in 1:15) {
    
    screen_fit <- fit_screening_fixed_d(
      y_train = y_train,
      X_train_raw = X_train_raw,
      X_test_raw = X_test_raw,
      d = d
    )
    
    forecasts[ii, paste0("Screening ", d)] <- screen_fit$forecast
    screen_fixed_selected[[paste0("Screening_", d)]][[ii]] <- x_cols[screen_fit$selected_idx]
    
    #Store selected variables for d = 5 and d = 10, since those are used later
    if (d == 5) {
      screen_5_selected[[ii]] <- x_cols[screen_fit$selected_idx]
      n_selected_screen_5[ii] <- screen_fit$n_selected
    }
    
    if (d == 10) {
      screen_10_selected[[ii]] <- x_cols[screen_fit$selected_idx]
      n_selected_screen_10[ii] <- screen_fit$n_selected
    }
  }
  

  #Screening with validation-selected d
  screen_cv_fit <- fit_screening_validation(
    y_train = y_train,
    X_train_raw = X_train_raw,
    X_test_raw = X_test_raw,
    d_grid = d_grid,
    validation_share = validation_share
  )
  
  forecasts[ii, "Screening CV"] <- screen_cv_fit$forecast
  screen_cv_selected[[ii]] <- x_cols[screen_cv_fit$selected_idx]
  screen_cv_d[ii] <- screen_cv_fit$best_d
  n_selected_screen_cv[ii] <- screen_cv_fit$n_selected
  
  #Componentwise boosting with fixed B values
  std_boost <- standardize_train_test(X_train_raw, X_test_raw)
  
  for (B_fixed in boost_B_grid) {
    
    boost_fixed_fit <- boosting_refit_predict(
      y_train = y_train,
      X_train = std_boost$X_train,
      X_test = std_boost$X_test,
      B = B_fixed,
      nu = boost_nu
    )
    
    forecasts[ii, paste0("Boosting ", B_fixed)] <- boost_fixed_fit$forecast
    boost_fixed_selected[[paste0("Boosting_", B_fixed)]][[ii]] <- x_cols[boost_fixed_fit$selected_idx]
  }
  

  #Componentwise boosting with validation-selected B
  boost_fit <- fit_boosting_validation(
    y_train = y_train,
    X_train_raw = X_train_raw,
    X_test_raw = X_test_raw,
    B_max = boost_B_max,
    nu = boost_nu,
    validation_share = validation_share
  )
  
  forecasts[ii, "Boosting CV"] <- boost_fit$forecast
  boost_selected[[ii]] <- x_cols[boost_fit$selected_idx]
  boost_B[ii] <- boost_fit$B_star
  n_selected_boost[ii] <- boost_fit$n_selected
  
  if (ii %% 12 == 0 || ii == q) {
    cat("Completed", ii, "of", q, "forecasts:",
        as.character(dates_oos[ii]), "\n")
  }
}


#Performance table
benchmark <- forecasts[, "Historical average"]

performance_greedy <- map_dfr(model_names, function(m) {
  
  f <- forecasts[, m]
  
  model_mspe <- mspe(actual_oos, f)
  bench_mspe <- mspe(actual_oos, benchmark)
  
  ros <- ifelse(
    m == "Historical average",
    NA_real_,
    100 * r2_os(actual_oos, f, benchmark)
  )
  
  mspe_ratio <- model_mspe / bench_mspe
  sr <- success_ratio(actual_oos, f)
  
  if (m == "Historical average") {
    cw_stat <- NA_real_
    cw_p <- NA_real_
    da_stat <- NA_real_
    da_p <- NA_real_
  } else {
    cw <- clark_west_test(
      actual = actual_oos,
      forecast_model = f,
      forecast_bench = benchmark,
      lags = cw_lags
    )
    
    da <- directional_accuracy_test(actual_oos, f)
    
    cw_stat <- cw$stat
    cw_p <- cw$p_value
    da_stat <- da$stat
    da_p <- da$p_value
  }
  
  tibble(
    Model = m,
    MSPE = model_mspe,
    MSPE_ratio = mspe_ratio,
    R2_OS_percent = ros,
    Success_ratio = sr,
    CW_stat = cw_stat,
    CW_p_value = cw_p,
    DA_stat = da_stat,
    DA_p_value = da_p
  )
})

performance_greedy_pretty <- performance_greedy %>%
  mutate(
    R2_OS_report = ifelse(
      is.na(R2_OS_percent),
      "-",
      map2_chr(R2_OS_percent, CW_p_value, add_stars)
    ),
    Success_report = ifelse(
      Model == "Historical average",
      sprintf("%.3f", Success_ratio),
      map2_chr(Success_ratio, DA_p_value, add_stars)
    ),
    MSPE = round(MSPE, 6),
    MSPE_ratio = round(MSPE_ratio, 3)
  ) %>%
  select(
    Model,
    R2_OS_report,
    Success_report,
    MSPE,
    MSPE_ratio,
    CW_p_value,
    DA_p_value
  )

print(performance_greedy_pretty, n = Inf)

write_csv(
  performance_greedy,
  file.path(output_dir, "greedy_performance_raw.csv")
)

write_csv(
  performance_greedy_pretty,
  file.path(output_dir, "greedy_performance_pretty.csv")
)


#Save forecasts
forecast_greedy_df <- tibble(
  date = dates_oos,
  actual = actual_oos
) %>%
  bind_cols(as_tibble(forecasts))

write_csv(
  forecast_greedy_df,
  file.path(output_dir, "greedy_forecasts.csv")
)


#CSPE data and plot
cspe_df <- tibble(date = dates_oos)

for (m in model_names[model_names != "Historical average"]) {
  
  e_b <- actual_oos - benchmark
  e_m <- actual_oos - forecasts[, m]
  
  cspe_df[[m]] <- cumsum(e_b^2 - e_m^2)
}

cspe_long <- cspe_df %>%
  pivot_longer(
    cols = -date,
    names_to = "Model",
    values_to = "CSPE"
  )

write_csv(
  cspe_long,
  file.path(output_dir, "cspe_greedy_long.csv")
)

cspe_plot <- ggplot(cspe_long, aes(x = date, y = CSPE, color = Model)) +
  geom_line(linewidth = 0.7) +
  labs(
    x = NULL,
    y = "CSPE difference",
    title = "Cumulative forecast performance of greedy methods"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(cspe_plot)

ggsave(
  filename = file.path(output_dir, "cspe_greedy_plot.pdf"),
  plot = cspe_plot,
  width = 10,
  height = 6
)


#Model size and tuning output
greedy_model_sizes <- tibble(
  date = dates_oos,
  screening_cv_d = screen_cv_d,
  boosting_B = boost_B,
  n_selected_screening_cv = n_selected_screen_cv,
  n_selected_screening_5 = n_selected_screen_5,
  n_selected_screening_10 = n_selected_screen_10,
  n_selected_boosting = n_selected_boost
)

write_csv(
  greedy_model_sizes,
  file.path(output_dir, "greedy_model_sizes.csv")
)

print(summary(greedy_model_sizes$screening_cv_d))
print(summary(greedy_model_sizes$boosting_B))
print(summary(greedy_model_sizes$n_selected_boosting))


#Selection frequencies
selection_frequency <- function(selected_list, predictor_names) {
  
  map_dfr(predictor_names, function(v) {
    tibble(
      variable = v,
      frequency = mean(map_lgl(selected_list, ~ v %in% .x), na.rm = TRUE)
    )
  }) %>%
    arrange(desc(frequency))
}

screen_cv_freq <- selection_frequency(screen_cv_selected, x_cols)
screen_5_freq <- selection_frequency(screen_5_selected, x_cols)
screen_10_freq <- selection_frequency(screen_10_selected, x_cols)
boost_freq <- selection_frequency(boost_selected, x_cols)

write_csv(
  screen_cv_freq,
  file.path(output_dir, "screening_cv_selection_frequency.csv")
)

write_csv(
  screen_5_freq,
  file.path(output_dir, "screening_5_selection_frequency.csv")
)

write_csv(
  screen_10_freq,
  file.path(output_dir, "screening_10_selection_frequency.csv")
)

write_csv(
  boost_freq,
  file.path(output_dir, "boosting_selection_frequency.csv")
)

cat("\nScreening CV selection frequency:\n")
print(screen_cv_freq, n = Inf)

cat("\nBoosting selection frequency:\n")
print(boost_freq, n = Inf)

#Selection frequencies for all fixed screening models
screen_fixed_freq_all <- map_dfr(names(screen_fixed_selected), function(model_name) {
  
  selection_frequency(screen_fixed_selected[[model_name]], x_cols) %>%
    mutate(model = model_name) %>%
    select(model, variable, frequency)
})

write_csv(
  screen_fixed_freq_all,
  file.path(output_dir, "screening_all_selection_frequency.csv")
)

for (model_name in names(screen_fixed_selected)) {
  
  freq <- selection_frequency(screen_fixed_selected[[model_name]], x_cols)
  
  write_csv(
    freq,
    file.path(output_dir, paste0(tolower(model_name), "_selection_frequency.csv"))
  )
}


#Selection frequencies for all fixed boosting models
boost_fixed_freq_all <- map_dfr(names(boost_fixed_selected), function(model_name) {
  
  selection_frequency(boost_fixed_selected[[model_name]], x_cols) %>%
    mutate(model = model_name) %>%
    select(model, variable, frequency)
})

write_csv(
  boost_fixed_freq_all,
  file.path(output_dir, "boosting_all_selection_frequency.csv")
)

for (model_name in names(boost_fixed_selected)) {
  
  freq <- selection_frequency(boost_fixed_selected[[model_name]], x_cols)
  
  write_csv(
    freq,
    file.path(output_dir, paste0(tolower(model_name), "_selection_frequency.csv"))
  )
}


#Optional: combine baseline and greedy performance tables
if (file.exists(baseline_performance_path)) {
  
  baseline_raw <- read_csv(baseline_performance_path, show_col_types = FALSE)
  
  #Remove duplicate Historical average from greedy table
  greedy_for_combined <- performance_greedy %>%
    filter(Model != "Historical average")
  
  combined_performance <- bind_rows(
    baseline_raw,
    greedy_for_combined
  )
  
  combined_performance_pretty <- combined_performance %>%
    mutate(
      R2_OS_report = ifelse(
        is.na(R2_OS_percent),
        "-",
        map2_chr(R2_OS_percent, CW_p_value, add_stars)
      ),
      Success_report = ifelse(
        Model == "Historical average",
        sprintf("%.3f", Success_ratio),
        map2_chr(Success_ratio, DA_p_value, add_stars)
      ),
      MSPE = round(MSPE, 6),
      MSPE_ratio = round(MSPE_ratio, 3)
    ) %>%
    select(
      Model,
      R2_OS_report,
      Success_report,
      MSPE,
      MSPE_ratio,
      CW_p_value,
      DA_p_value
    )
  
  print(combined_performance_pretty, n = Inf)
  
  write_csv(
    combined_performance,
    file.path(output_dir, "combined_baseline_greedy_performance_raw.csv")
  )
  
  write_csv(
    combined_performance_pretty,
    file.path(output_dir, "combined_baseline_greedy_performance_pretty.csv")
  )
  
  combined_latex <- combined_performance_pretty %>%
    select(Model, R2_OS_report, Success_report, MSPE_ratio) %>%
    kable(
      format = "latex",
      booktabs = TRUE,
      caption = "Out-of-sample forecasting performance including greedy methods",
      col.names = c(
        "Model",
        "$R^2_{OS}$ (\\%)",
        "Success ratio",
        "MSPE ratio"
      ),
      escape = FALSE
    )
  
  cat(combined_latex)
  
  writeLines(
    combined_latex,
    file.path(output_dir, "combined_baseline_greedy_table.tex")
  )
}

read_csv("/cloud/project/output_greedy/greedy_model_sizes.csv") %>%
  count(screening_cv_d) %>%
  arrange(desc(n))


#Boosting-only CSPE plot
boosting_order <- c(
  "Boosting 1",
  "Boosting 5",
  "Boosting 10",
  "Boosting 25",
  "Boosting 50",
  "Boosting 100",
  "Boosting 200",
  "Boosting 300",
  "Boosting CV"
)

cspe_boosting_long <- cspe_long %>%
  filter(Model %in% boosting_order) %>%
  mutate(Model = factor(Model, levels = boosting_order))

boosting_cspe_plot <- ggplot(cspe_boosting_long, aes(x = date, y = CSPE, color = Model)) +
  geom_line(linewidth = 0.8) +
  labs(
    x = NULL,
    y = "CSPE difference",
    title = "Cumulative forecast performance of boosting models"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(boosting_cspe_plot)

ggsave(
  filename = file.path(output_dir, "cspe_boosting_plot.pdf"),
  plot = boosting_cspe_plot,
  width = 10,
  height = 6
)


