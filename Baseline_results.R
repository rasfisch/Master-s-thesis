
#Corrected baseline forecasting script
#Dataset: forecast_data_full32_yahoo_volume_period.csv

rm(list = ls())


#Packages
packages <- c(
  "tidyverse",
  "lubridate",
  "glmnet",
  "pls",
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
output_dir <- "/cloud/project/output_corrected"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#Since the usable sample starts in 2001-09, Zhang's original
#1986-2000 estimation window cannot be used.
#Baseline: 60-month expanding initial estimation window.
initial_window <- 60

#Clark-West HAC lags.
#For one-month-ahead non-overlapping forecasts, 0 is defensible.
cw_lags <- 0


#Load final dataset
df <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

glimpse(df)

#Important:
#r is already defined as lead(wti_return) in the data construction step.
#Therefore, row t contains predictors observed at month t, while r is the
#crude oil return from t to t+1.

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
if (initial_window >= n) {
  stop("Initial window is too large relative to sample size.")
}

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
  mean(sign(actual) == sign(forecast), na.rm = TRUE)
}

safe_lm_predict <- function(y_train, X_train, X_test) {
  
  train_df <- data.frame(y = y_train, X_train)
  test_df <- data.frame(t(X_test))
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

#Clark-West test for nested forecast comparison
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

#Directional accuracy z-test.
#This is a simple sign-prediction test, not a full Pesaran-Timmermann statistic.
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


#Corrected penalised regression validation
#This function receives raw X values.
#Inner validation scaling is based only on the inner training sample.
#After lambda/alpha selection, the model is refit on the full recursive
#Training window using scaling is based only on the full recursive training window.

fit_glmnet_validation_raw <- function(y_train, X_train_raw, X_test_raw,
                                      alpha_grid,
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
  
  #Scale inner train and validation using inner train statistics only
  mu_ins <- apply(X_ins_raw, 2, mean, na.rm = TRUE)
  sig_ins <- apply(X_ins_raw, 2, sd, na.rm = TRUE)
  sig_ins[is.na(sig_ins) | sig_ins == 0] <- 1
  
  X_ins <- scale(X_ins_raw, center = mu_ins, scale = sig_ins)
  X_val <- scale(X_val_raw, center = mu_ins, scale = sig_ins)
  
  best_mspe <- Inf
  best_alpha <- NA_real_
  best_lambda <- NA_real_
  
  for (a in alpha_grid) {
    
    fit_path_initial <- glmnet(
      x = X_ins,
      y = y_ins,
      alpha = a,
      standardize = FALSE,
      intercept = TRUE
    )
    
    lambda_grid <- fit_path_initial$lambda
    
    fit_path <- glmnet(
      x = X_ins,
      y = y_ins,
      alpha = a,
      lambda = lambda_grid,
      standardize = FALSE,
      intercept = TRUE
    )
    
    pred_val <- predict(fit_path, newx = X_val)
    
    y_val_mat <- matrix(
      y_val,
      nrow = length(y_val),
      ncol = ncol(pred_val)
    )
    
    mspe_vals <- colMeans((y_val_mat - pred_val)^2)
    idx <- which.min(mspe_vals)
    
    if (mspe_vals[idx] < best_mspe) {
      best_mspe <- mspe_vals[idx]
      best_alpha <- a
      best_lambda <- lambda_grid[idx]
    }
  }
  
  #Refit on full recursive training window
  mu_full <- apply(X_train_raw, 2, mean, na.rm = TRUE)
  sig_full <- apply(X_train_raw, 2, sd, na.rm = TRUE)
  sig_full[is.na(sig_full) | sig_full == 0] <- 1
  
  X_train <- scale(X_train_raw, center = mu_full, scale = sig_full)
  X_test <- scale(matrix(X_test_raw, nrow = 1), center = mu_full, scale = sig_full)
  
  final_fit <- glmnet(
    x = X_train,
    y = y_train,
    alpha = best_alpha,
    lambda = best_lambda,
    standardize = FALSE,
    intercept = TRUE
  )
  
  forecast <- as.numeric(
    predict(final_fit, newx = X_test)
  )
  
  coef_mat <- as.matrix(coef(final_fit))
  selected <- rownames(coef_mat)[coef_mat[, 1] != 0]
  selected <- setdiff(selected, "(Intercept)")
  
  list(
    forecast = forecast,
    alpha = best_alpha,
    lambda = best_lambda,
    selected = selected,
    n_selected = length(selected)
  )
}


#PCR
fit_pcr_adjR2_capped <- function(y_train, X_train, X_test, max_k_cap = 5) {
  
  n_train <- length(y_train)
  max_k <- min(max_k_cap, ncol(X_train), n_train - 5)
  
  pc <- prcomp(X_train, center = FALSE, scale. = FALSE)
  
  scores <- pc$x[, 1:max_k, drop = FALSE]
  test_scores <- as.numeric(
    matrix(X_test, nrow = 1) %*% pc$rotation[, 1:max_k, drop = FALSE]
  )
  
  adj_r2 <- rep(NA_real_, max_k)
  
  for (k in 1:max_k) {
    fit <- lm(y_train ~ scores[, 1:k, drop = FALSE])
    adj_r2[k] <- summary(fit)$adj.r.squared
  }
  
  best_k <- which.max(adj_r2)
  
  final_scores <- scores[, 1:best_k, drop = FALSE]
  final_fit <- lm(y_train ~ final_scores)
  
  forecast <- as.numeric(
    cbind(1, matrix(test_scores[1:best_k], nrow = 1)) %*% coef(final_fit)
  )
  
  list(
    forecast = forecast,
    k = best_k
  )
}




#PLS
fit_pls_onefactor <- function(y_train, X_train, X_test) {
  
  train_df <- data.frame(y = y_train, X_train)
  test_df <- data.frame(t(X_test))
  names(test_df) <- names(train_df)[-1]
  
  fit <- tryCatch(
    plsr(y ~ ., data = train_df, ncomp = 1, scale = FALSE, validation = "none"),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(forecast = mean(y_train, na.rm = TRUE)))
  }
  
  forecast <- tryCatch(
    as.numeric(predict(fit, newdata = test_df, ncomp = 1)),
    error = function(e) mean(y_train, na.rm = TRUE)
  )
  
  list(forecast = forecast)
}


#Storage
model_names <- c(
  "Historical average",
  "Elastic Net",
  "Lasso",
  "Ridge",
  "PC",
  "PLS",
  "Mean",
  "Median",
  "Trimmed mean",
  "DMSPE(1)",
  "DMSPE(0.9)"
)

q <- length(oos_index)

forecasts <- matrix(NA_real_, nrow = q, ncol = length(model_names))
colnames(forecasts) <- model_names

individual_forecasts_store <- matrix(NA_real_, nrow = q, ncol = p)
colnames(individual_forecasts_store) <- x_cols

selected_lasso <- vector("list", q)
selected_en <- vector("list", q)

n_selected_lasso <- rep(NA_integer_, q)
n_selected_en <- rep(NA_integer_, q)

alpha_en_store <- rep(NA_real_, q)
lambda_en_store <- rep(NA_real_, q)
lambda_lasso_store <- rep(NA_real_, q)

pc_k_store <- rep(NA_integer_, q)


#Recursive expanding-window forecasting loop
for (ii in seq_along(oos_index)) {
  
  t_idx <- oos_index[ii]
  
  train_idx <- 1:(t_idx - 1)
  
  y_train <- y[train_idx]
  X_train_raw <- X[train_idx, , drop = FALSE]
  X_test_raw <- X[t_idx, , drop = FALSE]
  
  #Standardised version for PCR, PLS and univariate regressions
  std <- standardize_train_test(X_train_raw, X_test_raw)
  X_train <- std$X_train
  X_test <- std$X_test
  
  #Historical average
  forecasts[ii, "Historical average"] <- mean(y_train, na.rm = TRUE)
  
  #LASSO
  lasso_fit <- fit_glmnet_validation_raw(
    y_train = y_train,
    X_train_raw = X_train_raw,
    X_test_raw = X_test_raw,
    alpha_grid = 1,
    validation_share = 0.40
  )
  
  forecasts[ii, "Lasso"] <- lasso_fit$forecast
  selected_lasso[[ii]] <- lasso_fit$selected
  n_selected_lasso[ii] <- lasso_fit$n_selected
  lambda_lasso_store[ii] <- lasso_fit$lambda
  
  #Elastic Net
  en_fit <- fit_glmnet_validation_raw(
    y_train = y_train,
    X_train_raw = X_train_raw,
    X_test_raw = X_test_raw,
    alpha_grid = seq(0.1, 0.9, by = 0.1),
    validation_share = 0.40
  )
  
  forecasts[ii, "Elastic Net"] <- en_fit$forecast
  selected_en[[ii]] <- en_fit$selected
  n_selected_en[ii] <- en_fit$n_selected
  alpha_en_store[ii] <- en_fit$alpha
  lambda_en_store[ii] <- en_fit$lambda
  
  #Ridge
  ridge_fit <- fit_glmnet_validation_raw(
    y_train = y_train,
    X_train_raw = X_train_raw,
    X_test_raw = X_test_raw,
    alpha_grid = 0,
    validation_share = 0.40
  )
  
  forecasts[ii, "Ridge"] <- ridge_fit$forecast
  
  #PCR
  pcr_fit <- fit_pcr_adjR2_capped(
    y_train = y_train,
    X_train = X_train,
    X_test = X_test,
    max_k_cap = 5
  )
  
  forecasts[ii, "PC"] <- pcr_fit$forecast
  pc_k_store[ii] <- pcr_fit$k
  
  #PLS
  pls_fit <- fit_pls_onefactor(
    y_train = y_train,
    X_train = X_train,
    X_test = X_test
  )
  
  forecasts[ii, "PLS"] <- pls_fit$forecast
  
  #Individual univariate forecasts
  individual_forecasts <- rep(NA_real_, p)
  
  for (j in 1:p) {
    
    Xj_train <- matrix(X_train[, j], ncol = 1)
    Xj_test <- X_test[j]
    
    individual_forecasts[j] <- safe_lm_predict(
      y_train = y_train,
      X_train = Xj_train,
      X_test = Xj_test
    )
  }
  
  individual_forecasts_store[ii, ] <- individual_forecasts
  
  #Forecast combinations
  forecasts[ii, "Mean"] <- mean(individual_forecasts, na.rm = TRUE)
  forecasts[ii, "Median"] <- median(individual_forecasts, na.rm = TRUE)
  
  valid_ind <- sort(individual_forecasts[!is.na(individual_forecasts)])
  
  if (length(valid_ind) > 2) {
    forecasts[ii, "Trimmed mean"] <- mean(valid_ind[2:(length(valid_ind) - 1)])
  } else {
    forecasts[ii, "Trimmed mean"] <- mean(valid_ind)
  }
  
  if (ii <= 2) {
    
    forecasts[ii, "DMSPE(1)"] <- mean(individual_forecasts, na.rm = TRUE)
    forecasts[ii, "DMSPE(0.9)"] <- mean(individual_forecasts, na.rm = TRUE)
    
  } else {
    
    past_actual <- actual_oos[1:(ii - 1)]
    past_individual <- individual_forecasts_store[1:(ii - 1), , drop = FALSE]
    
    dmspe_forecast <- function(theta) {
      
      errors <- sweep(past_individual, 1, past_actual, "-")
      squared_errors <- errors^2
      
      age <- rev(seq_along(past_actual)) - 1
      time_weights <- theta^age
      
      phi <- colSums(
        sweep(squared_errors, 1, time_weights, "*"),
        na.rm = TRUE
      )
      
      eps <- 1e-12
      
      if (any(phi == 0, na.rm = TRUE)) {
        w <- as.numeric(phi == 0)
        w <- w / sum(w)
      } else {
        inv_phi <- 1 / (phi + eps)
        inv_phi[!is.finite(inv_phi)] <- 0
        
        if (sum(inv_phi) == 0) {
          w <- rep(1 / length(inv_phi), length(inv_phi))
        } else {
          w <- inv_phi / sum(inv_phi)
        }
      }
      
      sum(w * individual_forecasts)
    }
    
    forecasts[ii, "DMSPE(1)"] <- dmspe_forecast(theta = 1.0)
    forecasts[ii, "DMSPE(0.9)"] <- dmspe_forecast(theta = 0.9)
  }
  
  if (ii %% 12 == 0 || ii == q) {
    cat("Completed", ii, "of", q, "forecasts:",
        as.character(dates_oos[ii]), "\n")
  }
}


#Performance table
benchmark <- forecasts[, "Historical average"]

performance <- map_dfr(model_names, function(m) {
  
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
  } else {
    cw <- clark_west_test(
      actual = actual_oos,
      forecast_model = f,
      forecast_bench = benchmark,
      lags = cw_lags
    )
    cw_stat <- cw$stat
    cw_p <- cw$p_value
  }
  
  if (m == "Historical average") {
    da_stat <- NA_real_
    da_p <- NA_real_
  } else {
    da <- directional_accuracy_test(actual_oos, f)
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

performance_pretty <- performance %>%
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

print(performance_pretty, n = Inf)

#Save tables
write_csv(
  performance,
  file.path(output_dir, "baseline_performance_raw_corrected.csv")
)

write_csv(
  performance_pretty,
  file.path(output_dir, "baseline_performance_pretty_corrected.csv")
)

baseline_latex <- performance_pretty %>%
  select(Model, R2_OS_report, Success_report, MSPE_ratio) %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    caption = "Baseline out-of-sample forecasting performance",
    col.names = c(
      "Model",
      "$R^2_{OS}$ (\\%)",
      "Success ratio",
      "MSPE ratio"
    ),
    escape = FALSE
  )

cat(baseline_latex)

writeLines(
  baseline_latex,
  file.path(output_dir, "baseline_performance_table_corrected.tex")
)


#Save forecasts
forecast_df <- tibble(
  date = dates_oos,
  actual = actual_oos
) %>%
  bind_cols(as_tibble(forecasts))

write_csv(
  forecast_df,
  file.path(output_dir, "baseline_forecasts_corrected.csv")
)

#CSPE plot
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
  file.path(output_dir, "cspe_baseline_long_corrected.csv")
)

cspe_plot <- ggplot(cspe_long, aes(x = date, y = CSPE, color = Model)) +
  geom_line(linewidth = 0.7) +
  labs(
    x = NULL,
    y = "CSPE difference",
    title = "Cumulative forecast performance relative to historical average"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(cspe_plot)

ggsave(
  filename = file.path(output_dir, "cspe_baseline_plot_corrected.pdf"),
  plot = cspe_plot,
  width = 10,
  height = 6
)


#Variable selection output
selection_df <- tibble(
  date = dates_oos,
  n_selected_lasso = n_selected_lasso,
  n_selected_en = n_selected_en,
  alpha_en = alpha_en_store,
  lambda_en = lambda_en_store,
  lambda_lasso = lambda_lasso_store,
  pc_k = pc_k_store
)

write_csv(
  selection_df,
  file.path(output_dir, "selection_model_sizes_corrected.csv")
)

selection_frequency <- function(selected_list, predictor_names) {
  
  map_dfr(predictor_names, function(v) {
    tibble(
      variable = v,
      frequency = mean(map_lgl(selected_list, ~ v %in% .x), na.rm = TRUE)
    )
  }) %>%
    arrange(desc(frequency))
}

lasso_selection_freq <- selection_frequency(selected_lasso, x_cols)
en_selection_freq <- selection_frequency(selected_en, x_cols)

write_csv(
  lasso_selection_freq,
  file.path(output_dir, "lasso_selection_frequency_corrected.csv")
)

write_csv(
  en_selection_freq,
  file.path(output_dir, "en_selection_frequency_corrected.csv")
)

print(lasso_selection_freq, n = Inf)
print(en_selection_freq, n = Inf)

