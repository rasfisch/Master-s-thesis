High-Dimensional Crude Oil Return Forecasting

Overview

This repository contains an econometric forecasting framework developed as part of my master's thesis in Mathematics-Economics at Aalborg University.

The project investigates whether high-dimensional variable selection methods can improve forecasts of monthly crude oil returns using macroeconomic predictors and technical indicators.

The original research implements penalised regression methods, dimension reduction, forecast combinations, screening, and componentwise boosting within a recursive out-of-sample evaluation framework.

The strongest baseline model, LASSO, achieved an out-of-sample R² of 14.23%, while the best-performing boosting specification improved this to 18.13%.

Future Development

The original thesis research is complete. The next development phase aims to extend the existing codebase into a more modular and reusable framework for financial time-series forecasting.

Planned Extensions

1. Hybrid Forecasting Models

Combine screening with LASSO, Elastic Net, and componentwise boosting.

Evaluate whether preliminary variable selection improves predictive accuracy.

2. Improved Hyperparameter Optimization

Implement alternative time-series cross-validation procedures.

Introduce joint optimization of boosting iterations and learning rates.

3. Stability Selection

Investigate predictor-selection stability across time-series subsamples.

Identify variables that remain informative across different estimation periods.

4. Extended Dataset and Robustness Analysis

Extend the forecasting sample beyond 2016.

Evaluate performance during COVID-19, the 2022 energy crisis, and more recent market conditions.

5. Regime-Dependent Forecasting

Investigate whether different models perform better under different market conditions.

Explore volatility-dependent and time-varying forecasting approaches.

6. Economic Performance Evaluation

Develop a backtesting framework to assess the economic value of forecasts.

Incorporate transaction costs, position sizing, and risk-adjusted performance measures.

Software Development

Alongside these research extensions, the codebase will be progressively refactored into reusable components, with improved documentation, reproducibility, and automated testing.

Development tasks and implementation progress will be tracked through GitHub Issues.

**Instructions**

Run the scripts in the following order:

1. Predictor_csv.R
| Downloads the FRED/Yahoo series, reads the EIA raw files, 
constructs the 14 macro predictors and forward return, 
and writes macro_predictors_zhang.csv

2. Nymex_csv.R
| downloads Yahoo CL=F,
 aggregates daily futures prices/volume to monthly frequency,
 and creates the monthly futures CSV

4. Indicator_csv.R
| Constructs all 18 MA, momentum and OBV indicators.
   
5. Final_Data_csv.R
| Merges the 14 macro predictors and 18 technical indicators into the final 32-predictor forecasting dataset.

6. Baseline_results.R
| Runs the entire baseline forecasting exercise:
LASSO, Elastic Net, ridge, PCR, PLS, combinations,
performance tables, CSPE, selection frequencies, etc.
It explicitly reads the final 32-predictor dataset.
  
7. Greedy_results.R
| Runs screening and componentwise boosting,
produces the greedy tables, frequencies,
full CSPE and boosting-only CSPE.

PCR_diagnostic.R contains an additional diagnostic for the number of
principal components selected in the baseline PCR model.
