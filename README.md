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
