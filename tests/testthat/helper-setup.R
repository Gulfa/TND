# helper-setup.R — sourced automatically by testthat before all test files.
# Sources the model and runs all 6 scenarios once; results are reused across tests.

suppressMessages({
  source("/fhi/TND/R/model.R")
  source("/fhi/TND/R/simulate.R")
  source("/fhi/TND/R/measures.R")

  # scenarios.R defines make_params(), make_flu_scenario(), .calibrate_mu_oa(),
  # and runs: mu_oa_cal, scenarios (list), all_results (list), results_df (df)
  old_wd <- setwd("/fhi/TND")
  source("/fhi/TND/scenarios.R")
  setwd(old_wd)
})
