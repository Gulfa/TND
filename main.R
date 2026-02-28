# main.R — Entry point for VE Measure Comparison Simulation
# Run via: Rscript /home/gunnar/fhi/TND/main.R
# or:      bash run_in_docker.sh [container_name]

# Support both host path and container path (/fhi/TND)
script_dir <- "/fhi/TND"
if (!dir.exists(script_dir)) script_dir <- "/home/gunnar/fhi/TND"
setwd(script_dir)
dir.create("outputs", showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load model and functions
# ---------------------------------------------------------------------------
message("Loading model...")
source("R/model.R")       # defines odin_sir (compiled at source time)

message("Loading simulation functions...")
source("R/simulate.R")    # run_simulation()

message("Loading VE measure functions...")
source("R/measures.R")    # compute_ve(), compute_ve_statistical()

message("Loading plot functions...")
source("R/plots.R")       # plot_*()

# ---------------------------------------------------------------------------
# Run all scenarios (defines all_results and results_df)
# ---------------------------------------------------------------------------
message("\nRunning scenarios...")
source("scenarios.R")

# ---------------------------------------------------------------------------
# Print summary table
# ---------------------------------------------------------------------------
message("\n========== VE ESTIMATOR COMPARISON ==========")
print(results_df[, c("scenario", "AR_ref", "AR_vax", "VE_RR",
                      "VE_HR_analytical", "VE_HR_empirical", "VE_TND")],
      digits = 3, row.names = FALSE)

# ---------------------------------------------------------------------------
# Statistical analysis: unadjusted vs adjusted for observed risk group
# (restricted to tested individuals; piecewise Poisson counting process for HR)
# ---------------------------------------------------------------------------
message("\n========== STATISTICAL ANALYSIS (tested individuals only) ==========")
stat_df <- do.call(rbind, lapply(all_results, function(res) {
  s <- res$ve_stat
  data.frame(
    scenario         = res$params$name,
    VE_RR_unadj      = s$VE_RR_unadj,
    VE_RR_adj        = s$VE_RR_adj,
    VE_TND_unadj     = s$VE_TND_unadj,
    VE_TND_adj       = s$VE_TND_adj,
    VE_TND_tv_unadj  = s$VE_TND_tv_unadj,
    VE_TND_tv_adj    = s$VE_TND_tv_adj,
    VE_HR_unadj      = s$VE_HR_unadj,
    VE_HR_adj        = s$VE_HR_adj,
    stringsAsFactors = FALSE
  )
}))
print(stat_df, digits = 3, row.names = FALSE)

# ---------------------------------------------------------------------------
# Validation report
# ---------------------------------------------------------------------------
message("\n========== VALIDATION CHECKS ==========")
for (res in all_results) {
  diag <- res$ve$diagnostics
  p    <- res$params
  cat(sprintf(
    "\n[%s]  (n_strat=%d)\n  HR_sd (≈0 when sus uniform within groups): %.4f\n  TND identity check: %s\n  AR_ref: %.3f  AR_vax: %.3f\n",
    p$name, p$n_strat,
    ifelse(is.na(diag$HR_sd), NaN, diag$HR_sd),
    ifelse(is.na(diag$tnd_identity_diff),
           "N/A (heterogeneous seek/mu_oa)",
           sprintf("%.2e", diag$tnd_identity_diff)),
    diag$AR_ref, diag$AR_vax
  ))
}

# ---------------------------------------------------------------------------
# Generate all plots
# ---------------------------------------------------------------------------
message("\n========== GENERATING PLOTS ==========")

for (res in all_results) {
  nm  <- res$params$name
  sim <- res$sim
  ve  <- res$ve

  plot_epidemic(sim, scenario_name = nm)
  plot_tnd_table(ve, scenario_name = nm)
}

# Combined time-series: VE_RR / VE_HR / VE_TND over time, all scenarios
plot_ve_over_time(all_results)

# Test positivity by vaccination status over time
plot_test_positivity(all_results)

# End-point bar comparison across all scenarios
plot_ve_comparison(results_df)

# Statistical VE: adjusted vs unadjusted comparison
plot_hazard_rates(all_results)

plot_ve_statistical(all_results)

# Time-varying VE_HR from piecewise Poisson (adjusted, with 95% CI)
plot_ve_hr_time(all_results)

# Theoretical HRR (ODE) vs observed weekly HRR (Poisson / N-cum.cases)
plot_hrr_comparison(all_results)

message("\nDone. All outputs saved to outputs/")
