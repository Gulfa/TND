# R/simulate.R
# run_simulation(params) → list(ode_out, params, N_strat, T_study)
#
# params fields (all required):
#   n_strat        — number of strata
#   N              — total population
#   gamma, beta    — SIR parameters
#   f              — length-n vector of population fractions (sums to 1)
#   sus            — length-n susceptibility multipliers (1 = reference)
#   inf_mod        — length-n infectiousness multipliers
#   seed           — length-n initial infected counts per stratum
#   p_seek_case    — length-n P(test | target-disease infection) per stratum
#   p_seek_control — length-n P(test | other-illness episode) per stratum
#   mu_oa          — length-n daily other-illness rate per person per stratum
#   is_vaccinated  — length-n logical: which strata are treated as "vaccinated"
#   T_study        — study duration in days

run_simulation <- function(params) {
  N_strat <- params$f * params$N
  n       <- params$n_strat

  model <- odin_sir$new(
    n_strat          = n,
    beta             = params$beta,
    gamma            = params$gamma,
    N                = params$N,
    sus              = params$sus,
    inf_mod          = params$inf_mod,
    N_strat          = N_strat,
    seed             = params$seed,
    p_seek_case      = params$p_seek_case,
    p_seek_control   = params$p_seek_control,
    mu_oa            = params$mu_oa,
    mu_oa_amp        = params$mu_oa_amp,   # now length-n vector
    T_study          = params$T_study
  )

  times   <- seq(0, params$T_study, by = 1)
  ode_out <- as.data.frame(model$run(times))

  # Normalise column names:
  #   odin produces "S[1]" which as.data.frame converts to "S.1." via make.names()
  #   Standardise to "S_1", "I_1", "Rec_1", "CI_1"
  nms <- names(ode_out)
  nms <- gsub("\\[(\\d+)\\]", "_\\1", nms)   # "S[1]"  → "S_1"
  nms <- gsub("\\.(\\d+)\\.",  "_\\1", nms)   # "S.1."  → "S_1"
  # Handle "step" → "t" for any discrete-time fallback
  nms[nms == "step"] <- "t"
  names(ode_out) <- nms

  # Keep only core compartments + time; drop intermediate odin vars (I_wt_*, I_eff, lambda)
  keep <- grepl("^t$|^S_\\d+$|^I_\\d+$|^Rec_\\d+$|^CI_\\d+$|^h_\\d+$|^lambda$|^T_pos_\\d+$|^T_neg_\\d+$", names(ode_out))
  ode_out <- ode_out[, keep, drop = FALSE]

  # Validation: mass balance per stratum (relative tolerance 0.1%)
  tol_rel <- 1e-3
  for (i in seq_len(n)) {
    mb  <- ode_out[[paste0("S_", i)]] +
           ode_out[[paste0("I_", i)]] +
           ode_out[[paste0("Rec_", i)]]
    err <- max(abs(mb - N_strat[i])) / N_strat[i]
    if (err > tol_rel)
      warning(sprintf("Mass balance violation in stratum %d (rel err = %.2e)", i, err))
  }

  list(
    ode_out = ode_out,
    params  = params,
    N_strat = N_strat,
    T_study = params$T_study
  )
}
