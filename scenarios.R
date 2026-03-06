# scenarios.R
# Influenza-like VE comparison scenarios.
# gamma = 1/2 (2-day gen time), T_study = 90 days.
#
# Base structure: 2 observed risk groups (low / high susceptibility, ratio 2)
# within each vaccination status.  An optional layer of UNOBSERVED within-group
# heterogeneity splits each risk×vax cell into two sub-strata (A and B), giving
# 4 strata per vaccination status (8 total) when enabled.
#
# Susceptibilities are normalised so that weighted-mean(sus_ref) = 1.0,
# keeping R_eff constant regardless of sus_ratio or risk-group shares.

# ---------------------------------------------------------------------------
# make_params(): general n-stratum parameter builder
# ---------------------------------------------------------------------------
make_params <- function(
    name,
    R0,
    f,               # population fractions (length n, must sum to 1)
    sus,             # susceptibility multipliers
    inf_mod,         # infectiousness multipliers
    seed,            # initial infected per stratum
    p_seek_case,     # P(test | target-disease infection) per stratum
    p_seek_control,  # P(test | other-illness episode) per stratum
    mu_oa,           # daily other-illness rate per person per stratum
    is_vaccinated,   # logical: reference (FALSE) vs vaccinated (TRUE)
    strat_names = NULL,
    risk_obs    = NULL,  # observed risk label per stratum ("low"/"high"); used for adjusted analysis
    N       = 100000,
    gamma     = 1 / 2,   # influenza: ~2-day generation time
    T_study   = 150,     # influenza: one flu season
    mu_oa_amp = 0.5      # seasonal amplitude: background peaks at 1.5× mid-season
) {
  n <- length(f)
  mu_oa_amp <- if (length(mu_oa_amp) == 1) rep(mu_oa_amp, n) else mu_oa_amp
  stopifnot(
    length(sus)            == n,
    length(inf_mod)        == n,
    length(seed)           == n,
    length(p_seek_case)    == n,
    length(p_seek_control) == n,
    length(mu_oa)          == n,
    length(mu_oa_amp)      == n,
    length(is_vaccinated)  == n,
    abs(sum(f) - 1) < 1e-10,
    any(!is_vaccinated),
    any( is_vaccinated),
    all(seed >= 0),
    all(seed <= f * N + 1e-6),         # seed cannot exceed stratum size
    all(p_seek_case    >= 0 & p_seek_case    <= 1),
    all(p_seek_control >= 0 & p_seek_control <= 1),
    all(mu_oa          >= 0),
    all(mu_oa_amp      >= 0 & mu_oa_amp      <= 1)  # amp > 1 → negative rates
  )
  list(
    name           = name,
    R0             = R0,
    n_strat        = n,
    N              = N,
    gamma          = gamma,
    beta           = R0 * gamma,
    f              = f,
    sus            = sus,
    inf_mod        = inf_mod,
    seed           = seed,
    p_seek_case    = p_seek_case,
    p_seek_control = p_seek_control,
    mu_oa          = mu_oa,
    is_vaccinated  = is_vaccinated,
    strat_names    = if (!is.null(strat_names)) strat_names
                     else paste0("Stratum ", seq_len(n)),
    risk_obs       = if (!is.null(risk_obs)) risk_obs
                     else rep("unknown", n),
    T_study        = T_study,
    mu_oa_amp      = mu_oa_amp   # already expanded to length-n above
  )
}

# ---------------------------------------------------------------------------
# make_flu_scenario()
#
# Builds a 4-strat (within_delta = 0) or 8-strat (within_delta > 0) scenario.
#
# Groups (before within-cell split):
#   [Unvax low-risk, Unvax high-risk, Vax low-risk, Vax high-risk]
#
# Susceptibility normalisation:
#   sus_low  = 1 / (1 + f_highrisk_unvax * (sus_ratio - 1))
#   sus_high = sus_low * sus_ratio
#   => weighted-mean(sus_ref) = 1.0 for any sus_ratio / f_highrisk_unvax
#
# Within-cell heterogeneity (when within_delta > 0):
#   sub-A sus = group_sus * (1 - within_delta)   [less susceptible half]
#   sub-B sus = group_sus * (1 + within_delta)   [more susceptible half]
#   Factor between sub-A and sub-B = (1+d)/(1-d);
#   e.g. within_delta = 1/3 gives exactly a factor of 2.
# ---------------------------------------------------------------------------
make_flu_scenario <- function(
    name, R0,
    f_v               = 0.2,
    f_highrisk_unvax  = 0.5,    # P(high-risk | unvax)
    f_highrisk_vax    = 0.5,    # P(high-risk | vax); differs from unvax → indication bias
    VE_s              = 0.6,
    sus_ratio         = 2.0,    # sus_high / sus_low (observed risk-group ratio)
    within_delta      = 0.0,    # within-cell sus spread: 0 = off, 1/3 = factor-2 sub-groups
    p_seek_case_unvax = 0.2,
    p_seek_case_vax   = 0.2,
    p_seek_control    = 0.4,    # P(test | other illness) for vaccinated
    p_seek_control_unvax = p_seek_control,  # defaults to same as vax; set proportionally for care-seeking bias
    mu_oa             = mu_oa_cal,
    mu_oa_ratio       = 1.0,    # mu_oa_high / mu_oa_low; >1 → high-risk have more other illness
    VE_oa             = 0.0,    # vaccine effect on other-illness rate (0 = none, >0 → fewer vax controls)
    mu_oa_amp_ref     = 0.5,    # seasonal amplitude for reference group
    mu_oa_amp_vax     = 0.5,    # seasonal amplitude for vaccinated group
    ...
) {
  stopifnot(VE_s >= 0, VE_s <= 1,
            VE_oa >= 0, VE_oa <= 1,
            sus_ratio > 0,
            f_highrisk_unvax >= 0, f_highrisk_unvax <= 1,
            f_highrisk_vax   >= 0, f_highrisk_vax   <= 1,
            f_v >= 0, f_v <= 1)

  # Normalise to weighted-mean(sus_ref) = 1
  sus_low  <- 1 / (1 + f_highrisk_unvax * (sus_ratio - 1))
  sus_high <- sus_low * sus_ratio

  # Per-group susceptibilities and fractions
  sus_g    <- c(sus_low,         sus_high,
                sus_low*(1-VE_s), sus_high*(1-VE_s))
  f_g      <- c((1-f_v)*(1-f_highrisk_unvax), (1-f_v)*f_highrisk_unvax,
                  f_v *(1-f_highrisk_vax),       f_v *f_highrisk_vax)
  p_case_g <- c(p_seek_case_unvax, p_seek_case_unvax,
                p_seek_case_vax,   p_seek_case_vax)
  p_ctrl_g <- c(p_seek_control_unvax, p_seek_control_unvax,
                p_seek_control,       p_seek_control)
  # Other-illness rate: low-risk = mu_oa, high-risk = mu_oa * mu_oa_ratio;
  # vaccinated strata reduced by (1 - VE_oa) if VE_oa > 0
  mu_oa_g      <- c(mu_oa,           mu_oa * mu_oa_ratio,
                    mu_oa*(1-VE_oa), mu_oa * mu_oa_ratio*(1-VE_oa))
  mu_oa_amp_g  <- c(mu_oa_amp_ref, mu_oa_amp_ref,
                    mu_oa_amp_vax, mu_oa_amp_vax)
  is_vax_g   <- c(FALSE,  FALSE,   TRUE,  TRUE)
  risk_obs_g <- c("low", "high", "low", "high")
  nm_g       <- c("Unvax low-risk", "Unvax high-risk",
                  "Vax low-risk",   "Vax high-risk")

  if (within_delta == 0) {
    # 4 strata — no within-group split
    make_params(
      name = name, R0 = R0,
      f              = f_g,
      sus            = sus_g,
      inf_mod        = rep(1, 4),
      seed           = c(5, 5, 0, 0),
      p_seek_case    = p_case_g,
      p_seek_control = p_ctrl_g,
      mu_oa          = mu_oa_g,
      is_vaccinated  = is_vax_g,
      strat_names    = nm_g,
      risk_obs       = risk_obs_g,
      mu_oa_amp      = mu_oa_amp_g,
      ...
    )
  } else {
    # 8 strata — each group split into sub-A (low sus) and sub-B (high sus)
    sus_8         <- as.vector(rbind(sus_g*(1-within_delta), sus_g*(1+within_delta)))
    f_8           <- as.vector(rbind(f_g/2, f_g/2))
    p_case_8      <- as.vector(rbind(p_case_g, p_case_g))
    p_ctrl_8      <- as.vector(rbind(p_ctrl_g, p_ctrl_g))
    mu_oa_8       <- as.vector(rbind(mu_oa_g,       mu_oa_g))
    mu_oa_amp_8   <- as.vector(rbind(mu_oa_amp_g,   mu_oa_amp_g))
    is_vax_8      <- as.vector(rbind(is_vax_g, is_vax_g))
    risk_obs_8    <- as.vector(rbind(risk_obs_g, risk_obs_g))
    nm_8          <- as.vector(rbind(paste0(nm_g, " (A)"), paste0(nm_g, " (B)")))
    seed_8        <- c(3, 2, 3, 2, 0, 0, 0, 0)

    make_params(
      name = name, R0 = R0,
      f              = f_8,
      sus            = sus_8,
      inf_mod        = rep(1, 8),
      seed           = seed_8,
      p_seek_case    = p_case_8,
      p_seek_control = p_ctrl_8,
      mu_oa          = mu_oa_8,
      is_vaccinated  = is_vax_8,
      strat_names    = nm_8,
      risk_obs       = risk_obs_8,
      mu_oa_amp      = mu_oa_amp_8,
      ...
    )
  }
}

# ---------------------------------------------------------------------------
# Calibrate mu_oa: find the rate of non-influenza-illness seeking such that
# at R0=1.5, f_v=0.2, equal risk groups (within_delta=0) the overall season-end
# test positivity [T_pos / (T_pos + T_neg)] = 35%.
#
# Formula (uniform p_seek_case, p_seek_control, mu_oa):
#   positivity = CI_total * p_seek_case
#              / (CI_total * p_seek_case + N * mu_oa * T_study * p_seek_control)
#   => mu_oa = CI_total * p_seek_case * (1 - target)
#            / (target * N * T_study * p_seek_control)
# ---------------------------------------------------------------------------
.calibrate_mu_oa <- function(
    R0_cal         = 1.5,
    f_v            = 0.2,
    f_highrisk     = 0.5,
    VE_s           = 0.6,
    sus_ratio      = 2.0,
    p_seek_case    = 0.2,
    p_seek_control = 0.4,
    N              = 100000,
    gamma          = 1 / 2,
    T_study        = 150,
    target_pos     = 0.35
) {
  p_cal <- make_flu_scenario(
    "calibration", R0 = R0_cal, f_v = f_v,
    f_highrisk_unvax = f_highrisk, f_highrisk_vax = f_highrisk,
    VE_s = VE_s, sus_ratio = sus_ratio, within_delta = 0,
    p_seek_case_unvax = p_seek_case, p_seek_case_vax = p_seek_case,
    p_seek_control = p_seek_control,
    mu_oa = 0, N = N, gamma = gamma, T_study = T_study
  )
  sim_cal  <- run_simulation(p_cal)
  n        <- p_cal$n_strat
  CI_last  <- as.numeric(tail(
    as.matrix(sim_cal$ode_out[, paste0("CI_", seq_len(n)), drop = FALSE]), 1
  ))
  CI_total  <- sum(CI_last)
  pos_tests <- CI_total * p_seek_case
  mu_oa     <- pos_tests * (1 - target_pos) /
               (target_pos * N * T_study * p_seek_control)
  message(sprintf(
    "  mu_oa_cal = %.5f  (CI_total/N = %.1f%%,  positivity target = %.0f%%)",
    mu_oa, 100 * CI_total / N, 100 * target_pos
  ))
  mu_oa
}

message("Calibrating mu_oa for 35% positivity at R0 = 1.5 ...")
mu_oa_cal <- .calibrate_mu_oa()

# ---------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------
scenarios <- list(

  # 1. R0=1.3 — lower flu season
  #    2 risk groups, no within-group heterogeneity, uniform test-seeking
  #    All estimators expected close to VE_s = 0.6
  s1 = make_flu_scenario(
    "1. R0 = 1.3", R0 = 1.3
  ),

  # 2. R0=1.5 — moderate flu season  [calibration target: 35% positivity]
  #    2 risk groups, no within-group heterogeneity, uniform test-seeking
  s2 = make_flu_scenario(
    "2. R0 = 1.5", R0 = 1.5
  ),

  # 3. Unobserved heterogeneity — R0=1.5, within_delta = 0.75 (factor-7 sub-groups)
  #    Each risk×vax cell split into sub-A (sus * 0.25) and sub-B (sus * 1.75).
  #    High-sus sub-B individuals are infected first → the pool of remaining
  #    susceptibles gradually shifts toward sub-A → instantaneous HR rises over
  #    time (apparent VE waning), even though per-stratum VE_s = 0.6 is fixed.
  #    TND identity still holds (uniform seek/mu_oa) → VE_TND = VE_RR.
  s3 = make_flu_scenario(
    "3. Unobs. heterog.", R0 = 1.5,
    within_delta = 0.75
  ),

  # 4. Indication confounding — R0=1.5, high-risk disproportionately vaccinated
  #    f_highrisk_vax=0.7 vs f_highrisk_unvax=0.3: vaccinated group has higher
  #    mean susceptibility → all population-level estimators biased toward zero.
  #    VE_HR_analytical = 1 - (1-VE_s)*(1+0.7)/(1+0.3) ≈ 0.48 (vs true 0.60)
  s4 = make_flu_scenario(
    "4. Indication bias", R0 = 1.5,
    f_highrisk_unvax = 0.3,
    f_highrisk_vax   = 0.7
  ),

  # 5. Care-seeking bias — R0=1.5, unvaccinated less likely to seek healthcare
  #    p_seek_case_unvax=0.4 vs p_seek_case_vax=0.7, with p_seek_control scaled
  #    proportionally so that p_seek_case/p_seek_control is the same for both groups.
  #    TND is designed for this scenario: the test-negative controls absorb the
  #    proportional seeking difference → VE_TND = VE_RR (unbiased).
  #    VE_HR (using N as denominator) remains biased by the seeking ratio (1.75×).
  s5 = make_flu_scenario(
    "5. Care-seeking bias", R0 = 1.5,
    p_seek_case_unvax    = 0.1,
    p_seek_case_vax      = 0.2,
    p_seek_control_unvax = 0.4 * (0.1 / 0.2)   # = p_seek_control * (p_case_unvax/p_case_vax)
  ),

  # 6. Other-illness confounding — R0=1.5, high-risk disproportionately vaccinated
  #    (same indication bias as s4) AND high-risk have 3× the other-illness rate.
  #    OR_TND = RR × (mu_oa_ref_bar / mu_oa_vax_bar)
  #           = RR × (0.7*1 + 0.3*3) / (0.3*1 + 0.7*3)  [unvax 30% high, vax 70% high]
  #           = RR × 1.6 / 2.4  = RR × 2/3
  #    → VE_TND_unadj biased in the OPPOSITE direction to VE_RR_unadj.
  #    Adjusting for observed risk group (risk_obs) removes the confounding:
  #    within each risk stratum mu_oa is the same → VE_TND_adj ≈ VE_s.
  s6 = make_flu_scenario(
    "6. Other-illness confounding", R0 = 1.5,
    f_highrisk_unvax = 0.3,
    f_highrisk_vax   = 0.7,
    mu_oa_ratio      = 3
  ),

  # 7. VE on other illness (VE_oa) — R0=1.5, vaccine halves background illness rate.
  #    Vaccinated individuals have 50% fewer non-influenza episodes → fewer TND controls
  #    in the vax arm → OR_TND inflated → VE_TND overestimates true VE_s = 0.6.
  #    VE_RR and VE_HR are unaffected (they don't use the control stream).
  #    The bias equals: VE_TND ≈ 1 - RR * (1 - VE_oa) = 1 - 0.4 * 0.5 = 0.80
  s7 = make_flu_scenario(
    "7. VE on other illness", R0 = 1.5,
    VE_oa = 0.5
  ),

  # 8. Differential seasonal amplitude — R0=1.5, reference group strongly seasonal
  #    (mu_oa_amp_ref = 0.8: background peaks at 1.8× mid-season) while vaccinated
  #    group is flat (mu_oa_amp_vax = 0.0).
  #    With equal amplitudes the seasonal factor F(T) cancels in the OR → OR = RR.
  #    Here F_ref(T) ≠ F_vax(T) so cumulative TND OR ≠ RR; the difference grows
  #    monotonically through the season.  Week-adjustment can partially correct this
  #    because within each week the within-week seasonal factor is approximately equal.
  s8 = make_flu_scenario(
    "8. Differential seasonal", R0 = 1.5,
    mu_oa_amp_ref = 0.8,
    mu_oa_amp_vax = 0.0
  )
)

# ---------------------------------------------------------------------------
# Run all scenarios
# ---------------------------------------------------------------------------
all_results <- lapply(scenarios, function(p) {
  message("\n=== Running scenario: ", p$name, " ===")
  sim     <- run_simulation(p)
  ve      <- compute_ve(sim)
  ve_stat <- compute_ve_statistical(sim)
  list(sim = sim, ve = ve, ve_stat = ve_stat, params = p)
})

# Collate to summary data frame (one row per scenario)
results_df <- do.call(rbind, lapply(all_results, function(res) {
  ve <- res$ve
  p  <- res$params
  data.frame(
    scenario         = p$name,
    R0               = p$R0,
    n_strat          = p$n_strat,
    f_vax            = sum(p$f[p$is_vaccinated]),
    AR_ref           = ve$AR_ref,
    AR_vax           = ve$AR_vax,
    VE_RR            = ve$VE_RR,
    VE_HR_analytical = ve$VE_HR_analytical,
    VE_HR_empirical  = ve$VE_HR_empirical,
    VE_TND           = ve$VE_TND,
    OR_TND           = ve$OR_TND,
    HR_sd            = ve$HR_sd,
    stringsAsFactors = FALSE
  )
}))
