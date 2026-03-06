# R/measures.R
#
# compute_ve(sim)  ── ORACLE / GROUND TRUTH (uses ODE internal state) ─────────
#   All quantities derived from ODE compartments (CI, S, h[]) that are NOT
#   observable in a real study.  Intended as simulation ground truth only.
#   Returns: VE_RR (from CI), VE_HR_analytical (from sus[]), VE_HR_empirical
#            (from h[]), VE_TND (from T_pos/T_neg — the one observed quantity),
#            AR_ref, AR_vax, diagnostics, ve_ts time-series.
#
# compute_ve_statistical(sim)  ── OBSERVED DATA ONLY ──────────────────────────
#   All quantities derived from T_pos, T_neg, and N — exactly what a real study
#   would observe.  This is the PRIMARY comparison of VE estimators.
#   VE_RR  : Poisson rate ratio  (cases = T_pos, denominator = N)
#   VE_HR  : Poisson hazard ratio (cases = ΔT_pos per bin, denominator = N)
#   VE_TND : Logistic OR          (cases = T_pos, controls = T_neg)
#   All three use the same numerator (T_pos) but different denominators.
#
# Strata are partitioned into reference (is_vaccinated=FALSE) and vaccinated
# (is_vaccinated=TRUE); all estimators compare the aggregate vaccinated group
# vs the aggregate reference group.

compute_ve <- function(sim) {
  ode_out <- sim$ode_out
  params  <- sim$params
  N_strat <- sim$N_strat
  T_study <- sim$T_study
  n       <- params$n_strat

  vax_idx <- which( params$is_vaccinated)
  ref_idx <- which(!params$is_vaccinated)

  N_vax <- sum(N_strat[vax_idx])
  N_ref <- sum(N_strat[ref_idx])

  # Extract CI and S matrices: rows = time steps, cols = strata
  CI_mat <- as.matrix(ode_out[, paste0("CI_",  seq_len(n)), drop = FALSE])
  S_mat  <- as.matrix(ode_out[, paste0("S_",   seq_len(n)), drop = FALSE])

  # -----------------------------------------------------------------------
  # VE_RR: 1 - cumulative attack rate ratio at end of study
  # -----------------------------------------------------------------------
  CI_last  <- as.numeric(tail(CI_mat, 1))
  AR_vax   <- sum(CI_last[vax_idx]) / N_vax
  AR_ref   <- sum(CI_last[ref_idx]) / N_ref
  VE_RR    <- 1 - AR_vax / AR_ref

  # -----------------------------------------------------------------------
  # VE_HR: hazard ratio from odin's per-stratum instantaneous hazard h[i] = sus[i]*lambda
  #
  # The aggregate group hazard at time t is the S-weighted mean of h[i](t)
  # over each group's strata:
  #   h_group(t) = sum_i( h[i](t) * S[i](t) ) / sum_i( S[i](t) )
  # This is the rate at which the group loses susceptibles per susceptible,
  # read directly from the ODE rather than approximated by finite differences.
  # -----------------------------------------------------------------------
  h_mat <- as.matrix(ode_out[, paste0("h_", seq_len(n)), drop = FALSE])
  S_vax <- rowSums(S_mat[, vax_idx, drop = FALSE])
  S_ref <- rowSums(S_mat[, ref_idx, drop = FALSE])
  eps   <- 1e-10

  # Numerators: sum of h[i](t)*S[i](t) per group
  hS_vax <- rowSums(h_mat[, vax_idx, drop = FALSE] * S_mat[, vax_idx, drop = FALSE])
  hS_ref <- rowSums(h_mat[, ref_idx, drop = FALSE] * S_mat[, ref_idx, drop = FALSE])

  h_vax_ts <- hS_vax / pmax(S_vax, eps)   # aggregate group hazard, vaccinated
  h_ref_ts <- hS_ref / pmax(S_ref, eps)   # aggregate group hazard, reference

  HR_ts <- h_vax_ts / pmax(h_ref_ts, eps)

  # Cox-style weighted average HR: weight each time point by the reference
  # group's event rate  h_ref(t) * S_ref(t), i.e. the -dS_ref/dt equivalent
  ref_event_rate  <- hS_ref   # = h_ref(t) * S_ref(t)
  pos             <- ref_event_rate > 0
  HR_weighted     <- if (any(pos)) sum(HR_ts[pos] * ref_event_rate[pos]) / sum(ref_event_rate[pos]) else NA_real_
  VE_HR_empirical <- 1 - HR_weighted
  HR_sd           <- if (any(pos)) sd(HR_ts[pos]) else NA_real_

  # Analytical VE_HR: N_strat-weighted susceptibility ratio between groups.
  # This equals HR(t) only when within-group susceptibility is uniform; with
  # within-group heterogeneity the empirical HR drifts over time (scenario 3)
  # because high-sus individuals deplete faster, shifting the composition of
  # remaining susceptibles.  VE_HR_analytical is therefore the "intended"
  # per-stratum VE translated to group level, not the realized average.
  sus_vax          <- weighted.mean(params$sus[vax_idx], N_strat[vax_idx])
  sus_ref          <- weighted.mean(params$sus[ref_idx], N_strat[ref_idx])
  VE_HR_analytical <- 1 - sus_vax / sus_ref

  # -----------------------------------------------------------------------
  # VE_TND: OR from 2×2 test-negative table, read directly from odin accumulators
  # T_pos[i]: cumulative positive tests (infections that sought care)
  # T_neg[i]: cumulative negative tests (other-illness episodes that sought care)
  # -----------------------------------------------------------------------
  T_pos_mat <- as.matrix(ode_out[, paste0("T_pos_", seq_len(n)), drop = FALSE])
  T_neg_mat <- as.matrix(ode_out[, paste0("T_neg_", seq_len(n)), drop = FALSE])

  T_pos_last <- as.numeric(tail(T_pos_mat, 1))
  T_neg_last <- as.numeric(tail(T_neg_mat, 1))

  cases_vax    <- sum(T_pos_last[vax_idx])
  cases_ref    <- sum(T_pos_last[ref_idx])
  controls_vax <- sum(T_neg_last[vax_idx])
  controls_ref <- sum(T_neg_last[ref_idx])

  OR_TND  <- (cases_vax / controls_vax) / (cases_ref / controls_ref)
  VE_TND  <- 1 - OR_TND

  # TND identity: when seeking behaviour is uniform OR_TND = AR_vax/AR_ref exactly
  all_seek_uniform <- (length(unique(params$p_seek_case))    == 1 &&
                       length(unique(params$p_seek_control)) == 1 &&
                       length(unique(params$mu_oa))          == 1)
  tnd_identity_diff <- if (all_seek_uniform) abs(OR_TND - AR_vax / AR_ref) else NA_real_

  # -----------------------------------------------------------------------
  # VE time series (per day)
  #
  # VE_RR(t) = 1 - AR_vax(t) / AR_ref(t)
  #
  # VE_TND(t): within each group, controls ∝ t and t cancels in the OR:
  #   OR_TND(t) = [Σ_vax CI_i(t)*p_case_i  /  ctrl_rate_vax] /
  #               [Σ_ref CI_i(t)*p_case_i  /  ctrl_rate_ref]
  # where ctrl_rate = Σ N_i * mu_oa_i * p_ctrl_i  (time-invariant).
  #
  # VE_HR(t) = 1 - HR_ts(t)  (length n_t-1, padded with NA at t=0)
  # -----------------------------------------------------------------------
  AR_vax_ts <- rowSums(CI_mat[, vax_idx, drop = FALSE]) / N_vax
  AR_ref_ts <- rowSums(CI_mat[, ref_idx, drop = FALSE]) / N_ref

  # Mask time points with fewer than 1 expected case in reference group
  active <- AR_ref_ts * N_ref >= 1

  VE_RR_ts <- ifelse(active, 1 - AR_vax_ts / AR_ref_ts, NA_real_)

  # VE_TND time series: read directly from odin's T_pos / T_neg accumulators
  T_pos_vax_ts <- rowSums(T_pos_mat[, vax_idx, drop = FALSE])
  T_pos_ref_ts <- rowSums(T_pos_mat[, ref_idx, drop = FALSE])
  T_neg_vax_ts <- rowSums(T_neg_mat[, vax_idx, drop = FALSE])
  T_neg_ref_ts <- rowSums(T_neg_mat[, ref_idx, drop = FALSE])

  OR_TND_ts <- (T_pos_vax_ts / pmax(T_neg_vax_ts, eps)) /
               (T_pos_ref_ts / pmax(T_neg_ref_ts, eps))
  VE_TND_ts <- ifelse(active, 1 - OR_TND_ts, NA_real_)

  # VE_HR: same length as other series; mask time points where reference hazard
  # is negligible (pre-epidemic and post-epidemic near-zero lambda)
  VE_HR_ts <- ifelse(h_ref_ts > eps, 1 - HR_ts, NA_real_)

  # Test positivity = T_pos / (T_pos + T_neg) per group
  # Mask t=0 where no tests have accumulated yet
  T_total_vax_ts <- T_pos_vax_ts + T_neg_vax_ts
  T_total_ref_ts <- T_pos_ref_ts + T_neg_ref_ts
  pos_vax_ts <- ifelse(T_total_vax_ts > 0, T_pos_vax_ts / T_total_vax_ts, NA_real_)
  pos_ref_ts <- ifelse(T_total_ref_ts > 0, T_pos_ref_ts / T_total_ref_ts, NA_real_)

  ve_ts <- data.frame(
    t       = ode_out$t,
    VE_RR   = VE_RR_ts,
    VE_HR   = VE_HR_ts,
    VE_TND  = VE_TND_ts,
    pos_vax = pos_vax_ts,
    pos_ref = pos_ref_ts,
    h_ref   = h_ref_ts,
    h_vax   = h_vax_ts
  )

  list(
    VE_RR            = VE_RR,
    VE_HR_analytical = VE_HR_analytical,
    VE_HR_empirical  = VE_HR_empirical,
    VE_TND           = VE_TND,
    AR_ref           = AR_ref,
    AR_vax           = AR_vax,
    OR_TND           = OR_TND,
    HR_ts            = HR_ts,
    HR_weighted      = HR_weighted,
    HR_sd            = HR_sd,
    cases_ref        = cases_ref,
    cases_vax        = cases_vax,
    controls_ref     = controls_ref,
    controls_vax     = controls_vax,
    ve_ts            = ve_ts,
    diagnostics      = list(
      HR_sd             = HR_sd,
      tnd_identity_diff = tnd_identity_diff,
      AR_ref            = AR_ref,
      AR_vax            = AR_vax
    )
  )
}

# ---------------------------------------------------------------------------
# compute_ve_statistical(sim)
#
# Estimates VE using statistical models restricted to TESTED individuals only.
# Three estimators, each unadjusted and adjusted for observed risk group:
#
#   VE_RR_tested  — Poisson regression on test-positivity rate
#                   (T_pos ~ vax [+ risk_obs] + offset(log(T_total)))
#                   Estimates positivity rate ratio among tested.
#
#   VE_TND        — Weighted logistic regression on case/control status
#                   (case ~ vax [+ risk_obs], weighted by counts)
#                   Estimates OR; equivalent to RR when positivity is low.
#
#   VE_HR         — Piecewise Poisson counting process (equivalent to Cox)
#                   Events: daily new positive tests per stratum
#                   At-risk: S_i(t) × p_seek_case_i  [susceptibles × test propensity]
#                   Baseline hazard absorbed by time-bin fixed effects (6 bins)
#                   Both an overall estimate and a per-bin time-varying series.
#
# Returns a list with point estimates, 95% Wald CIs, and VE_HR_tv data frame.
# ---------------------------------------------------------------------------
compute_ve_statistical <- function(sim) {
  ode <- sim$ode_out
  p   <- sim$params
  Ns  <- sim$N_strat
  n   <- p$n_strat

  risk_obs <- p$risk_obs        # "low"/"high" per stratum
  vax      <- p$is_vaccinated   # logical, length n

  # ── End-of-season cumulative test counts per stratum ──────────────────────
  last <- function(prefix)
    as.numeric(tail(
      as.matrix(ode[, paste0(prefix, seq_len(n)), drop = FALSE]), 1))

  T_pos <- last("T_pos_")
  T_neg <- last("T_neg_")
  T_tot <- T_pos + T_neg

  cs <- data.frame(
    vax      = vax,
    risk_obs = factor(risk_obs, c("low", "high")),
    T_pos    = T_pos,
    T_neg    = T_neg,
    T_tot    = T_tot,
    N        = Ns,
    stringsAsFactors = FALSE
  )
  cs <- cs[cs$T_tot > 0, ]   # drop strata with no tests

  # Helper: Wald VE + 95% CI from a fitted glm
  .ve <- function(fit) {
    b  <- coef(fit)["vaxTRUE"]
    se <- sqrt(vcov(fit)["vaxTRUE", "vaxTRUE"])
    list(ve = unname(1 - exp(b)),
         lo = unname(1 - exp(b + 1.96 * se)),
         hi = unname(1 - exp(b - 1.96 * se)))
  }

  # ODE outputs are continuous; Poisson/binomial GLMs warn about non-integer
  # counts — suppress these harmless warnings for the duration of this function
  oldwarn <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = oldwarn), add = TRUE)

  # ── 1. Poisson: cumulative incidence rate ratio (VE_RR) ───────────────────
  # Cases = positive tests; at-risk = population size (N per stratum).
  # Estimates cases/person: proportional to true AR when p_seek_case is uniform,
  # biased by p_seek_case ratio when it differs between groups.
  m_rr_u <- glm(T_pos ~ vax             + offset(log(N)), poisson, cs)
  m_rr_a <- glm(T_pos ~ vax + risk_obs  + offset(log(N)), poisson, cs)
  rr_u   <- .ve(m_rr_u)
  rr_a   <- .ve(m_rr_a)

  # ── 2. Logistic: OR — VE_TND ──────────────────────────────────────────────
  long <- rbind(
    data.frame(case = 1L, vax = cs$vax, risk_obs = cs$risk_obs, w = cs$T_pos),
    data.frame(case = 0L, vax = cs$vax, risk_obs = cs$risk_obs, w = cs$T_neg)
  )
  long <- long[long$w > 0, ]
  m_tnd_u <- glm(case ~ vax,            binomial, long, weights = w)
  m_tnd_a <- glm(case ~ vax + risk_obs, binomial, long, weights = w)
  tnd_u   <- .ve(m_tnd_u)
  tnd_a   <- .ve(m_tnd_a)

  # ── 3. Counting-process models: VE_HR (Poisson/N) and VE_TND_tv (logistic) ─
  # Events:    daily new positive tests per stratum            [from TND data]
  # Controls:  daily new negative tests per stratum            [from TND data]
  # At-risk:   known population size per stratum               [assumed known]
  #
  # VE_HR  — Poisson, offset = log(N): rate of positive tests per person per day.
  #           Time-bin FEs absorb baseline; estimates hazard ratio (biased by
  #           p_seek_case ratio when it differs between groups).
  # VE_TND_tv — logistic on incident cases/controls with time-bin FEs (week-
  #           adjusted TND).  Controls absorb proportional seeking; week FEs
  #           make it approximate the instantaneous hazard ratio rather than the
  #           cumulative incidence ratio.
  Tpos_mat <- as.matrix(ode[, paste0("T_pos_", seq_len(n)), drop = FALSE])
  Tneg_mat <- as.matrix(ode[, paste0("T_neg_", seq_len(n)), drop = FALSE])

  dT_pos <- apply(Tpos_mat, 2, diff)   # daily new positives, T_study × n
  dT_neg <- apply(Tneg_mat, 2, diff)   # daily new negatives, T_study × n
  t_days <- ode$t[-1]

  # Equal time bins for piecewise Poisson (baseline hazard absorbed by FEs)
  n_bins   <- 6L
  t_breaks <- seq(0, sim$T_study, length.out = n_bins + 1L)
  t_bin    <- cut(t_days, t_breaks, include.lowest = TRUE)
  t_mids   <- (t_breaks[-length(t_breaks)] + t_breaks[-1L]) / 2
  bins     <- levels(t_bin)

  # Build unified counting-process data frame: one row per stratum per day
  cp <- do.call(rbind, lapply(seq_len(n), function(i) {
    data.frame(
      t_bin    = t_bin,
      new_pos  = dT_pos[, i],
      new_neg  = dT_neg[, i],
      N        = Ns[i],
      vax      = vax[i],
      risk_obs = factor(risk_obs[i], c("low", "high"))
    )
  }))
  cp <- cp[!is.na(cp$t_bin), ]

  # Overall VE_HR (Poisson, N offset, time-bin FEs)
  m_hr_u <- glm(new_pos ~ vax            + t_bin + offset(log(N)), poisson, cp)
  m_hr_a <- glm(new_pos ~ vax + risk_obs + t_bin + offset(log(N)), poisson, cp)
  hr_u   <- .ve(m_hr_u)
  hr_a   <- .ve(m_hr_a)

  # Overall VE_TND time-adjusted (week-adjusted logistic on incident counts)
  long_cp <- rbind(
    data.frame(case = 1L, cp[, c("t_bin", "vax", "risk_obs")], w = cp$new_pos),
    data.frame(case = 0L, cp[, c("t_bin", "vax", "risk_obs")], w = cp$new_neg)
  )
  long_cp <- long_cp[long_cp$w > 0, ]
  m_tnd_tv_u <- glm(case ~ vax            + t_bin, binomial, long_cp, weights = w)
  m_tnd_tv_a <- glm(case ~ vax + risk_obs + t_bin, binomial, long_cp, weights = w)
  tnd_tv_u   <- .ve(m_tnd_tv_u)
  tnd_tv_a   <- .ve(m_tnd_tv_a)

  # Per-bin time-varying VE_HR and VE_TND_tv (both adjusted for risk_obs)
  tv <- do.call(rbind, lapply(seq_along(bins), function(b) {
    d <- cp[cp$t_bin == bins[b], ]
    na_row <- data.frame(t_mid  = t_mids[b],
                         VE_HR  = NA_real_, HR_lo  = NA_real_, HR_hi  = NA_real_,
                         VE_TND = NA_real_, TND_lo = NA_real_, TND_hi = NA_real_)
    if (sum(d$new_pos) < 0.5 || length(unique(d$vax)) < 2) return(na_row)
    tryCatch({
      # VE_HR per bin
      mh <- glm(new_pos ~ vax + risk_obs + offset(log(N)), poisson, d)
      vh <- .ve(mh)
      # VE_TND per bin (logistic, incident counts, no extra time FE within bin)
      dl <- rbind(
        data.frame(case = 1L, vax = d$vax, risk_obs = d$risk_obs, w = d$new_pos),
        data.frame(case = 0L, vax = d$vax, risk_obs = d$risk_obs, w = d$new_neg)
      )
      dl <- dl[dl$w > 0, ]
      mt <- glm(case ~ vax + risk_obs, binomial, dl, weights = w)
      vt <- .ve(mt)
      data.frame(t_mid  = t_mids[b],
                 VE_HR  = vh$ve, HR_lo  = vh$lo, HR_hi  = vh$hi,
                 VE_TND = vt$ve, TND_lo = vt$lo, TND_hi = vt$hi)
    }, error = function(e) na_row)
  }))

  # ── 4. Weekly estimates (adjusted for risk_obs; point estimates only) ───────
  # For each calendar week of the season:
  #   VE_RR_week  — Poisson on cumulative T_pos ~ vax + risk_obs + offset(log(N))
  #   VE_TND_cum  — logistic OR on cumulative T_pos / T_neg up to end of week
  #   VE_TND_inc  — logistic OR on incident T_pos / T_neg for that week only
  #   VE_HR_week  — Poisson, new T_pos ~ vax + risk_obs + offset(log(N − cum.T_pos))

  t_week_breaks <- c(seq(0, sim$T_study, by = 7))
  if (tail(t_week_breaks, 1) < sim$T_study)
    t_week_breaks <- c(t_week_breaks, sim$T_study)
  n_weeks     <- length(t_week_breaks) - 1
  t_week_mids <- (t_week_breaks[-length(t_week_breaks)] + t_week_breaks[-1]) / 2

  # Point-estimate only extractor (no vcov needed for deterministic model)
  .ve_pt <- function(fit) unname(1 - exp(coef(fit)["vaxTRUE"]))

  weekly <- do.call(rbind, lapply(seq_len(n_weeks), function(w) {
    t_start <- t_week_breaks[w]
    t_end   <- t_week_breaks[w + 1]
    t_mid   <- t_week_mids[w]
    row_s   <- which(ode$t == t_start)
    row_e   <- which(ode$t == t_end)

    na_row <- data.frame(t_mid = t_mid,
                         VE_RR_week = NA_real_,
                         VE_TND_cum = NA_real_,
                         VE_TND_inc = NA_real_,
                         VE_HR_week = NA_real_)

    if (!length(row_s) || !length(row_e)) return(na_row)

    cum_pos      <- as.numeric(Tpos_mat[row_e, ])
    cum_neg      <- as.numeric(Tneg_mat[row_e, ])
    inc_pos      <- cum_pos - as.numeric(Tpos_mat[row_s, ])
    inc_neg      <- cum_neg - as.numeric(Tneg_mat[row_s, ])
    # Observable at-risk: N minus cumulative positive tests at week start.
    # Uses only TND data (no ODE susceptibles); underestimates true removed
    # (since T_pos = CI * p_seek_case < CI), so at_risk_obs > true S.
    at_risk_obs  <- pmax(Ns - as.numeric(Tpos_mat[row_s, ]), 1)

    dw <- data.frame(
      vax      = vax,
      risk_obs = factor(risk_obs, c("low", "high")),
      cum_pos  = cum_pos,  cum_neg = cum_neg,
      inc_pos  = inc_pos,  inc_neg = inc_neg,
      at_risk  = at_risk_obs,
      N        = Ns
    )
    if (length(unique(dw$vax)) < 2) return(na_row)

    tryCatch({
      # Cumulative RR (Poisson, N offset, adjusted for risk group)
      vr  <- if (sum(dw$cum_pos[dw$vax]) < 0.5 || sum(dw$cum_pos[!dw$vax]) < 0.5) NA_real_ else
               .ve_pt(glm(cum_pos ~ vax + risk_obs + offset(log(N)), poisson, dw))

      # Cumulative OR (adjusted)
      lc  <- rbind(data.frame(case = 1L, vax = dw$vax, risk_obs = dw$risk_obs, w = dw$cum_pos),
                   data.frame(case = 0L, vax = dw$vax, risk_obs = dw$risk_obs, w = dw$cum_neg))
      lc  <- lc[lc$w > 0, ]
      vc  <- if (length(unique(lc$vax)) < 2) NA_real_ else
               .ve_pt(glm(case ~ vax + risk_obs, binomial, lc, weights = w))

      # Incident OR (adjusted)
      li  <- rbind(data.frame(case = 1L, vax = dw$vax, risk_obs = dw$risk_obs, w = dw$inc_pos),
                   data.frame(case = 0L, vax = dw$vax, risk_obs = dw$risk_obs, w = dw$inc_neg))
      li  <- li[li$w > 0, ]
      vi  <- if (length(unique(li$vax)) < 2 ||
                 sum(dw$inc_pos[dw$vax]) < 0.5 || sum(dw$inc_pos[!dw$vax]) < 0.5) NA_real_ else
               .ve_pt(glm(case ~ vax + risk_obs, binomial, li, weights = w))

      # Weekly HRR (Poisson, at-risk = N − cumulative positive tests, adjusted)
      vh  <- if (sum(dw$inc_pos[dw$vax]) < 0.5 || sum(dw$inc_pos[!dw$vax]) < 0.5) NA_real_ else
               .ve_pt(glm(inc_pos ~ vax + risk_obs + offset(log(at_risk)), poisson, dw))

      data.frame(t_mid = t_mid, VE_RR_week = vr,
                 VE_TND_cum = vc, VE_TND_inc = vi, VE_HR_week = vh)
    }, error = function(e) na_row)
  }))

  list(
    VE_RR_unadj      = rr_u$ve,      VE_RR_unadj_lo    = rr_u$lo,      VE_RR_unadj_hi    = rr_u$hi,
    VE_RR_adj        = rr_a$ve,      VE_RR_adj_lo      = rr_a$lo,      VE_RR_adj_hi      = rr_a$hi,
    VE_TND_unadj     = tnd_u$ve,     VE_TND_unadj_lo   = tnd_u$lo,     VE_TND_unadj_hi   = tnd_u$hi,
    VE_TND_adj       = tnd_a$ve,     VE_TND_adj_lo     = tnd_a$lo,     VE_TND_adj_hi     = tnd_a$hi,
    VE_TND_tv_unadj  = tnd_tv_u$ve,  VE_TND_tv_unadj_lo = tnd_tv_u$lo, VE_TND_tv_unadj_hi = tnd_tv_u$hi,
    VE_TND_tv_adj    = tnd_tv_a$ve,  VE_TND_tv_adj_lo  = tnd_tv_a$lo,  VE_TND_tv_adj_hi  = tnd_tv_a$hi,
    VE_HR_unadj      = hr_u$ve,      VE_HR_unadj_lo    = hr_u$lo,      VE_HR_unadj_hi    = hr_u$hi,
    VE_HR_adj        = hr_a$ve,      VE_HR_adj_lo      = hr_a$lo,      VE_HR_adj_hi      = hr_a$hi,
    tv               = tv,           # per-bin VE_HR and VE_TND_tv (adjusted)
    weekly           = weekly,       # weekly: VE_TND_cum, VE_TND_inc, VE_HR_week (adjusted)
    df_cs            = cs
  )
}
