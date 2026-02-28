# test-scenarios.R
# End-to-end tests asserting epidemiologically expected behaviour for each scenario.
#
# Important: In a finite-time SIR, VE_RR (cumulative AR ratio) and VE_TND
# (OR from test-negative counts) are NOT equal to VE_s = 0.6 even without
# confounding.  Only VE_HR_analytical = 1 - sus_vax/sus_ref is guaranteed to
# equal VE_s by construction.  Tests therefore check relationships between
# estimators and directional biases, not absolute equality to VE_s.

# ── s1: R0 = 1.3, no confounding ─────────────────────────────────────────────
# Baseline: uniform susceptibility within groups, no indication or care-seeking bias.
# VE_HR_analytical = 0.6 (tested in test-measures.R).
# VE_RR ≈ VE_TND (TND identity, tested in test-measures.R).
# VE_RR < VE_HR_analytical due to SIR depletion dynamics.

test_that("s1: VE_RR is in plausible range and below VE_HR_analytical", {
  ve <- all_results$s1$ve
  expect_gt(ve$VE_RR, 0.45)
  expect_lt(ve$VE_RR, 0.65)
  expect_lt(ve$VE_RR, ve$VE_HR_analytical)
})

test_that("s1: statistical estimators agree (no confounding)", {
  vs <- all_results$s1$ve_stat
  # Adjusting for risk group should change nothing when groups are unconfounded
  expect_equal(vs$VE_RR_unadj,  vs$VE_RR_adj,  tolerance = 0.01)
  expect_equal(vs$VE_TND_unadj, vs$VE_TND_adj, tolerance = 0.01)
  expect_equal(vs$VE_HR_unadj,  vs$VE_HR_adj,  tolerance = 0.01)
})

# ── s2: R0 = 1.5, calibration baseline ───────────────────────────────────────
# Same structural properties as s1, but higher R0 → more depletion → lower VE_RR.

test_that("s2: VE_RR < s1 VE_RR (higher R0 means more selective depletion)", {
  expect_lt(all_results$s2$ve$VE_RR, all_results$s1$ve$VE_RR)
})

test_that("s2: VE_HR_analytical = 0.6 (instantaneous hazard ratio equals VE_s)", {
  expect_equal(all_results$s2$ve$VE_HR_analytical, 0.6, tolerance = 1e-6)
})

test_that("s2: statistical unadjusted == adjusted (no confounding to remove)", {
  vs <- all_results$s2$ve_stat
  expect_equal(vs$VE_RR_unadj,  vs$VE_RR_adj,  tolerance = 0.01)
  expect_equal(vs$VE_TND_unadj, vs$VE_TND_adj, tolerance = 0.01)
})

# ── s3: unobserved within-group heterogeneity ─────────────────────────────────
# Each risk×vax cell split into low-sus (A) and high-sus (B) sub-strata.
# High-sus sub-groups are depleted first → the S-weighted hazard ratio drifts
# over time (apparent waning) even though per-stratum VE_s = 0.6 is fixed.

test_that("s3: TND identity holds (uniform seeking) → VE_TND = VE_RR", {
  ve <- all_results$s3$ve
  expect_equal(ve$VE_TND, ve$VE_RR, tolerance = 1e-3)
})

test_that("s3: HR_sd > 0 (hazard ratio varies over time due to heterogeneity)", {
  expect_gt(all_results$s3$ve$HR_sd, 0.02)
})

test_that("s3: HR_sd exceeds s2 HR_sd (more time-variation with heterogeneity)", {
  expect_gt(all_results$s3$ve$HR_sd, all_results$s2$ve$HR_sd)
})

test_that("s3: VE_HR_empirical < VE_HR_analytical (selective depletion biases weighted HR down)", {
  ve <- all_results$s3$ve
  expect_lt(ve$VE_HR_empirical, ve$VE_HR_analytical)
})

# ── s4: indication bias (high-risk disproportionately vaccinated) ─────────────
# f_highrisk_vax = 0.7 vs f_highrisk_unvax = 0.3.
# The vaccinated group has higher mean susceptibility → all ODE estimators
# biased toward the null relative to s2 (no confounding baseline).
# Statistical adjustment for observed risk group (risk_obs) partly corrects this.

test_that("s4: VE_HR_analytical ≈ 0.48 (derived from shifted group-sus ratio)", {
  # 1 - sus_vax_mean/sus_ref_mean with f_highrisk_vax = 0.7, f_highrisk_unvax = 0.3
  expect_equal(all_results$s4$ve$VE_HR_analytical, 0.48, tolerance = 0.03)
})

test_that("s4: VE_RR < s2 VE_RR (indication bias reduces apparent efficacy)", {
  expect_lt(all_results$s4$ve$VE_RR, all_results$s2$ve$VE_RR)
})

test_that("s4: statistical VE_RR_adj > VE_RR_unadj (risk adjustment corrects some confounding)", {
  vs <- all_results$s4$ve_stat
  expect_gt(vs$VE_RR_adj, vs$VE_RR_unadj)
})

test_that("s4: statistical VE_TND_adj > VE_TND_unadj (same directional correction)", {
  vs <- all_results$s4$ve_stat
  expect_gt(vs$VE_TND_adj, vs$VE_TND_unadj)
})

test_that("s4: adjusted VE_TND closer to s2 VE_TND than unadjusted is", {
  s2_ve  <- all_results$s2$ve$VE_TND
  vs     <- all_results$s4$ve_stat
  expect_lt(abs(vs$VE_TND_adj  - s2_ve),
            abs(vs$VE_TND_unadj - s2_ve))
})

# ── s5: care-seeking bias ─────────────────────────────────────────────────────
# p_seek_case_unvax = 0.4, p_seek_case_vax = 0.7.
# p_seek_case / p_seek_control is the SAME in both groups (proportional).
# ODE VE_TND absorbs the seeking ratio → VE_TND ≈ s2 VE_TND (unbiased).
# Statistical VE_HR (N offset) is biased down; statistical VE_TND is unbiased.

test_that("s5: ODE VE_TND ≈ s2 VE_TND (TND absorbs proportional care-seeking)", {
  expect_equal(all_results$s5$ve$VE_TND, all_results$s2$ve$VE_TND, tolerance = 0.02)
})

test_that("s5: statistical VE_TND_unadj ≈ s2 stat VE_TND_unadj (unbiased)", {
  expect_equal(all_results$s5$ve_stat$VE_TND_unadj,
               all_results$s2$ve_stat$VE_TND_unadj, tolerance = 0.02)
})

test_that("s5: statistical VE_HR_unadj notably below s2 (care-seeking biases N-denominator HR)", {
  # Unvaccinated seek 0.4/0.7 ≈ 57% as often → positive test rate per person
  # understates their true incidence → apparent VE_HR is deflated.
  expect_lt(all_results$s5$ve_stat$VE_HR_unadj,
            all_results$s2$ve_stat$VE_HR_unadj - 0.10)
})

# ── s6: other-illness confounding + indication bias ───────────────────────────
# High-risk have 3× other-illness rate (mu_oa_ratio = 3) AND are
# disproportionately vaccinated (same as s4: f_highrisk_vax = 0.7).
# Effect: vaccinated group generates more negative tests from high-risk →
# OR_TND inflated relative to AR_vax/AR_ref → VE_TND pushed ABOVE VE_RR.
# Adjusting for risk group removes the mu_oa_ratio confounding, leaving only
# the indication bias (same as s4).

test_that("s6: ODE VE_TND > ODE VE_RR (other-illness confounding inflates OR_TND)", {
  ve <- all_results$s6$ve
  expect_gt(ve$VE_TND, ve$VE_RR)
})

test_that("s6: stat VE_TND_unadj > stat VE_RR_unadj (same direction in statistical analysis)", {
  vs <- all_results$s6$ve_stat
  expect_gt(vs$VE_TND_unadj, vs$VE_RR_unadj)
})

test_that("s6: stat VE_TND_adj < stat VE_TND_unadj (adjustment removes other-illness inflation)", {
  vs <- all_results$s6$ve_stat
  expect_lt(vs$VE_TND_adj, vs$VE_TND_unadj)
})

test_that("s6: stat VE_TND_adj ≈ s4 stat VE_TND_adj (same remaining indication bias)", {
  expect_equal(all_results$s6$ve_stat$VE_TND_adj,
               all_results$s4$ve_stat$VE_TND_adj, tolerance = 0.02)
})
