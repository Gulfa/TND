# test-measures.R
# Tests for compute_ve() and compute_ve_statistical() in R/measures.R.
# Uses pre-computed all_results from helper-setup.R.

# ── Output structure ─────────────────────────────────────────────────────────

test_that("compute_ve returns expected fields", {
  ve <- all_results$s1$ve
  expected <- c("VE_RR", "VE_HR_analytical", "VE_HR_empirical", "VE_TND",
                "AR_ref", "AR_vax", "OR_TND", "HR_ts", "HR_weighted", "HR_sd",
                "cases_ref", "cases_vax", "controls_ref", "controls_vax",
                "ve_ts", "diagnostics")
  expect_true(all(expected %in% names(ve)))
})

test_that("compute_ve_statistical returns expected fields", {
  vs <- all_results$s1$ve_stat
  expected <- c("VE_RR_unadj", "VE_RR_adj",
                "VE_TND_unadj", "VE_TND_adj",
                "VE_HR_unadj",  "VE_HR_adj",
                "tv", "df_cs")
  expect_true(all(expected %in% names(vs)))
})

# ── Attack rate bounds ────────────────────────────────────────────────────────

test_that("AR_ref and AR_vax are in [0, 1] for all scenarios", {
  for (nm in names(all_results)) {
    ve <- all_results[[nm]]$ve
    expect_gte(ve$AR_ref, 0, label = paste("AR_ref >=0 in", nm))
    expect_lte(ve$AR_ref, 1, label = paste("AR_ref <=1 in", nm))
    expect_gte(ve$AR_vax, 0, label = paste("AR_vax >=0 in", nm))
    expect_lte(ve$AR_vax, 1, label = paste("AR_vax <=1 in", nm))
  }
})

test_that("AR_ref > AR_vax for all scenarios (vaccine reduces infection)", {
  for (nm in names(all_results)) {
    ve <- all_results[[nm]]$ve
    expect_gt(ve$AR_ref, ve$AR_vax, label = paste("AR_ref > AR_vax in", nm))
  }
})

# ── VE_RR bounds ─────────────────────────────────────────────────────────────

test_that("VE_RR <= 1 for all scenarios", {
  for (nm in names(all_results)) {
    expect_lte(all_results[[nm]]$ve$VE_RR, 1,
               label = paste("VE_RR <= 1 in", nm))
  }
})

# ── Test counts are positive ─────────────────────────────────────────────────

test_that("cases and controls are positive for all scenarios", {
  for (nm in names(all_results)) {
    ve <- all_results[[nm]]$ve
    expect_gt(ve$cases_ref,    0, label = paste("cases_ref > 0 in",    nm))
    expect_gt(ve$cases_vax,    0, label = paste("cases_vax > 0 in",    nm))
    expect_gt(ve$controls_ref, 0, label = paste("controls_ref > 0 in", nm))
    expect_gt(ve$controls_vax, 0, label = paste("controls_vax > 0 in", nm))
  }
})

# ── TND identity ─────────────────────────────────────────────────────────────
# When all strata share the same p_seek_case, p_seek_control, and mu_oa,
# OR_TND = AR_vax / AR_ref exactly (up to ODE integration error).
# Scenarios where this holds: s1, s2, s3, s4 (uniform seeking).
# Scenarios where it does NOT hold (tnd_identity_diff = NA): s5, s6.

test_that("TND identity holds for scenarios with uniform seeking (s1, s2, s3, s4)", {
  for (nm in c("s1", "s2", "s3", "s4")) {
    diff <- all_results[[nm]]$ve$diagnostics$tnd_identity_diff
    expect_false(is.na(diff),
                 label = paste("tnd_identity_diff not NA for", nm))
    expect_lt(diff, 1e-3,
              label = paste("TND identity |OR_TND - AR_vax/AR_ref| < 1e-3 in", nm))
  }
})

test_that("tnd_identity_diff is NA for non-uniform seeking (s5, s6)", {
  for (nm in c("s5", "s6")) {
    diff <- all_results[[nm]]$ve$diagnostics$tnd_identity_diff
    expect_true(is.na(diff),
                label = paste("tnd_identity_diff is NA for", nm))
  }
})

# ── VE_HR_analytical = 1 - sus_vax / sus_ref ─────────────────────────────────
# This is a closed-form property: the analytical HR at t=0 equals
# the ratio of population-weighted mean susceptibilities.

test_that("VE_HR_analytical matches 1 - sus_vax/sus_ref for all scenarios", {
  for (nm in names(all_results)) {
    res     <- all_results[[nm]]
    p       <- res$params
    Ns      <- res$sim$N_strat
    vax_idx <- which( p$is_vaccinated)
    ref_idx <- which(!p$is_vaccinated)
    sus_vax <- weighted.mean(p$sus[vax_idx], Ns[vax_idx])
    sus_ref <- weighted.mean(p$sus[ref_idx], Ns[ref_idx])
    expected <- 1 - sus_vax / sus_ref
    expect_equal(res$ve$VE_HR_analytical, expected, tolerance = 1e-10,
                 label = paste("VE_HR_analytical in", nm))
  }
})

# ── OR_TND is positive ────────────────────────────────────────────────────────

test_that("OR_TND is positive for all scenarios", {
  for (nm in names(all_results)) {
    expect_gt(all_results[[nm]]$ve$OR_TND, 0,
              label = paste("OR_TND > 0 in", nm))
  }
})

# ── ve_ts data frame ─────────────────────────────────────────────────────────

test_that("ve_ts has T_study + 1 rows and expected columns", {
  for (nm in names(all_results)) {
    ve_ts <- all_results[[nm]]$ve$ve_ts
    T_s   <- all_results[[nm]]$sim$T_study
    expect_equal(nrow(ve_ts), T_s + 1L,
                 label = paste("ve_ts nrow in", nm))
    expect_true(all(c("t", "VE_RR", "VE_HR", "VE_TND") %in% names(ve_ts)),
                label = paste("ve_ts columns in", nm))
  }
})

# ── HR_sd: zero when within-group susceptibility is uniform ──────────────────
# Scenarios without within-group heterogeneity (s1, s2, s4, s5, s6):
# HR_sd should be close to zero because each group's hazard ratio is constant.

test_that("HR_sd is near-zero for scenarios without within-group heterogeneity", {
  for (nm in c("s1", "s2", "s4", "s5", "s6")) {
    hr_sd <- all_results[[nm]]$ve$HR_sd
    expect_lt(hr_sd, 0.05,
              label = paste("HR_sd near-zero in", nm))
  }
})

test_that("HR_sd is notably positive for s3 (within-group heterogeneity)", {
  expect_gt(all_results$s3$ve$HR_sd, 0.02)
})
