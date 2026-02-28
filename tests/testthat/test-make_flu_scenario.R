# test-make_flu_scenario.R
# Tests for make_flu_scenario(): the main parameter builder for flu VE scenarios.
# All tests are pure-function: no ODE required.

# ── Structural invariants ────────────────────────────────────────────────────

test_that("population fractions sum to 1 (4-stratum)", {
  p <- make_flu_scenario("t", R0 = 1.5)
  expect_equal(sum(p$f), 1.0, tolerance = 1e-10)
})

test_that("population fractions sum to 1 (8-stratum)", {
  p <- make_flu_scenario("t", R0 = 1.5, within_delta = 0.5)
  expect_equal(sum(p$f), 1.0, tolerance = 1e-10)
})

test_that("4 strata when within_delta = 0", {
  p <- make_flu_scenario("t", R0 = 1.5, within_delta = 0)
  expect_equal(p$n_strat, 4L)
})

test_that("8 strata when within_delta > 0", {
  p <- make_flu_scenario("t", R0 = 1.5, within_delta = 1/3)
  expect_equal(p$n_strat, 8L)
})

test_that("all susceptibilities are positive", {
  p <- make_flu_scenario("t", R0 = 1.5)
  expect_true(all(p$sus > 0))
})

test_that("all fractions are positive", {
  p <- make_flu_scenario("t", R0 = 1.5)
  expect_true(all(p$f > 0))
})

test_that("is_vaccinated has both TRUE and FALSE", {
  p <- make_flu_scenario("t", R0 = 1.5)
  expect_true(any( p$is_vaccinated))
  expect_true(any(!p$is_vaccinated))
})

# ── Susceptibility normalisation ─────────────────────────────────────────────
# Weighted mean susceptibility of the reference group must equal 1.0,
# regardless of sus_ratio or f_highrisk_unvax.

test_that("weighted-mean sus of reference group = 1 (default params)", {
  p <- make_flu_scenario("t", R0 = 1.5)
  ref  <- !p$is_vaccinated
  wmean <- weighted.mean(p$sus[ref], (p$f * p$N)[ref])
  expect_equal(wmean, 1.0, tolerance = 1e-10)
})

test_that("weighted-mean sus of reference group = 1 (sus_ratio = 3)", {
  p <- make_flu_scenario("t", R0 = 1.5, sus_ratio = 3.0)
  ref  <- !p$is_vaccinated
  wmean <- weighted.mean(p$sus[ref], (p$f * p$N)[ref])
  expect_equal(wmean, 1.0, tolerance = 1e-10)
})

test_that("weighted-mean sus of reference group = 1 (f_highrisk_unvax = 0.3)", {
  p <- make_flu_scenario("t", R0 = 1.5, f_highrisk_unvax = 0.3)
  ref  <- !p$is_vaccinated
  wmean <- weighted.mean(p$sus[ref], (p$f * p$N)[ref])
  expect_equal(wmean, 1.0, tolerance = 1e-10)
})

test_that("weighted-mean sus of reference group = 1 (8-stratum)", {
  p <- make_flu_scenario("t", R0 = 1.5, within_delta = 0.5)
  ref  <- !p$is_vaccinated
  wmean <- weighted.mean(p$sus[ref], (p$f * p$N)[ref])
  expect_equal(wmean, 1.0, tolerance = 1e-10)
})

# ── VE_s encodes correctly ───────────────────────────────────────────────────
# sus_vax = sus_ref * (1 - VE_s) for matching risk strata.
# => weighted-mean sus_vax / weighted-mean sus_ref = (1 - VE_s).

test_that("ratio of mean sus (vax/ref) equals (1 - VE_s)", {
  VE_s <- 0.6
  p    <- make_flu_scenario("t", R0 = 1.5, VE_s = VE_s,
                             f_highrisk_unvax = 0.5, f_highrisk_vax = 0.5)
  # When risk-group proportions are equal, mean_sus_vax / mean_sus_ref = (1 - VE_s)
  vax <- p$is_vaccinated; ref <- !vax
  w_vax <- (p$f * p$N)[vax]; w_ref <- (p$f * p$N)[ref]
  ratio <- weighted.mean(p$sus[vax], w_vax) / weighted.mean(p$sus[ref], w_ref)
  expect_equal(ratio, 1 - VE_s, tolerance = 1e-10)
})

test_that("VE_s = 0 gives same sus in vax and ref (same risk groups, equal proportions)", {
  p <- make_flu_scenario("t", R0 = 1.5, VE_s = 0,
                          f_highrisk_unvax = 0.5, f_highrisk_vax = 0.5)
  vax <- p$is_vaccinated; ref <- !vax
  w_vax <- (p$f * p$N)[vax]; w_ref <- (p$f * p$N)[ref]
  expect_equal(weighted.mean(p$sus[vax], w_vax),
               weighted.mean(p$sus[ref], w_ref), tolerance = 1e-10)
})

# ── Within-cell split symmetry ───────────────────────────────────────────────
# The 8-stratum split must preserve group mean susceptibility exactly.

test_that("within-cell split preserves group-mean sus (within_delta = 1/3)", {
  p4 <- make_flu_scenario("t", R0 = 1.5, within_delta = 0.0)
  p8 <- make_flu_scenario("t", R0 = 1.5, within_delta = 1/3)
  # Each group of 2 sub-strata (A, B) should average back to the 4-stratum value.
  # Groups: [1,2] unvax-low, [3,4] unvax-high, [5,6] vax-low, [7,8] vax-high
  for (g in 1:4) {
    idx8 <- c(2*g - 1, 2*g)
    mean8 <- weighted.mean(p8$sus[idx8], p8$f[idx8])
    expect_equal(mean8, p4$sus[g], tolerance = 1e-10)
  }
})

# ── Vaccination fraction ─────────────────────────────────────────────────────

test_that("vaccinated fraction matches f_v", {
  f_v <- 0.3
  p <- make_flu_scenario("t", R0 = 1.5, f_v = f_v)
  expect_equal(sum(p$f[p$is_vaccinated]), f_v, tolerance = 1e-10)
})
