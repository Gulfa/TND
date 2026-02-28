# test-simulate.R
# Tests for run_simulation(): ODE mass balance, monotonicity, and output structure.
# Uses pre-computed all_results from helper-setup.R (no extra ODE runs).

# ── Output structure ─────────────────────────────────────────────────────────

test_that("run_simulation returns list with ode_out, params, N_strat, T_study", {
  sim <- all_results$s1$sim
  expect_named(sim, c("ode_out", "params", "N_strat", "T_study"), ignore.order = TRUE)
})

test_that("ode_out has columns t, S_i, I_i, Rec_i, CI_i, h_i, T_pos_i, T_neg_i for all strata", {
  sim <- all_results$s2$sim
  n   <- sim$params$n_strat
  expected_prefixes <- c("t", paste0("S_",     seq_len(n)),
                              paste0("I_",     seq_len(n)),
                              paste0("Rec_",   seq_len(n)),
                              paste0("CI_",    seq_len(n)),
                              paste0("h_",     seq_len(n)),
                              paste0("T_pos_", seq_len(n)),
                              paste0("T_neg_", seq_len(n)))
  expect_true(all(expected_prefixes %in% names(sim$ode_out)))
})

test_that("time column runs from 0 to T_study in steps of 1", {
  sim <- all_results$s1$sim
  expect_equal(sim$ode_out$t[1], 0)
  expect_equal(tail(sim$ode_out$t, 1), sim$T_study)
  expect_equal(diff(sim$ode_out$t), rep(1, sim$T_study))
})

test_that("N_strat sums to total N", {
  sim <- all_results$s2$sim
  expect_equal(sum(sim$N_strat), sim$params$N)
})

# ── Mass balance: S + I + Rec = N_strat for every stratum at every time ──────

test_that("mass balance holds in s1 (4-stratum, no heterogeneity)", {
  sim <- all_results$s1$sim
  n   <- sim$params$n_strat
  for (i in seq_len(n)) {
    mass <- sim$ode_out[[paste0("S_",   i)]] +
            sim$ode_out[[paste0("I_",   i)]] +
            sim$ode_out[[paste0("Rec_", i)]]
    rel_err <- abs(mass - sim$N_strat[i]) / sim$N_strat[i]
    expect_lt(max(rel_err), 1e-3,
              label = sprintf("mass balance stratum %d in s1", i))
  }
})

test_that("mass balance holds in s3 (8-stratum, within-group heterogeneity)", {
  sim <- all_results$s3$sim
  n   <- sim$params$n_strat
  for (i in seq_len(n)) {
    mass <- sim$ode_out[[paste0("S_",   i)]] +
            sim$ode_out[[paste0("I_",   i)]] +
            sim$ode_out[[paste0("Rec_", i)]]
    rel_err <- abs(mass - sim$N_strat[i]) / sim$N_strat[i]
    expect_lt(max(rel_err), 1e-3,
              label = sprintf("mass balance stratum %d in s3", i))
  }
})

# ── Monotonicity ─────────────────────────────────────────────────────────────

test_that("S is non-increasing in every stratum (s2)", {
  sim <- all_results$s2$sim
  n   <- sim$params$n_strat
  for (i in seq_len(n)) {
    dS <- diff(sim$ode_out[[paste0("S_", i)]])
    expect_lte(max(dS), 1e-6,
               label = sprintf("S non-increasing, stratum %d, s2", i))
  }
})

test_that("CI (cumulative incidence) is non-decreasing in every stratum (s2)", {
  sim <- all_results$s2$sim
  n   <- sim$params$n_strat
  for (i in seq_len(n)) {
    dCI <- diff(sim$ode_out[[paste0("CI_", i)]])
    expect_gte(min(dCI), -1e-6,
               label = sprintf("CI non-decreasing, stratum %d, s2", i))
  }
})

test_that("T_pos is non-decreasing in every stratum (s2)", {
  sim <- all_results$s2$sim
  n   <- sim$params$n_strat
  for (i in seq_len(n)) {
    d <- diff(sim$ode_out[[paste0("T_pos_", i)]])
    expect_gte(min(d), -1e-6,
               label = sprintf("T_pos non-decreasing, stratum %d, s2", i))
  }
})

test_that("T_neg is non-decreasing in every stratum (s2)", {
  sim <- all_results$s2$sim
  n   <- sim$params$n_strat
  for (i in seq_len(n)) {
    d <- diff(sim$ode_out[[paste0("T_neg_", i)]])
    expect_gte(min(d), -1e-6,
               label = sprintf("T_neg non-decreasing, stratum %d, s2", i))
  }
})

# ── Seeding: vaccinated strata start with 0 infected ─────────────────────────

test_that("vaccinated strata have I = 0 at t = 0 in all scenarios", {
  for (nm in names(all_results)) {
    sim <- all_results[[nm]]$sim
    vax_idx <- which(sim$params$is_vaccinated)
    for (i in vax_idx) {
      expect_equal(sim$ode_out[[paste0("I_", i)]][1], 0,
                   label = sprintf("I_%d at t=0 in %s", i, nm))
    }
  }
})

# ── Stochastic consistency: higher R0 → higher total attack rate ─────────────

test_that("AR_total(s2, R0=1.5) > AR_total(s1, R0=1.3)", {
  ar <- function(nm) {
    sim <- all_results[[nm]]$sim
    n   <- sim$params$n_strat
    CI_last <- vapply(seq_len(n),
                      function(i) tail(sim$ode_out[[paste0("CI_", i)]], 1),
                      numeric(1))
    sum(CI_last) / sim$params$N
  }
  expect_gt(ar("s2"), ar("s1"))
})
