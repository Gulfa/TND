# test-make_params.R
# Unit tests for make_params() in scenarios.R.
# All tests are pure-function: no ODE required.

# A minimal valid argument list used as the baseline for mutation tests.
valid <- list(
  name           = "test",
  R0             = 1.5,
  f              = c(0.4, 0.4, 0.1, 0.1),
  sus            = c(1.0, 1.0, 0.4, 0.4),
  inf_mod        = c(1.0, 1.0, 1.0, 1.0),
  seed           = c(5L, 5L, 0L, 0L),
  p_seek_case    = c(0.7, 0.7, 0.7, 0.7),
  p_seek_control = c(0.4, 0.4, 0.4, 0.4),
  mu_oa          = c(1e-3, 1e-3, 1e-3, 1e-3),
  is_vaccinated  = c(FALSE, FALSE, TRUE, TRUE)
)

test_that("make_params returns list with all required fields", {
  p <- do.call(make_params, valid)
  expected_fields <- c("name", "R0", "n_strat", "N", "gamma", "beta",
                        "f", "sus", "inf_mod", "seed",
                        "p_seek_case", "p_seek_control", "mu_oa",
                        "is_vaccinated", "strat_names", "risk_obs",
                        "T_study", "mu_oa_amp")
  expect_true(all(expected_fields %in% names(p)))
})

test_that("beta equals R0 * gamma", {
  p <- do.call(make_params, valid)
  expect_equal(p$beta, p$R0 * p$gamma)
})

test_that("n_strat equals length of f", {
  p <- do.call(make_params, valid)
  expect_equal(p$n_strat, length(valid$f))
})

test_that("strat_names defaults to 'Stratum i' when not supplied", {
  p <- do.call(make_params, valid)
  expect_equal(p$strat_names, paste0("Stratum ", 1:4))
})

test_that("risk_obs defaults to 'unknown' when not supplied", {
  p <- do.call(make_params, valid)
  expect_true(all(p$risk_obs == "unknown"))
})

test_that("custom strat_names and risk_obs are preserved", {
  args <- valid
  args$strat_names <- c("A", "B", "C", "D")
  args$risk_obs    <- c("low", "high", "low", "high")
  p <- do.call(make_params, args)
  expect_equal(p$strat_names, c("A", "B", "C", "D"))
  expect_equal(p$risk_obs,    c("low", "high", "low", "high"))
})

# ── Validation: wrong-length vectors ────────────────────────────────────────

test_that("wrong-length sus triggers error", {
  args <- valid; args$sus <- c(1.0, 1.0, 0.4)   # length 3, not 4
  expect_error(do.call(make_params, args))
})

test_that("wrong-length inf_mod triggers error", {
  args <- valid; args$inf_mod <- c(1.0, 1.0, 1.0)
  expect_error(do.call(make_params, args))
})

test_that("wrong-length seed triggers error", {
  args <- valid; args$seed <- c(5L, 0L, 0L)
  expect_error(do.call(make_params, args))
})

test_that("wrong-length p_seek_case triggers error", {
  args <- valid; args$p_seek_case <- c(0.7, 0.7, 0.7)
  expect_error(do.call(make_params, args))
})

test_that("wrong-length p_seek_control triggers error", {
  args <- valid; args$p_seek_control <- c(0.4, 0.4, 0.4)
  expect_error(do.call(make_params, args))
})

test_that("wrong-length mu_oa triggers error", {
  args <- valid; args$mu_oa <- c(1e-3, 1e-3, 1e-3)
  expect_error(do.call(make_params, args))
})

test_that("wrong-length is_vaccinated triggers error", {
  args <- valid; args$is_vaccinated <- c(FALSE, FALSE, TRUE)
  expect_error(do.call(make_params, args))
})

# ── Validation: f does not sum to 1 ─────────────────────────────────────────

test_that("fractions not summing to 1 triggers error", {
  args <- valid; args$f <- c(0.4, 0.4, 0.1, 0.2)   # sum = 1.1
  expect_error(do.call(make_params, args))
})

# ── Validation: missing vaccination groups ───────────────────────────────────

test_that("all-unvaccinated triggers error", {
  args <- valid; args$is_vaccinated <- c(FALSE, FALSE, FALSE, FALSE)
  expect_error(do.call(make_params, args))
})

test_that("all-vaccinated triggers error", {
  args <- valid; args$is_vaccinated <- c(TRUE, TRUE, TRUE, TRUE)
  expect_error(do.call(make_params, args))
})
