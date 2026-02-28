# test-finalsize.R
# Compares ODE cumulative attack rates to the analytical epidemic final size.
#
# For an n-stratum SIR with frequency-dependent, homogeneous mixing and
# per-stratum susceptibility multipliers sus[i], the final-size equations are:
#
#   z_i = 1 - exp(-sus[i] * R0 * sum_j(f[j] * z[j]))           [1]
#
# where z_i is the fraction of stratum i ultimately infected.
# This is solved numerically below via fixed-point iteration.
#
# Additionally, if the `finalsize` package is installed, its final_size()
# function is used as an independent cross-check.
#
# Tolerance: ±0.5 pp (absolute), reflecting small discrepancies from
# (a) finite ODE integration horizon (T_study = 150 days may not fully exhaust
#     the epidemic) and (b) seeded initial infections (seed > 0 are already
#     counted in CI but shift the fixed-point slightly).

# ── Helpers ──────────────────────────────────────────────────────────────────

# Solve final-size equations [1] by fixed-point iteration.
# Returns z: numeric vector of length n (fraction infected per stratum).
.solve_final_size <- function(R0, sus, f, tol = 1e-12, max_iter = 10000L) {
  n <- length(f)
  z <- rep(0.5, n)
  for (iter in seq_len(max_iter)) {
    force_of_inf <- R0 * sum(f * z)
    z_new        <- 1 - exp(-sus * force_of_inf)
    if (max(abs(z_new - z)) < tol) return(z_new)
    z <- z_new
  }
  warning("final_size iteration did not converge")
  z
}

# Extract per-stratum attack rate from an ODE result (fraction infected).
.ode_attack_rate <- function(sim) {
  n       <- sim$params$n_strat
  CI_last <- vapply(seq_len(n),
                    function(i) tail(sim$ode_out[[paste0("CI_", i)]], 1),
                    numeric(1))
  CI_last / sim$N_strat
}

# ── Fixed-point final-size vs ODE ────────────────────────────────────────────
# Use scenarios with uniform infectiousness (inf_mod = 1) and no mu_oa seeding
# artefact. Scenarios s1 and s2 have 4 strata and uniform parameters.

test_that("ODE AR per stratum matches fixed-point final size in s1 (R0=1.3)", {
  sim  <- all_results$s1$sim
  p    <- sim$params
  z_fp <- .solve_final_size(p$R0, p$sus, p$f)
  z_ode <- .ode_attack_rate(sim)
  expect_equal(z_ode, z_fp, tolerance = 5e-3,
               label = "per-stratum AR: ODE vs fixed-point, s1")
})

test_that("ODE AR per stratum matches fixed-point final size in s2 (R0=1.5)", {
  sim  <- all_results$s2$sim
  p    <- sim$params
  z_fp <- .solve_final_size(p$R0, p$sus, p$f)
  z_ode <- .ode_attack_rate(sim)
  expect_equal(z_ode, z_fp, tolerance = 5e-3,
               label = "per-stratum AR: ODE vs fixed-point, s2")
})

test_that("ODE total AR matches fixed-point final size in s3 (8-stratum)", {
  sim   <- all_results$s3$sim
  p     <- sim$params
  z_fp  <- .solve_final_size(p$R0, p$sus, p$f)
  z_ode <- .ode_attack_rate(sim)
  # Compare total (population-weighted) AR rather than per-stratum,
  # as small per-stratum differences can arise from the within-group seeds.
  expect_equal(sum(p$f * z_ode), sum(p$f * z_fp), tolerance = 5e-3,
               label = "total AR: ODE vs fixed-point, s3")
})

test_that("fixed-point solution is self-consistent (residual < 1e-10)", {
  for (nm in c("s1", "s2")) {
    p    <- all_results[[nm]]$params
    z    <- .solve_final_size(p$R0, p$sus, p$f)
    foi  <- p$R0 * sum(p$f * z)
    resid <- max(abs(z - (1 - exp(-p$sus * foi))))
    expect_lt(resid, 1e-10, label = paste("fixed-point residual in", nm))
  }
})

# ── Cross-check with the `finalsize` package ─────────────────────────────────
# finalsize::final_size() solves the same fixed-point equations.
# Homogeneous mixing maps to contact_matrix = matrix(1, n, n), r0 = R0.
# finalsize validates max(Re(eigen(contact_matrix * demography_vector))) == 1;
# for the all-ones matrix and sum(f) == 1, this equals sum(f) = 1 exactly.

test_that("finalsize package gives same per-stratum AR as fixed-point (s1, s2)", {
  skip_if_not_installed("finalsize")

  for (nm in c("s1", "s2")) {
    p <- all_results[[nm]]$params
    n <- p$n_strat

    # Build NGM inputs for finalsize::final_size()
    #
    # Our equation:  z_i = 1 - exp(-sus[i]*R0*sum_j(f[j]*z[j]))
    # finalsize:     z_i = 1 - exp(-r0*sus[i]*sum_j(C[i,j]*f[j]*z[j]))
    #
    # Homogeneous mixing: C[i,j] = 1 for all i,j, r0 = R0.
    # finalsize checks max(Re(eigen(C * f))) == 1:
    #   C * f (recycling) = outer(f, rep(1,n)), dominant eigenvalue = sum(f) = 1
    contact_mat <- matrix(1, nrow = n, ncol = n)
    r0_fs       <- p$R0

    fs <- finalsize::final_size(
      r0                = r0_fs,
      contact_matrix    = contact_mat,
      demography_vector = p$f,
      susceptibility    = matrix(p$sus, ncol = 1),
      p_susceptibility  = matrix(1,     nrow = n, ncol = 1)
    )

    z_pkg <- fs$p_infected   # fraction infected per group
    z_fp  <- .solve_final_size(p$R0, p$sus, p$f)

    expect_equal(z_pkg, z_fp, tolerance = 1e-4,
                 label = paste("finalsize pkg vs fixed-point, per-stratum, in", nm))
  }
})
