# R/plots.R
# Visualization functions for VE comparison simulation.
# All return ggplot objects and save to outputs/.

library(ggplot2)
library(tidyr)
library(dplyr)

.ensure_outdir <- function() {
  dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
}

# ---------------------------------------------------------------------------
# 1. plot_epidemic(sim) — SIR curves aggregated by reference / vaccinated group
# ---------------------------------------------------------------------------
plot_epidemic <- function(sim, scenario_name = "") {
  .ensure_outdir()
  ode     <- sim$ode_out
  params  <- sim$params
  N_strat <- sim$N_strat
  n       <- params$n_strat
  vax_idx <- which( params$is_vaccinated)
  ref_idx <- which(!params$is_vaccinated)
  N_ref   <- sum(N_strat[ref_idx])
  N_vax   <- sum(N_strat[vax_idx])

  agg <- function(idx, prefix) rowSums(
    as.matrix(ode[, paste0(prefix, idx), drop = FALSE])
  )

  ref_df <- data.frame(
    t     = ode$t,
    S     = agg(ref_idx, "S_"),
    I     = agg(ref_idx, "I_"),
    R     = agg(ref_idx, "Rec_"),
    group = sprintf("Reference (N=%s)", format(round(N_ref), big.mark = ","))
  )
  vax_df <- data.frame(
    t     = ode$t,
    S     = agg(vax_idx, "S_"),
    I     = agg(vax_idx, "I_"),
    R     = agg(vax_idx, "Rec_"),
    group = sprintf("Vaccinated (N=%s)", format(round(N_vax), big.mark = ","))
  )

  long <- bind_rows(ref_df, vax_df) %>%
    pivot_longer(c(S, I, R), names_to = "Compartment", values_to = "Count") %>%
    mutate(Compartment = factor(Compartment, c("S", "I", "R")))

  subtitle <- sprintf(
    "R0=%.1f  n_strat=%d  f_vax=%.2f  VE_HR(t=0)=%.2f",
    params$beta / params$gamma, n,
    sum(params$f[vax_idx]),
    1 - weighted.mean(params$sus[vax_idx], N_strat[vax_idx]) /
        weighted.mean(params$sus[ref_idx], N_strat[ref_idx])
  )

  p <- ggplot(long, aes(x = t, y = Count, colour = Compartment)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~group, scales = "free_y") +
    scale_colour_manual(values = c(S = "#2196F3", I = "#F44336", R = "#4CAF50")) +
    labs(
      title    = if (nchar(scenario_name) > 0)
                   paste("Epidemic curves —", scenario_name) else "Epidemic curves",
      subtitle = subtitle,
      x = "Day", y = "Individuals"
    ) +
    theme_bw(base_size = 12)

  fname <- file.path("outputs", paste0("epidemic_", gsub("\\W+", "_", scenario_name), ".png"))
  ggsave(fname, p, width = 9, height = 4, dpi = 150)
  message("Saved: ", fname)
  invisible(p)
}

# ---------------------------------------------------------------------------
# 2. plot_ve_comparison(results_df) — grouped bar chart, one bar per estimator.
#    All three estimators use observed data only (T_pos, T_neg, N).
#    Diamond marker shows oracle VE_true (N-weighted susceptibility ratio from ODE).
# ---------------------------------------------------------------------------
plot_ve_comparison <- function(results_df) {
  .ensure_outdir()

  sc_levels <- unique(results_df$scenario)

  long <- results_df %>%
    select(scenario, VE_RR, VE_HR, VE_TND) %>%
    pivot_longer(c(VE_RR, VE_HR, VE_TND), names_to = "Estimator", values_to = "VE") %>%
    mutate(
      Estimator = factor(Estimator, c("VE_RR", "VE_HR", "VE_TND")),
      scenario  = factor(scenario, sc_levels)
    )

  ref_df <- results_df %>%
    mutate(scenario = factor(scenario, sc_levels))

  p <- ggplot(long, aes(x = scenario, y = VE, fill = Estimator)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    geom_point(
      data        = ref_df,
      aes(x = scenario, y = VE_true),
      inherit.aes = FALSE,
      shape = 18, size = 4, colour = "black"
    ) +
    scale_fill_manual(values = c(VE_RR = "#E53935", VE_HR = "#1E88E5", VE_TND = "#43A047")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(ylim = c(-0.1, 0.8)) +
    labs(
      title    = "VE estimator comparison — observed data only",
      subtitle = "All bars: T_pos / T_neg / N  |  \u25c6 = oracle VE_true (N-weighted susceptibility ratio, not observable)",
      x = "Scenario", y = "Estimated VE", fill = "Estimator"
    ) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  fname <- "outputs/ve_comparison.png"
  ggsave(fname, p, width = 11, height = 5, dpi = 150)
  message("Saved: ", fname)
  invisible(p)
}

# ---------------------------------------------------------------------------
# 3. plot_ve_over_time(all_results) — VE_RR / VE_HR / VE_TND time series,
#    2-column facet grid (one panel per scenario).
# ---------------------------------------------------------------------------
plot_ve_over_time <- function(all_results) {
  .ensure_outdir()

  # Weekly estimates (adjusted for risk group; no uncertainty — deterministic model)
  # VE_TND_cum : cumulative OR from season start through each week
  # VE_TND_inc : incident OR for that week only (cases / controls in the window)
  # VE_HR_week : Poisson HRR, at-risk = N − cumulative infections at week start
  combined <- do.call(rbind, lapply(all_results, function(res) {
    wk <- res$ve_stat$weekly
    nm <- res$params$name
    rbind(
      data.frame(t_mid = wk$t_mid, VE = wk$VE_RR_week,
                 Estimator = "Cumul. RR",      scenario = nm),
      data.frame(t_mid = wk$t_mid, VE = wk$VE_TND_cum,
                 Estimator = "Cumul. TND OR",  scenario = nm),
      data.frame(t_mid = wk$t_mid, VE = wk$VE_TND_inc,
                 Estimator = "Weekly TND OR",  scenario = nm),
      data.frame(t_mid = wk$t_mid, VE = wk$VE_HR_week,
                 Estimator = "Weekly HRR (S)", scenario = nm)
    )
  })) %>%
    mutate(
      Estimator = factor(Estimator,
                         c("Cumul. RR", "Cumul. TND OR",
                           "Weekly TND OR", "Weekly HRR (S)")),
      scenario  = factor(scenario, levels = sapply(all_results, function(r) r$params$name))
    )

  cols <- c("Cumul. RR"      = "#FB8C00",
            "Cumul. TND OR"  = "#1E88E5",
            "Weekly TND OR"  = "#43A047",
            "Weekly HRR (S)" = "#E53935")

  p <- ggplot(combined, aes(x = t_mid, y = VE, colour = Estimator)) +
    geom_line(linewidth = 0.85, na.rm = TRUE) +
    geom_point(size = 1.6, na.rm = TRUE) +
    facet_wrap(~scenario, ncol = 2) +
    scale_colour_manual(
      values = cols,
      labels = c(
        "Cumul. RR"      = "Cumulative RR",
        "Cumul. TND OR"  = "Cumulative TND OR",
        "Weekly TND OR"  = "Weekly TND OR",
        "Weekly HRR (S)" = "Weekly HRR (N \u2212 cum. cases)"
      )
    ) +
    scale_x_continuous(breaks = seq(0, 150, by = 30)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, 0.75)) +
    labs(
      title    = "Weekly VE estimates over the epidemic season (adjusted for risk group)",
      subtitle = "Adjusted for observed risk group",
      x = "Day (week midpoint)", y = "Estimated VE", colour = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9))

  fname <- "outputs/ve_over_time.png"
  ggsave(fname, p, width = 11, height = 9, dpi = 150)
  message("Saved: ", fname)
  invisible(p)
}

# ---------------------------------------------------------------------------
# 4. plot_test_positivity(all_results) — test positivity over time by group,
#    2-column facet grid (one panel per scenario).
# ---------------------------------------------------------------------------
plot_test_positivity <- function(all_results) {
  .ensure_outdir()

  combined <- do.call(rbind, lapply(all_results, function(res) {
    ts <- res$ve$ve_ts
    rbind(
      data.frame(t = ts$t, positivity = ts$pos_vax, group = "Vaccinated",
                 scenario = res$params$name),
      data.frame(t = ts$t, positivity = ts$pos_ref,  group = "Reference",
                 scenario = res$params$name)
    )
  })) %>%
    mutate(
      group    = factor(group, c("Reference", "Vaccinated")),
      scenario = factor(scenario, levels = sapply(all_results, function(r) r$params$name))
    )

  p <- ggplot(combined, aes(x = t, y = positivity, colour = group)) +
    geom_line(linewidth = 0.85, na.rm = TRUE) +
    facet_wrap(~scenario, ncol = 2, scales = "free_x") +
    scale_colour_manual(values = c(Reference = "#E53935", Vaccinated = "#1E88E5")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title   = "Test positivity over time by vaccination status",
      subtitle = "Positivity = cumulative positive tests / cumulative total tests",
      x = "Day", y = "Test positivity", colour = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  fname <- "outputs/test_positivity.png"
  ggsave(fname, p, width = 11, height = 8, dpi = 150)
  message("Saved: ", fname)
  invisible(p)
}

# ---------------------------------------------------------------------------
# 5. plot_tnd_table(ve_result) — stacked bars showing case/control split
# ---------------------------------------------------------------------------
plot_tnd_table <- function(ve_result, scenario_name = "") {
  .ensure_outdir()

  df <- data.frame(
    Group    = c("Reference", "Reference", "Vaccinated", "Vaccinated"),
    Category = c("Cases", "Controls", "Cases", "Controls"),
    Count    = c(ve_result$cases_ref, ve_result$controls_ref,
                 ve_result$cases_vax, ve_result$controls_vax)
  ) %>%
    group_by(Group) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup() %>%
    mutate(
      Group    = factor(Group, c("Reference", "Vaccinated")),
      Category = factor(Category, c("Controls", "Cases"))
    )

  label <- sprintf("OR = %.3f\nVE_TND = %.1f%%",
                   ve_result$OR_TND, ve_result$VE_TND * 100)

  p <- ggplot(df, aes(x = Group, y = Proportion, fill = Category)) +
    geom_col(width = 0.5) +
    scale_fill_manual(values = c(Cases = "#E53935", Controls = "#90CAF9")) +
    scale_y_continuous(labels = scales::percent_format()) +
    annotate("text", x = 2.4, y = 0.5, label = label, hjust = 0, size = 4) +
    coord_cartesian(xlim = c(0.5, 3)) +
    labs(
      title    = if (nchar(scenario_name) > 0)
                   paste("TND 2\u00d72 table \u2014", scenario_name) else "TND 2\u00d72 table",
      subtitle = "Proportion of test-positives (cases) among test-seekers",
      x = NULL, y = "Proportion", fill = NULL
    ) +
    theme_bw(base_size = 12)

  fname <- file.path("outputs", paste0("tnd_table_", gsub("\\W+", "_", scenario_name), ".png"))
  ggsave(fname, p, width = 6, height = 4, dpi = 150)
  message("Saved: ", fname)
  invisible(p)
}

# ---------------------------------------------------------------------------
# 6. plot_ve_statistical(all_results) — adjusted vs unadjusted bar chart for
#    all three statistical estimators, faceted by estimator, grouped by scenario.
# ---------------------------------------------------------------------------
plot_ve_statistical <- function(all_results) {
  .ensure_outdir()

  rows <- do.call(rbind, lapply(all_results, function(res) {
    s  <- res$ve_stat
    nm <- res$params$name
    data.frame(
      scenario  = nm,
      estimator = rep(c("VE_RR", "VE_TND", "VE_HR"), each = 2),
      adjusted  = rep(c("Unadjusted", "Adjusted"), times = 3),
      VE        = c(s$VE_RR_unadj,  s$VE_RR_adj,
                    s$VE_TND_unadj, s$VE_TND_adj,
                    s$VE_HR_unadj,  s$VE_HR_adj),
      lo        = c(s$VE_RR_unadj_lo,  s$VE_RR_adj_lo,
                    s$VE_TND_unadj_lo, s$VE_TND_adj_lo,
                    s$VE_HR_unadj_lo,  s$VE_HR_adj_lo),
      hi        = c(s$VE_RR_unadj_hi,  s$VE_RR_adj_hi,
                    s$VE_TND_unadj_hi, s$VE_TND_adj_hi,
                    s$VE_HR_unadj_hi,  s$VE_HR_adj_hi),
      stringsAsFactors = FALSE
    )
  }))

  rows <- rows %>%
    mutate(
      estimator = factor(estimator, c("VE_RR", "VE_TND", "VE_HR")),
      adjusted  = factor(adjusted,  c("Unadjusted", "Adjusted")),
      scenario  = factor(scenario,  unique(sapply(all_results, function(r) r$params$name)))
    )

  p <- ggplot(rows, aes(x = scenario, y = VE, fill = adjusted,
                        ymin = lo, ymax = hi)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    geom_errorbar(position = position_dodge(width = 0.7), width = 0.25,
                  linewidth = 0.6) +
    geom_hline(yintercept = 0.6, linetype = "dashed", colour = "black", linewidth = 0.6) +
    facet_wrap(~estimator, ncol = 3) +
    scale_fill_manual(values = c(Unadjusted = "#90CAF9", Adjusted = "#1565C0")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(ylim = c(-0.1, 0.75)) +
    labs(
      title    = "Statistical VE estimates: unadjusted vs adjusted for observed risk group",
      subtitle = "Dashed line = true VE_s = 60%  |  Error bars = 95% Wald CI",
      x = NULL, y = "Estimated VE", fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "top")

  fname <- "outputs/ve_statistical.png"
  ggsave(fname, p, width = 13, height = 5, dpi = 150)
  message("Saved: ", fname)
  invisible(p)
}

# ---------------------------------------------------------------------------
# 7. plot_ve_hr_time(all_results) — time-varying VE_HR per bin with CIs,
#    faceted by scenario (shows whether VE changes over the epidemic course).
# ---------------------------------------------------------------------------
plot_ve_hr_time <- function(all_results) {
  .ensure_outdir()

  combined <- do.call(rbind, lapply(all_results, function(res) {
    tv <- res$ve_stat$tv
    data.frame(t_mid = tv$t_mid, VE = tv$VE_HR, lo = tv$HR_lo, hi = tv$HR_hi,
               scenario = res$params$name)
  })) %>%
    mutate(scenario = factor(scenario,
                             levels = sapply(all_results, function(r) r$params$name)))

  p <- ggplot(combined, aes(x = t_mid, y = VE, ymin = lo, ymax = hi)) +
    geom_hline(yintercept = 0.6, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
    geom_hline(yintercept = 0,   linetype = "solid",  colour = "grey80", linewidth = 0.4) +
    geom_ribbon(fill = "#BBDEFB", alpha = 0.6, na.rm = TRUE) +
    geom_line(colour = "#1565C0", linewidth = 0.9, na.rm = TRUE) +
    geom_point(colour = "#1565C0", size = 2, na.rm = TRUE) +
    facet_wrap(~scenario, ncol = 2) +
    scale_x_continuous(breaks = seq(0, 150, by = 30)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(ylim = c(-0.1, 0.75)) +
    labs(
      title    = "Time-varying VE_HR from piecewise Poisson counting process",
      subtitle = "Adjusted for observed risk group  |  At-risk = daily total tests  |  Dashed = true VE_s = 60%",
      x = "Day (bin midpoint)", y = "VE_HR"
    ) +
    theme_bw(base_size = 12)

  fname <- "outputs/ve_hr_time.png"
  ggsave(fname, p, width = 11, height = 8, dpi = 150)
  message("Saved: ", fname)
  invisible(p)
}

# ---------------------------------------------------------------------------
# 8. plot_hazard_rates(all_results) — instantaneous group hazard rates from ODE,
#    reference vs vaccinated, faceted by scenario.  A second row per scenario
#    shows the hazard ratio HR(t) = h_vax / h_ref = 1 - VE_HR(t).
# ---------------------------------------------------------------------------
plot_hazard_rates <- function(all_results) {
  .ensure_outdir()

  hz <- do.call(rbind, lapply(all_results, function(res) {
    ts  <- res$ve$ve_ts
    eps <- 1e-10
    rbind(
      data.frame(t = ts$t, value = ts$h_ref, series = "Reference",  type = "Hazard rate",
                 scenario = res$params$name),
      data.frame(t = ts$t, value = ts$h_vax, series = "Vaccinated", type = "Hazard rate",
                 scenario = res$params$name),
      data.frame(t = ts$t,
                 value   = ifelse(ts$h_ref > eps, ts$h_vax / ts$h_ref, NA_real_),
                 series  = "HR (vax/ref)", type = "Hazard ratio",
                 scenario = res$params$name)
    )
  })) %>%
    mutate(
      series   = factor(series, c("Reference", "Vaccinated", "HR (vax/ref)")),
      type     = factor(type,   c("Hazard rate", "Hazard ratio")),
      scenario = factor(scenario, levels = sapply(all_results, function(r) r$params$name))
    )

  cols <- c(Reference = "#E53935", Vaccinated = "#1E88E5", `HR (vax/ref)` = "#7B1FA2")

  p <- ggplot(hz, aes(x = t, y = value, colour = series)) +
    geom_line(linewidth = 0.85, na.rm = TRUE) +
    geom_hline(
      data = data.frame(type = factor("Hazard ratio", c("Hazard rate", "Hazard ratio")),
                        yint = 0.4),
      aes(yintercept = yint), linetype = "dashed", colour = "grey50", linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    facet_grid(type ~ scenario, scales = "free_y") +
    scale_colour_manual(values = cols) +
    scale_x_continuous(breaks = seq(0, 150, by = 30)) +
    labs(
      title    = "Instantaneous hazard rates and hazard ratio over time [oracle / theoretical]",
      subtitle = "ODE internal state: S-weighted mean of h_i(t)  |  Not observable in practice  |  Dashed = true HR = 1 \u2212 VE_s = 0.4",
      x = "Day", y = NULL, colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text.x = element_text(size = 8))

  fname <- "outputs/hazard_rates.png"
  ggsave(fname, p, width = 14, height = 6, dpi = 150)
  message("Saved: ", fname)
  invisible(p)
}

# ---------------------------------------------------------------------------
# 9. plot_hrr_comparison(all_results) — theoretical HRR from ODE vs observed
#    weekly HRR from Poisson model (at-risk = N − cumulative positive tests).
# ---------------------------------------------------------------------------
plot_hrr_comparison <- function(all_results) {
  .ensure_outdir()

  theory <- do.call(rbind, lapply(all_results, function(res) {
    ts <- res$ve$ve_ts
    data.frame(t = ts$t, VE = ts$VE_HR, series = "Theoretical (ODE)",
               scenario = res$params$name)
  }))

  observed <- do.call(rbind, lapply(all_results, function(res) {
    wk <- res$ve_stat$weekly
    data.frame(t = wk$t_mid, VE = wk$VE_HR_week, series = "Observed (weekly)",
               scenario = res$params$name)
  }))

  combined <- rbind(theory, observed) %>%
    mutate(
      series   = factor(series, c("Theoretical (ODE)", "Observed (weekly)")),
      scenario = factor(scenario, levels = sapply(all_results, function(r) r$params$name))
    )

  cols <- c("Theoretical (ODE)" = "#1E88E5", "Observed (weekly)" = "#E53935")

  p <- ggplot(combined, aes(x = t, y = VE, colour = series)) +
    geom_line(data  = \(d) d[d$series == "Theoretical (ODE)", ],
              linewidth = 0.8, na.rm = TRUE) +
    geom_line(data  = \(d) d[d$series == "Observed (weekly)", ],
              linewidth = 0.6, linetype = "dashed", na.rm = TRUE) +
    geom_point(data = \(d) d[d$series == "Observed (weekly)", ],
               size = 1.8, na.rm = TRUE) +
    facet_wrap(~scenario, ncol = 2) +
    scale_colour_manual(values = cols) +
    scale_x_continuous(breaks = seq(0, 150, by = 30)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, 0.75)) +
    labs(
      title    = "Theoretical vs observed hazard rate ratio",
      subtitle = "Theoretical = ODE h_vax/h_ref  |  Observed = Poisson, at-risk = N \u2212 cum. positive tests",
      x = "Day", y = "VE (= 1 \u2212 HRR)", colour = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  fname <- "outputs/hrr_comparison.png"
  ggsave(fname, p, width = 11, height = 9, dpi = 150)
  message("Saved: ", fname)
  invisible(p)
}
