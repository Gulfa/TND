# TND — Vaccine Effectiveness Estimator Comparison

A deterministic simulation study comparing three common vaccine effectiveness (VE) estimators across epidemiological scenarios where they are expected to agree or diverge.

## Estimators

| Estimator | Definition | Based on |
|-----------|-----------|----------|
| **VE_RR** | 1 − cumulative attack-rate ratio | Full population incidence |
| **VE_HR** | 1 − hazard ratio | Instantaneous infection rates |
| **VE_TND** | 1 − odds ratio from test-negative design | Tested individuals only |

## Model

An *n*-stratum deterministic SIR model implemented with [odin](https://mrc-ide.github.io/odin/). Each stratum carries its own susceptibility, infectiousness, test-seeking behaviour, and background (non-target) illness rate. A parallel accumulator tracks TND-style test counts (positive and negative) with an optional seasonal background rate.

Force of infection is frequency-dependent with homogeneous mixing across strata.

## Scenarios

| # | Name | Key feature | Expected bias |
|---|------|-------------|---------------|
| 1 | R0 = 1.3 | Low epidemic, uniform groups | All estimators ≈ true VE |
| 2 | R0 = 1.5 | Moderate epidemic | VE_RR < VE_HR (susceptible depletion) |
| 3 | Unobserved heterogeneity | Within-group frailty (factor-7 spread) | HR non-constant; VE_HR wanes spuriously |
| 4 | Indication bias | High-risk disproportionately vaccinated | All estimators biased; adjustment partially corrects |
| 5 | Care-seeking bias | Unvaccinated less likely to seek care | VE_HR biased; VE_TND self-corrects |
| 6 | Other-illness confounding | Indication bias + high-risk have 3× background illness | VE_TND_unadj biased opposite direction to VE_RR |
| 7 | VE on other illness | Vaccine halves background illness rate | VE_TND severely overestimates VE |
| 8 | Differential seasonal | Reference group strongly seasonal, vaccinated flat | Cumulative TND OR ≠ RR; week-adjustment partially corrects |

## File Structure

```
TND/
├── R/
│   ├── model.R       # odin SIR model (compiled once at source time)
│   ├── simulate.R    # run_simulation(params) → list(ode_out, ...)
│   ├── measures.R    # compute_ve() + compute_ve_statistical()
│   └── plots.R       # plot_*() functions
├── scenarios.R       # make_flu_scenario(), calibration, 8 scenarios
├── main.R            # entry point: source all, run, save outputs/
├── run_in_docker.sh  # run inside FHI R Docker container
└── tests/            # testthat unit tests
```

## Running

### In Docker (recommended)

```bash
bash run_in_docker.sh          # uses container named "fhi" by default
bash run_in_docker.sh my_cont  # specify container name
```

The container must have the host directory mounted at `/fhi` and the following R packages installed: `odin`, `dplyr`, `tidyr`, `ggplot2`, `scales`, `purrr`.

Install packages (first time):

```bash
docker exec fhi Rscript -e \
  "install.packages(c('odin','dplyr','tidyr','ggplot2','scales','purrr'), repos='https://cloud.r-project.org')"
```

### Directly with Rscript

```bash
Rscript main.R
```

Outputs are saved to `outputs/` (one PNG per scenario + combined plots).

## Key Results

All scenarios use `gamma = 0.5` (2-day generation time), `T_study = 150` days, `VE_s = 0.6` (true susceptibility VE), and `p_seek_case = 0.2`. The background illness rate is calibrated to 35% test positivity at R0 = 1.5.

The simulation demonstrates:

- **VE_RR** always underestimates VE_HR when the epidemic is large (attack-rate ratio ≠ hazard ratio under susceptible depletion).
- **VE_TND** equals VE_RR when test-seeking rates are proportional between groups, but diverges when the background illness rate differs between vaccinated and unvaccinated (scenarios 6–8).
- **Risk-group adjustment** corrects indication bias (scenarios 4, 6) but cannot fix differential test-seeking (scenario 5) or VE on background illness (scenario 7).
- **Week-adjustment** in TND only helps when there is group-differential time-varying confounding in the background rate (scenario 8); it has no effect when the seasonal pattern is identical across groups.
