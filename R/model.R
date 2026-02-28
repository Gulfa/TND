# R/model.R
# n-stratum SIR odin model with cumulative incidence and TND testing trackers.
# Each stratum i has its own susceptibility (sus[i]), infectiousness (inf_mod[i]),
# and testing behaviour (p_seek_case[i], p_seek_control[i], mu_oa[i]).
#
# odin 1.x array convention:
#   LHS uses [] (implicit loop over all indices), RHS uses [i] (loop variable).
#   e.g.  deriv(S[]) <- -sus[i] * lambda * S[i]
#         I_wt[]     <- inf_mod[i] * I[i]         (intermediate equation)

odin_sir <- odin::odin({
  # ---- Dimensions ----
  n_strat <- user(integer = TRUE)
  dim(sus)           <- n_strat
  dim(inf_mod)       <- n_strat
  dim(N_strat)       <- n_strat
  dim(seed)          <- n_strat
  dim(p_seek_case)   <- n_strat
  dim(p_seek_control)<- n_strat
  dim(mu_oa)         <- n_strat
  dim(mu_oa_amp)     <- n_strat
  dim(S)             <- n_strat
  dim(I)             <- n_strat
  dim(Rec)           <- n_strat
  dim(CI)            <- n_strat
  dim(T_pos)         <- n_strat
  dim(T_neg)         <- n_strat
  dim(I_wt)          <- n_strat
  dim(h)             <- n_strat

  # ---- Force of infection (frequency-dependent, homogeneous mixing) ----
  I_wt[]  <- inf_mod[i] * I[i]     # weighted infectious per stratum
  I_eff   <- sum(I_wt[])
  lambda  <- beta * I_eff / N

  # ---- Per-stratum instantaneous hazard and total hazard ----
  # h[i](t) = sus[i] * lambda(t): individual risk rate for stratum i
  h[]   <- sus[i] * lambda
  output(h[])    <- TRUE
  output(lambda) <- TRUE

  # ---- Stratum dynamics ----
  deriv(S[])   <- -sus[i] * lambda * S[i]
  deriv(I[])   <-  sus[i] * lambda * S[i] - gamma * I[i]
  deriv(Rec[]) <-  gamma * I[i]
  deriv(CI[])  <-  sus[i] * lambda * S[i]

  # ---- TND testing accumulators ----
  # T_pos[i]: cumulative positive tests — infected individuals who sought care
  # T_neg[i]: cumulative negative tests — other-illness episodes who sought care
  #           Rate varies seasonally: mu_oa[i] * (1 + mu_oa_amp * sin(pi*t/T_study))
  #           peaking at mid-season (t = T_study/2).  mu_oa_amp=0 → constant rate.
  deriv(T_pos[]) <-  sus[i] * lambda * S[i] * p_seek_case[i]
  deriv(T_neg[]) <-  N_strat[i] * mu_oa[i] * (1 + mu_oa_amp[i] * sin(3.14159265 * t / T_study)) * p_seek_control[i]

  # ---- Initial conditions ----
  initial(S[])     <- N_strat[i] - seed[i]
  initial(I[])     <- seed[i]
  initial(Rec[])   <- 0
  initial(CI[])    <- 0
  initial(T_pos[]) <- 0
  initial(T_neg[]) <- 0

  # ---- Scalar parameters ----
  beta      <- user()
  gamma     <- user()
  N         <- user()
  T_study        <- user()   # study duration in days (used in seasonal formula)

  # ---- Array parameters ----
  sus[]            <- user()
  inf_mod[]        <- user()
  N_strat[]        <- user()
  seed[]           <- user()
  p_seek_case[]    <- user()
  p_seek_control[] <- user()
  mu_oa[]          <- user()
  mu_oa_amp[]      <- user()
}, verbose = FALSE)
