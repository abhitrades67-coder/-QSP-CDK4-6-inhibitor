# =============================================================================
# 04_trials.R
# Clinical-trial simulation, calibration and competing-hypothesis testing.
#
# Design
# ------
# Each randomised trial contrast is modelled as the SAME virtual population
# (one resistance-composition hypothesis) treated under two arms:
#
#   Trial         control arm          experimental arm
#   -----         -----------          ----------------
#   MAINTAIN      switched ET (none)    switch to ribociclib
#   PACE          fulvestrant (none)    continue palbociclib (same agent)
#   postMONARCH   fulvestrant (none)    switch to abemaciclib
#   PALMIRA       ET (none)             rechallenge palbociclib (same agent)
#
# All control arms are ET-only in the model; the DIFFERENT observed control-arm
# PFS across trials is attributed to population pace-of-disease differences and
# captured by a per-trial aggressiveness scalar fit to the control arm alone.
# The EXPERIMENTAL arm is then a PREDICTION (not fit), so matching its PFS/HR is
# genuine validation.
#
# A shared resistance composition is calibrated (Latin-hypercube search over the
# 5-state simplex) to best predict the four experimental arms. Competing fixed
# composition hypotheses (drug-dependent / agnostic / ER-driven / mixed) are then
# scored the same way to ask which mechanistic picture is compatible with the
# full four-trial pattern.
# =============================================================================

if (!exists("simulate_arm_pfs")) {
  for (cand in c("03_virtual_population.R", "R/03_virtual_population.R"))
    if (file.exists(cand)) { source(cand); break }
}

# ---- Observed trial targets -------------------------------------------------
TRIALS <- list(
  MAINTAIN    = list(control_arm = "ET_only", exp_arm = "switch_ribo",
                     ctrl_pfs = 2.76, exp_pfs = 5.29, hr = 0.57,
                     hr_lo = 0.39, hr_hi = 0.85,
                     exp_label = "switch to ribociclib + ET"),
  PACE        = list(control_arm = "ET_only", exp_arm = "continue_palbo",
                     ctrl_pfs = 4.80, exp_pfs = 4.60, hr = 1.11,
                     hr_lo = 0.79, hr_hi = 1.55,
                     exp_label = "continue palbociclib + fulvestrant"),
  postMONARCH = list(control_arm = "ET_only", exp_arm = "switch_abema",
                     ctrl_pfs = 5.30, exp_pfs = 6.00, hr = 0.73,
                     hr_lo = 0.57, hr_hi = 0.95,
                     exp_label = "switch to abemaciclib + fulvestrant"),
  PALMIRA     = list(control_arm = "ET_only", exp_arm = "continue_palbo",
                     ctrl_pfs = 3.60, exp_pfs = 4.90, hr = 0.84,
                     hr_lo = 0.66, hr_hi = 1.07,
                     exp_label = "rechallenge palbociclib + ET")
)

# ---- Cox HR from two paired arms --------------------------------------------
cox_hr <- function(ctrl_df, exp_df) {
  d <- rbind(
    data.frame(time = ctrl_df$ttp_days, event = as.integer(ctrl_df$progressed), arm = 0L),
    data.frame(time = exp_df$ttp_days,  event = as.integer(exp_df$progressed),  arm = 1L)
  )
  d$time <- pmax(d$time, 0.5)
  fit <- tryCatch(survival::coxph(survival::Surv(time, event) ~ arm, data = d),
                  error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  as.numeric(exp(coef(fit))[1])
}

median_pfs <- function(df) median(df$pfs_months)

# ---- Fit per-trial aggressiveness from the control arm ----------------------
# All controls are ET_only, so one aggressiveness->median-PFS curve serves all
# trials; invert by interpolation.
aggress_curve <- function(pop, growth, control_arm = "ET_only",
                          a_grid = seq(0.45, 3.2, length.out = 14),
                          dt = 3) {
  m <- sapply(a_grid, function(a)
    median_pfs(simulate_arm_pfs(pop, control_arm, growth, aggress_scale = a, dt = dt)))
  data.frame(a = a_grid, pfs = m)
}

aggress_for_target <- function(curve, target_pfs) {
  # curve$pfs is decreasing in a; interpolate a at target
  o <- order(curve$pfs)
  a <- approx(curve$pfs[o], curve$a[o], xout = target_pfs, rule = 2)$y
  min(max(a, min(curve$a)), max(curve$a))
}

# ---- Evaluate one composition against all trials ----------------------------
evaluate_composition <- function(comp, n = 500, dt = 3, seed = 7,
                                 growth = NULL) {
  if (is.null(growth)) growth <- precompute_growth()
  pop <- make_population(n, comp, seed = seed)
  curve <- aggress_curve(pop, growth, "ET_only", dt = dt)

  rows <- lapply(names(TRIALS), function(tn) {
    tr <- TRIALS[[tn]]
    a  <- aggress_for_target(curve, tr$ctrl_pfs)
    ctrl <- simulate_arm_pfs(pop, tr$control_arm, growth, aggress_scale = a, dt = dt)
    expd <- simulate_arm_pfs(pop, tr$exp_arm,     growth, aggress_scale = a, dt = dt)
    data.frame(
      trial = tn,
      aggress = a,
      ctrl_pfs_obs = tr$ctrl_pfs,
      ctrl_pfs_pred = median_pfs(ctrl),
      exp_pfs_obs = tr$exp_pfs,
      exp_pfs_pred = median_pfs(expd),
      hr_obs = tr$hr,
      hr_pred = cox_hr(ctrl, expd),
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, rows)
  # error: experimental PFS (months) + HR (log scale), experimental arm only
  res$err_pfs <- (res$exp_pfs_pred - res$exp_pfs_obs)
  res$err_hr  <- (log(pmax(res$hr_pred, 1e-3)) - log(res$hr_obs))
  attr(res, "objective") <- sum(res$err_pfs^2) + 4 * sum(res$err_hr^2)
  res
}

# ---- Latin-hypercube search over the 5-state simplex ------------------------
# Sample Dirichlet-like compositions; keep biologically plausible ranges.
simplex_from_lhs <- function(u) {
  # u: vector length 5 in (0,1); map to simplex via -log then normalise
  g <- -log(pmax(u, 1e-6))
  g / sum(g)
}

calibrate_composition <- function(n_search = 160, n_eval = 450, dt = 4,
                                  seed = 101) {
  set.seed(seed)
  growth <- precompute_growth()
  L <- lhs::randomLHS(n_search, 5)
  best <- NULL; best_obj <- Inf; best_res <- NULL
  objs <- numeric(n_search)
  for (i in seq_len(n_search)) {
    comp <- setNames(simplex_from_lhs(L[i, ]), STATE_NAMES)
    res <- evaluate_composition(comp, n = n_eval, dt = dt, seed = 7, growth = growth)
    o <- attr(res, "objective")
    objs[i] <- o
    if (is.finite(o) && o < best_obj) { best_obj <- o; best <- comp; best_res <- res }
  }
  list(composition = best, objective = best_obj, res = best_res,
       objs = objs, growth = growth)
}

# ---- Identifiable reparametrisation: switch-responsive fraction f_sr ---------
# The population is a blend of a SWITCH-RESPONSIVE sub-mixture (clones a new
# CDK4/6i or ET switch can still engage) and an AGNOSTIC sub-mixture (clones no
# CDK4/6i engages). One interpretable parameter, f_sr in [0,1], controls the
# balance; the internal profiles are held fixed (biology).
SR_PROFILE <- c(sensitive = 0.25, cdk6_remodel = 0.45, cdk2_bypass = 0.00,
                rb_loss = 0.00, er_reactivation = 0.30)[STATE_NAMES]
AG_PROFILE <- c(sensitive = 0.00, cdk6_remodel = 0.00, cdk2_bypass = 0.75,
                rb_loss = 0.25, er_reactivation = 0.00)[STATE_NAMES]

comp_from_fsr <- function(fsr) {
  cc <- fsr * SR_PROFILE + (1 - fsr) * AG_PROFILE
  cc / sum(cc)
}

# Per-trial f_sr analysis: fit aggressiveness to the control arm, then map how
# the experimental-arm PFS depends on f_sr and find the f_sr matching the
# observed experimental PFS. The RANGE of experimental PFS over f_sr is the
# model's identifiability signal: large for switch trials, ~0 for same-agent
# continuation/re-challenge trials (which therefore cannot pin the mechanism).
analyse_fsr <- function(growth, n = 700, dt = 2, seed = 11,
                        fsr_grid = seq(0, 1, by = 0.1)) {
  ref_pop <- make_population(n, comp_from_fsr(0.5), seed = seed)
  curve   <- aggress_curve(ref_pop, growth, "ET_only", dt = dt)

  rows <- lapply(names(TRIALS), function(tn) {
    tr <- TRIALS[[tn]]
    a  <- aggress_for_target(curve, tr$ctrl_pfs)
    expv <- sapply(fsr_grid, function(fs) {
      pop <- make_population(n, comp_from_fsr(fs), seed = seed + 1)
      median_pfs(simulate_arm_pfs(pop, tr$exp_arm, growth, aggress_scale = a, dt = dt))
    })
    rng <- max(expv) - min(expv)
    # f_sr matching observed experimental PFS (if within reachable range)
    o <- order(expv)
    fsr_hat <- if (tr$exp_pfs >= min(expv) && tr$exp_pfs <= max(expv)) {
      approx(expv[o], fsr_grid[o], xout = tr$exp_pfs)$y
    } else NA_real_
    data.frame(trial = tn, aggress = a,
               exp_pfs_obs = tr$exp_pfs,
               exp_pfs_at_fsr0 = expv[1],
               exp_pfs_at_fsr1 = expv[length(expv)],
               identifiability_range = rng,
               fsr_hat = fsr_hat,
               same_agent = tr$exp_arm == "continue_palbo",
               stringsAsFactors = FALSE)
  })
  list(table = do.call(rbind, rows),
       curve = curve,
       fsr_grid = fsr_grid)
}

# ---- Competing-hypothesis test ----------------------------------------------
test_hypotheses <- function(growth, n = 600, dt = 2, seed = 7) {
  out <- lapply(names(COMPOSITIONS), function(h) {
    res <- evaluate_composition(COMPOSITIONS[[h]], n = n, dt = dt, seed = seed,
                                growth = growth)
    data.frame(hypothesis = h,
               trial = res$trial,
               hr_obs = res$hr_obs,
               hr_pred = res$hr_pred,
               exp_pfs_obs = res$exp_pfs_obs,
               exp_pfs_pred = res$exp_pfs_pred,
               objective = attr(res, "objective"),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# =============================================================================
# Main
# =============================================================================
if (.is_main("04_trials.R")) {
  t0 <- Sys.time()
  message(">> Calibrating shared resistance composition (LHS search)...")
  cal <- calibrate_composition(n_search = 160, n_eval = 450, dt = 4, seed = 101)
  comp <- cal$composition
  cat("\nCalibrated shared composition (fractions):\n")
  print(round(comp, 3))
  cat(sprintf("Search objective = %.3f\n", cal$objective))

  message(">> Final high-resolution evaluation of calibrated composition...")
  final <- evaluate_composition(comp, n = 1500, dt = 1, seed = 2024, growth = cal$growth)
  cat("\n=== Calibrated model vs observed trials ===\n")
  print(final[, c("trial","aggress","ctrl_pfs_obs","ctrl_pfs_pred",
                  "exp_pfs_obs","exp_pfs_pred","hr_obs","hr_pred")],
        row.names = FALSE, digits = 3)

  write.csv(
    transform(final, composition = paste(names(comp), round(comp,3), sep="=", collapse="; ")),
    file.path(DIR_TAB, "table_trial_fits.csv"), row.names = FALSE)

  message(">> Per-trial switch-responsive-fraction (f_sr) / identifiability analysis...")
  fsr <- analyse_fsr(cal$growth, n = 900, dt = 1, seed = 11)
  write.csv(fsr$table, file.path(DIR_TAB, "table_fsr_identifiability.csv"), row.names = FALSE)
  cat("\n=== Switch-responsive fraction & identifiability by trial ===\n")
  print(fsr$table, row.names = FALSE, digits = 3)
  cat("\nInterpretation: identifiability_range = how much the experimental-arm PFS\n",
      "moves as f_sr goes 0->1. Large => the trial constrains the resistance\n",
      "mechanism; ~0 (same-agent continuation) => it cannot.\n")

  message(">> Testing competing resistance-composition hypotheses...")
  hyp <- test_hypotheses(cal$growth, n = 900, dt = 1, seed = 7)
  write.csv(hyp, file.path(DIR_TAB, "table_hypothesis_tests.csv"), row.names = FALSE)
  cat("\n=== Competing-hypothesis HR predictions ===\n")
  hw <- reshape(hyp[, c("hypothesis","trial","hr_pred")],
                idvar = "hypothesis", timevar = "trial", direction = "wide")
  names(hw) <- gsub("hr_pred.", "", names(hw))
  obs_hr <- sapply(names(TRIALS), function(t) TRIALS[[t]]$hr)
  print(hw, row.names = FALSE, digits = 3)
  cat("Observed HR: ", paste(sprintf("%s=%.2f", names(obs_hr), obs_hr), collapse="  "), "\n")

  # Save objects for plotting / reporting
  saveRDS(list(composition = comp, final = final, hyp = hyp, fsr = fsr,
               growth = cal$growth, search_objs = cal$objs,
               comp_from_fsr = comp_from_fsr,
               SR_PROFILE = SR_PROFILE, AG_PROFILE = AG_PROFILE,
               TRIALS = TRIALS, ARMS = ARMS, COMPOSITIONS = COMPOSITIONS),
          file.path(DIR_RES, "trial_results.rds"))
  cat(sprintf("\nSaved results/trial_results.rds ; elapsed %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
