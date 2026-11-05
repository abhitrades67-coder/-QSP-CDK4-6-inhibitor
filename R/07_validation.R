# =============================================================================
# 07_validation.R
# Methodological completion: (A) model verification, (B) leave-one-trial-out
# cross-validation (genuine out-of-sample prediction), (C) bootstrap prediction
# intervals on predicted HR / median PFS, (D) formal validation-criteria table.
# =============================================================================

if (!exists("simulate_arm_pfs")) {
  for (f in c("00_config.R","01_model.R","03_virtual_population.R","04_trials.R"))
    for (cand in c(f, file.path("R", f))) if (file.exists(cand)) { source(cand); break }
}
suppressWarnings(RNGkind())

# =============================================================================
# (A) MODEL VERIFICATION
# =============================================================================
verify_model <- function() {
  checks <- list()
  add <- function(name, pass, detail) checks[[length(checks)+1]] <<-
    data.frame(check = name, pass = pass, detail = detail, stringsAsFactors = FALSE)

  # 1. closed-form net rate == late-time ODE growth rate
  maxdiff <- 0
  for (s in STATE_NAMES) for (rx in list(c("none","fulvestrant"), c("abemaciclib","fulvestrant"))) {
    sim <- simulate_tumour(s, rx[1], rx[2], horizon = 400, dt = 1)
    n <- nrow(sim); ode_rate <- (log(sim$V[n]) - log(sim$V[n-100]))/100
    cf <- net_growth_rate(s, rx[1], rx[2])$r_eff
    maxdiff <- max(maxdiff, abs(ode_rate - cf))
  }
  add("ODE vs closed-form net growth rate", maxdiff < 1e-6,
      sprintf("max |diff| = %.2e (tol 1e-6)", maxdiff))

  # 2. ER ODE converges to analytic er_steady
  sim <- simulate_tumour("sensitive","none","fulvestrant", horizon = 80, dt = 1)
  p <- build_context("sensitive","none","fulvestrant")
  erdiff <- abs(tail(sim$ER,1) - er_steady(p))
  add("ER ODE reaches analytic steady state", erdiff < 1e-4,
      sprintf("|diff| = %.2e", erdiff))

  # 3. fixed-point is seed-independent (no hidden hysteresis)
  p <- build_context("cdk6_remodel","abemaciclib","fulvestrant"); ER <- er_steady(p)
  solve_e2f <- function(seed) {
    E2F <- seed
    for (i in 1:300){CycE<-p$kE_syn*E2F;a2<-p$cdk2_cap*(CycE+p$cdk2_basal)
      a46<-p$cdk46_cap*cyclinD_level(p,ER)*(1-p$p16)*(1-cdk46_inhibition(p))
      Krb<-p$kp46*a46+p$kp2*a2;pRB<-Krb/(Krb+p$kdp_RB)
      E2F<-0.5*E2F+0.5*hill(pRB,p$E2F_K,p$E2F_n)}
    E2F
  }
  es <- sapply(c(0.02,0.5,0.98), solve_e2f)
  add("G1/S fixed point is seed-independent", diff(range(es)) < 1e-4,
      sprintf("range over seeds = %.2e", diff(range(es))))

  # 4. RB-loss clone is CDK4/6i-agnostic (no CDK4/6i changes its rate)
  base <- net_growth_rate("rb_loss","none","fulvestrant", prior_cdk46i="palbociclib")$r_eff
  dmax <- max(sapply(c("palbociclib","ribociclib","abemaciclib"), function(d)
    abs(net_growth_rate("rb_loss",d,"fulvestrant", prior_cdk46i="palbociclib")$r_eff - base)))
  add("RB-loss clone is CDK4/6i-agnostic", dmax < 1e-6,
      sprintf("max rate change by any CDK4/6i = %.2e", dmax))

  # 5. prior-drug rule nullifies palbociclib re-challenge (continue ~ ET)
  dmax2 <- max(sapply(STATE_NAMES, function(s) {
    et <- net_growth_rate(s,"none","fulvestrant", prior_cdk46i="palbociclib")$r_eff
    cp <- net_growth_rate(s,"palbociclib","fulvestrant", prior_cdk46i="palbociclib")$r_eff
    abs(cp - et) }))
  add("palbociclib re-challenge ~ ET (prior-drug rule)", dmax2 < 5e-4,
      sprintf("max |continue - ET| rate = %.2e", dmax2))

  # 6. monotonicity: median PFS decreases with aggressiveness
  g <- precompute_growth(prior_cdk46i="palbociclib")
  pop <- make_population(400, COMPOSITIONS$mixed_base, seed = 1)
  ms <- sapply(c(0.6,1.0,1.6,2.4), function(a)
    median(simulate_arm_pfs(pop, "ET_only", g, aggress_scale = a, dt = 2)$pfs_months))
  add("PFS monotonically decreasing in aggressiveness", all(diff(ms) <= 1e-9),
      sprintf("PFS = %s", paste(round(ms,2), collapse=" > ")))

  # 7. patient volume starts at 1 and progression threshold is RECIST +20% diam
  add("RECIST threshold = +20%% diameter (1.728x volume)",
      abs(RECIST$pd_vol_ratio - 1.2^3) < 1e-9,
      sprintf("vol ratio = %.3f", RECIST$pd_vol_ratio))

  do.call(rbind, checks)
}

# =============================================================================
# (B) LEAVE-ONE-TRIAL-OUT CROSS-VALIDATION
# Calibrate the shared composition on 3 trials; predict the held-out 4th.
# =============================================================================
loo_cv <- function(n_search = 140, n_eval = 450, n_final = 1500, dt = 4,
                   growth = NULL) {
  if (is.null(growth)) growth <- precompute_growth(prior_cdk46i = "palbociclib")
  tnames <- names(TRIALS)
  rows <- lapply(tnames, function(held) {
    train <- setdiff(tnames, held)
    cal <- calibrate_composition(n_search = n_search, n_eval = n_eval, dt = dt,
                                 seed = 202, which_trials = train, growth = growth)
    # predict held-out trial at high resolution
    pred <- evaluate_composition(cal$composition, n = n_final, dt = 1, seed = 2024,
                                 growth = growth, which_trials = held)
    data.frame(held_out = held,
               trained_on = paste(train, collapse = "+"),
               exp_pfs_obs = pred$exp_pfs_obs,
               exp_pfs_pred = pred$exp_pfs_pred,
               hr_obs = pred$hr_obs,
               hr_pred = pred$hr_pred,
               hr_lo = TRIALS[[held]]$hr_lo,
               hr_hi = TRIALS[[held]]$hr_hi,
               composition = paste(sprintf("%s=%.2f", names(cal$composition),
                                           cal$composition), collapse="; "),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$hr_in_CI <- (out$hr_pred >= out$hr_lo & out$hr_pred <= out$hr_hi)
  out
}

# =============================================================================
# (C) BOOTSTRAP PREDICTION INTERVALS
# Resample the virtual population to get CIs on predicted median PFS and HR.
# =============================================================================
bootstrap_predictions <- function(comp, growth, B = 600, n = 1500, dt = 1,
                                  seed = 2024) {
  set.seed(seed)
  pop <- make_population(n, comp, seed = seed)
  curve <- aggress_curve(pop, growth, "ET_only", dt = dt)
  rows <- lapply(names(TRIALS), function(tn) {
    tr <- TRIALS[[tn]]
    a  <- aggress_for_target(curve, tr$ctrl_pfs)
    ctrl <- simulate_arm_pfs(pop, tr$control_arm, growth, aggress_scale = a, dt = dt)
    expd <- simulate_arm_pfs(pop, tr$exp_arm,     growth, aggress_scale = a, dt = dt)
    med <- numeric(B); hr <- numeric(B)
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      med[b] <- median(expd$pfs_months[idx])
      hr[b]  <- cox_hr(ctrl[idx, ], expd[idx, ])
    }
    data.frame(trial = tn,
               exp_pfs_pred = median(expd$pfs_months),
               exp_pfs_lo = quantile(med, 0.025, names = FALSE),
               exp_pfs_hi = quantile(med, 0.975, names = FALSE),
               exp_pfs_obs = tr$exp_pfs,
               hr_pred = cox_hr(ctrl, expd),
               hr_pred_lo = quantile(hr, 0.025, names = FALSE),
               hr_pred_hi = quantile(hr, 0.975, names = FALSE),
               hr_obs = tr$hr, hr_obs_lo = tr$hr_lo, hr_obs_hi = tr$hr_hi,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# =============================================================================
# Main
# =============================================================================
if (.is_main("07_validation.R")) {
  t0 <- Sys.time()
  message(">> (A) Model verification ...")
  V <- verify_model()
  write.csv(V, file.path(DIR_TAB, "table_verification.csv"), row.names = FALSE)
  print(V, row.names = FALSE, right = FALSE)
  cat(sprintf("\nAll verification checks passed: %s\n", all(V$pass)))

  RES <- readRDS(file.path(DIR_RES, "trial_results.rds"))
  comp <- RES$composition; growth <- RES$growth

  message("\n>> (C) Bootstrap prediction intervals ...")
  BP <- bootstrap_predictions(comp, growth, B = 600, n = 1500, dt = 1)
  write.csv(BP, file.path(DIR_TAB, "table_prediction_intervals.csv"), row.names = FALSE)
  print(BP[, c("trial","exp_pfs_obs","exp_pfs_pred","exp_pfs_lo","exp_pfs_hi",
               "hr_obs","hr_pred","hr_pred_lo","hr_pred_hi")], row.names = FALSE, digits = 3)

  message("\n>> (B) Leave-one-trial-out cross-validation ...")
  LOO <- loo_cv(n_search = 140, n_eval = 450, n_final = 1500, dt = 4, growth = growth)
  write.csv(LOO, file.path(DIR_TAB, "table_loo_cv.csv"), row.names = FALSE)
  print(LOO[, c("held_out","trained_on","exp_pfs_obs","exp_pfs_pred",
                "hr_obs","hr_pred","hr_lo","hr_hi","hr_in_CI")],
        row.names = FALSE, digits = 3)

  # (D) Validation summary
  rmse_pfs <- sqrt(mean((BP$exp_pfs_pred - BP$exp_pfs_obs)^2))
  rmse_loghr <- sqrt(mean((log(BP$hr_pred) - log(BP$hr_obs))^2))
  in_ci <- sum(BP$hr_pred >= BP$hr_obs_lo & BP$hr_pred <= BP$hr_obs_hi)
  sign_ok <- sum(sign(BP$hr_pred - 1) == sign(BP$hr_obs - 1))
  loo_in_ci <- sum(LOO$hr_in_CI)
  summary_tab <- data.frame(
    criterion = c("Experimental-arm median PFS RMSE (months)",
                  "log-HR RMSE",
                  "Predicted HR within observed CI (of 4)",
                  "Correct HR direction vs 1.0 (of 4)",
                  "LOO-CV: held-out HR within observed CI (of 4)"),
    value = c(sprintf("%.2f", rmse_pfs),
              sprintf("%.3f", rmse_loghr),
              sprintf("%d / 4", in_ci),
              sprintf("%d / 4", sign_ok),
              sprintf("%d / 4", loo_in_ci)))
  write.csv(summary_tab, file.path(DIR_TAB, "table_validation_summary.csv"), row.names = FALSE)
  cat("\n=== Validation summary ===\n"); print(summary_tab, row.names = FALSE, right = FALSE)

  saveRDS(list(verification = V, bootstrap = BP, loo = LOO, summary = summary_tab),
          file.path(DIR_RES, "validation.rds"))
  cat(sprintf("\n07_validation.R complete; elapsed %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units="mins"))))
}
