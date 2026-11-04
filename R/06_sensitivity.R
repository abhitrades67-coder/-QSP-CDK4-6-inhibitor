# =============================================================================
# 06_sensitivity.R
# Global sensitivity analysis (Latin-hypercube sampling + partial rank
# correlation, PRCC) of the predicted switch-trial hazard ratios and the
# switch-vs-continue PFS gap to the model's key parameters.
# =============================================================================

if (!exists("simulate_arm_pfs")) {
  for (f in c("00_config.R","01_model.R","03_virtual_population.R","04_trials.R"))
    for (cand in c(f, file.path("R", f))) if (file.exists(cand)) { source(cand); break }
}

set.seed(424242)

# ---- Parameters and ranges (uniform LHS) ------------------------------------
PARAMS <- list(
  g_vol         = c(0.0045, 0.0085),
  death_bas     = c(0.0004, 0.0016),
  kE_syn        = c(0.55,   0.90),
  cdk2_cap      = c(0.35,   0.70),
  kp46          = c(1.25,   1.80),
  kp2           = c(0.95,   1.45),
  er_dependence = c(0.40,   0.70),
  palbo_penalty = c(6.0,    12.0),
  ribo_penalty  = c(1.8,    3.6),
  abema_penalty = c(1.2,    2.2),
  prior_residual= c(0.05,   0.25),
  f_sr          = c(0.15,   0.65)
)

OUTPUTS <- c("HR_MAINTAIN", "HR_postMONARCH", "PFS_gap_switch_vs_continue")

# ---- One model evaluation under a parameter draw ----------------------------
eval_sample <- function(pv, n = 350, dt = 3) {
  core <- CORE
  core$g_vol        <- pv["g_vol"]
  core$death_bas    <- pv["death_bas"]
  core$kE_syn       <- pv["kE_syn"]
  core$cdk2_cap     <- pv["cdk2_cap"]
  core$kp46         <- pv["kp46"]
  core$kp2          <- pv["kp2"]
  core$er_dependence<- pv["er_dependence"]
  drugs <- DRUGS
  drugs$palbociclib$cdk6_penalty <- pv["palbo_penalty"]
  drugs$ribociclib$cdk6_penalty  <- pv["ribo_penalty"]
  drugs$abemaciclib$cdk6_penalty <- pv["abema_penalty"]

  # temporarily install globals the model functions read
  old_core <- CORE; old_drugs <- DRUGS; old_pr <- PRIOR_RESIDUAL
  assign("CORE", core, envir = .GlobalEnv)
  assign("DRUGS", drugs, envir = .GlobalEnv)
  assign("PRIOR_RESIDUAL", pv["prior_residual"], envir = .GlobalEnv)
  on.exit({ assign("CORE", old_core, envir = .GlobalEnv)
            assign("DRUGS", old_drugs, envir = .GlobalEnv)
            assign("PRIOR_RESIDUAL", old_pr, envir = .GlobalEnv) }, add = TRUE)

  growth <- precompute_growth(core = core, prior_cdk46i = "palbociclib")
  comp <- comp_from_fsr(pv["f_sr"])
  pop <- make_population(n, comp, seed = 99)
  curve <- aggress_curve(pop, growth, "ET_only", dt = dt)

  hr_for <- function(ctrl_pfs, exp_arm) {
    a <- aggress_for_target(curve, ctrl_pfs)
    ctrl <- simulate_arm_pfs(pop, "ET_only", growth, aggress_scale = a, dt = dt)
    expd <- simulate_arm_pfs(pop, exp_arm,  growth, aggress_scale = a, dt = dt)
    list(hr = cox_hr(ctrl, expd),
         gap = median(expd$pfs_months) - median(ctrl$pfs_months))
  }
  m <- hr_for(2.76, "switch_ribo")
  p <- hr_for(5.30, "switch_abema")
  c(HR_MAINTAIN = m$hr, HR_postMONARCH = p$hr,
    PFS_gap_switch_vs_continue = p$gap)
}

# ---- PRCC -------------------------------------------------------------------
prcc <- function(X, y) {
  Xr <- apply(X, 2, rank); yr <- rank(y)
  p <- ncol(Xr); out <- numeric(p); names(out) <- colnames(X)
  for (j in seq_len(p)) {
    others <- Xr[, -j, drop = FALSE]
    rx <- residuals(lm(Xr[, j] ~ others))
    ry <- residuals(lm(yr ~ others))
    out[j] <- cor(rx, ry)
  }
  out
}

run_sensitivity <- function(N = 220, n_pat = 350, dt = 3) {
  pnames <- names(PARAMS)
  L <- lhs::randomLHS(N, length(pnames))
  colnames(L) <- pnames
  X <- sapply(pnames, function(nm) PARAMS[[nm]][1] + L[, nm] * diff(PARAMS[[nm]]))
  Y <- matrix(NA_real_, N, length(OUTPUTS), dimnames = list(NULL, OUTPUTS))
  for (i in seq_len(N)) {
    Y[i, ] <- tryCatch(eval_sample(X[i, ], n = n_pat, dt = dt),
                       error = function(e) rep(NA_real_, length(OUTPUTS)))
    if (i %% 40 == 0) message("  sensitivity sample ", i, "/", N)
  }
  ok <- stats::complete.cases(Y)
  X <- X[ok, , drop = FALSE]; Y <- Y[ok, , drop = FALSE]
  P <- sapply(OUTPUTS, function(o) prcc(X, Y[, o]))
  list(X = X, Y = Y, prcc = P)
}

if (.is_main("06_sensitivity.R")) {
  message(">> Running global sensitivity analysis (LHS + PRCC)...")
  S <- run_sensitivity(N = 220, n_pat = 350, dt = 3)
  Pdf <- as.data.frame(S$prcc); Pdf$parameter <- rownames(Pdf)
  write.csv(Pdf, file.path(DIR_TAB, "table_sensitivity_prcc.csv"), row.names = FALSE)
  cat("\n=== PRCC (partial rank correlation) of outputs on parameters ===\n")
  print(round(S$prcc, 3))

  # Tornado plot for the two switch-trial HRs
  long <- do.call(rbind, lapply(c("HR_MAINTAIN","HR_postMONARCH"), function(o)
    data.frame(parameter = rownames(S$prcc), prcc = S$prcc[, o], output = o)))
  pretty <- c(g_vol="tumour growth rate", death_bas="baseline cell loss",
              kE_syn="cyclin E synthesis", cdk2_cap="CDK2 capacity",
              kp46="RB phos. by CDK4/6", kp2="RB phos. by CDK2",
              er_dependence="ER-dependence of growth",
              palbo_penalty="palbociclib CDK6 penalty",
              ribo_penalty="ribociclib CDK6 penalty",
              abema_penalty="abemaciclib CDK6 penalty",
              prior_residual="prior-drug residual activity",
              f_sr="switch-responsive fraction")
  long$parameter <- factor(pretty[long$parameter],
                           levels = pretty[order(abs(S$prcc[,"HR_postMONARCH"]))])
  long$output <- factor(long$output, labels = c("MAINTAIN HR","postMONARCH HR"))
  g <- ggplot(long, aes(parameter, prcc, fill = prcc > 0)) +
    geom_col() + coord_flip() +
    facet_wrap(~output) +
    scale_fill_manual(values = c("TRUE"="#d73027","FALSE"="#4575b4"),
                      labels = c("TRUE"="increases HR (less benefit)",
                                 "FALSE"="decreases HR (more benefit)"),
                      name = NULL) +
    labs(title = "Global sensitivity of predicted switch-trial HRs (PRCC)",
         subtitle = "Partial rank correlation of each parameter with the predicted hazard ratio across 220 LHS samples.",
         x = NULL, y = "PRCC") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face="bold", size=11.5),
          plot.subtitle = element_text(size=9.2, colour="grey30"),
          panel.grid.minor = element_blank())
  ggsave(file.path(DIR_FIG, "fig10_sensitivity_prcc.png"), g, width = 9, height = 5.2, dpi = 150)
  cat("\nSaved results/figures/fig10_sensitivity_prcc.png\n")
  cat("Saved results/tables/table_sensitivity_prcc.csv\n")
  saveRDS(S, file.path(DIR_RES, "sensitivity.rds"))
}
