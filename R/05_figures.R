# =============================================================================
# 05_figures.R
# Figures and summary tables for the QSP CDK4/6i resistance & rechallenge study.
# Reads results/trial_results.rds (produced by 04_trials.R).
# =============================================================================

if (!exists("simulate_arm_pfs")) {
  for (f in c("00_config.R","01_model.R","03_virtual_population.R","04_trials.R"))
    for (cand in c(f, file.path("R", f))) if (file.exists(cand)) { source(cand); break }
}

RES <- readRDS(file.path(DIR_RES, "trial_results.rds"))
comp   <- RES$composition
growth <- RES$growth
final  <- RES$final
fsr    <- RES$fsr
hyp    <- RES$hyp
TR     <- RES$TRIALS

theme_qsp <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(size = 9.3, colour = "grey30"),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

trial_levels <- c("MAINTAIN","PACE","postMONARCH","PALMIRA")
arm_cols <- c(control = "#4575b4", experimental = "#d73027")

# -----------------------------------------------------------------------------
# Helper: re-simulate patient-level data for a trial at calibrated composition
# -----------------------------------------------------------------------------
simulate_trial_patients <- function(trial_name, n = 1500, dt = 1, seed = 2024) {
  tr <- TR[[trial_name]]
  pop <- make_population(n, comp, seed = seed)
  curve <- aggress_curve(pop, growth, "ET_only", dt = dt)
  a <- aggress_for_target(curve, tr$ctrl_pfs)
  ctrl <- simulate_arm_pfs(pop, tr$control_arm, growth, aggress_scale = a, dt = dt)
  expd <- simulate_arm_pfs(pop, tr$exp_arm,     growth, aggress_scale = a, dt = dt)
  ctrl$arm <- "control"; expd$arm <- "experimental"
  list(data = rbind(ctrl, expd), aggress = a, tr = tr)
}

# -----------------------------------------------------------------------------
# FIG 1: Mechanism map - net growth rate heatmap (state x arm)
# -----------------------------------------------------------------------------
rmat <- growth$prolif - matrix(growth$death, nrow(growth$prolif), ncol(growth$prolif)) - growth$kill
rdf <- as.data.frame(as.table(rmat))
names(rdf) <- c("state","arm","r_eff")
rdf$state <- factor(rdf$state, levels = rev(STATE_NAMES),
                    labels = rev(c("Sensitive /\nER-dependent","CDK6 remodelling\n(drug-dependent)",
                               "CDK2 / cyclin E\nbypass (agnostic)","RB loss\n(agnostic)",
                               "ER reactivation\n(ESR1)")))
rdf$arm <- factor(rdf$arm, levels = c("ET_only","continue_palbo","switch_ribo","switch_abema"),
                  labels = c("ET alone","continue/\nrechallenge\npalbociclib",
                             "switch to\nribociclib","switch to\nabemaciclib"))
f1 <- ggplot(rdf, aes(arm, state, fill = r_eff)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%+.4f", r_eff)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                       midpoint = 0, name = "net growth\nrate r (per day)") +
  labs(title = "Mechanistic basis of the rechallenge dilemma",
       subtitle = "Net tumour-subpopulation growth rate by resistance state and treatment. Red = grows, blue = controlled.\nAgnostic clones (CDK2-bypass, RB-loss) are unaffected by any CDK4/6i; CDK6-remodelling responds only to a potent switch.",
       x = NULL, y = NULL) +
  theme_qsp + theme(legend.position = "right")
ggsave(file.path(DIR_FIG, "fig01_mechanism_heatmap.png"), f1, width = 9, height = 5, dpi = 150)

# -----------------------------------------------------------------------------
# FIG 3: Calibration scatter - observed vs predicted median PFS
# -----------------------------------------------------------------------------
pfs_long <- rbind(
  data.frame(trial = final$trial, arm = "control",
             obs = final$ctrl_pfs_obs, pred = final$ctrl_pfs_pred),
  data.frame(trial = final$trial, arm = "experimental",
             obs = final$exp_pfs_obs, pred = final$exp_pfs_pred))
f3 <- ggplot(pfs_long, aes(obs, pred, colour = arm, shape = trial)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3.4, stroke = 1) +
  scale_colour_manual(values = arm_cols, name = NULL) +
  scale_shape_manual(values = c(MAINTAIN=16,PACE=17,postMONARCH=15,PALMIRA=18), name = NULL) +
  coord_equal(xlim = c(2, 7), ylim = c(2, 7)) +
  labs(title = "Calibrated model reproduces observed median PFS",
       subtitle = "Control arms fit via per-trial aggressiveness; experimental arms are predictions. Dashed = identity.",
       x = "Observed median PFS (months)", y = "Model-predicted median PFS (months)") +
  theme_qsp
ggsave(file.path(DIR_FIG, "fig03_calibration_pfs.png"), f3, width = 6.6, height = 5.4, dpi = 150)

# -----------------------------------------------------------------------------
# FIG 4: Forest plot - observed HR (with CI) vs predicted HR
# -----------------------------------------------------------------------------
hrdf <- data.frame(
  trial = factor(final$trial, levels = rev(trial_levels)),
  hr_obs = final$hr_obs,
  hr_lo = sapply(final$trial, function(t) TR[[t]]$hr_lo),
  hr_hi = sapply(final$trial, function(t) TR[[t]]$hr_hi),
  hr_pred = final$hr_pred,
  regimen = sapply(final$trial, function(t) TR[[t]]$exp_label))
f4 <- ggplot(hrdf, aes(y = trial)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = hr_lo, xmax = hr_hi), height = 0.18, colour = "#4575b4") +
  geom_point(aes(x = hr_obs), colour = "#4575b4", size = 3) +
  geom_point(aes(x = hr_pred), colour = "#d73027", shape = 18, size = 4) +
  geom_text(aes(x = hr_obs, label = sprintf("obs %.2f", hr_obs)),
            vjust = -1.1, size = 2.9, colour = "#4575b4") +
  geom_text(aes(x = hr_pred, label = sprintf("pred %.2f", hr_pred)),
            vjust = 1.9, size = 2.9, colour = "#d73027") +
  scale_x_log10(breaks = c(0.4,0.6,0.8,1.0,1.3)) +
  labs(title = "Progression-free survival hazard ratios: observed vs QSP-predicted",
       subtitle = "Blue = observed HR (90-95% CI); red diamond = model prediction. HR<1 favours the experimental arm.",
       x = "Hazard ratio (experimental vs control, log scale)", y = NULL) +
  theme_qsp
ggsave(file.path(DIR_FIG, "fig04_forest_hr.png"), f4, width = 8, height = 4.6, dpi = 150)

# -----------------------------------------------------------------------------
# FIG 5: Identifiability - experimental PFS vs f_sr per trial
# -----------------------------------------------------------------------------
fsr_grid <- fsr$fsr_grid
curve_list <- lapply(names(TR), function(tn) {
  tr <- TR[[tn]]; a <- fsr$table$aggress[fsr$table$trial == tn]
  expv <- sapply(fsr_grid, function(fs) {
    pop <- make_population(700, RES$comp_from_fsr(fs), seed = 12)
    median_pfs(simulate_arm_pfs(pop, tr$exp_arm, growth, aggress_scale = a, dt = 2))
  })
  data.frame(trial = tn, fsr = fsr_grid, exp_pfs = pmin(expv, 24),
             obs = tr$exp_pfs, same_agent = tr$exp_arm == "continue_palbo")
})
fsr_curves <- do.call(rbind, curve_list)
fsr_curves$trial <- factor(fsr_curves$trial, levels = trial_levels)
f5 <- ggplot(fsr_curves, aes(fsr, exp_pfs)) +
  geom_line(aes(colour = same_agent), linewidth = 1.1) +
  geom_hline(aes(yintercept = obs), linetype = "dotted", colour = "grey30") +
  facet_wrap(~trial, scales = "free_y") +
  scale_colour_manual(values = c("FALSE"="#1b7837","TRUE"="#762a83"),
                      labels = c("switch to new agent","same-agent continuation"),
                      name = NULL) +
  labs(title = "Identifiability: how experimental-arm PFS depends on the switch-responsive fraction",
       subtitle = "Switch trials (green) are highly sensitive to f_sr and constrain the resistance mechanism;\nsame-agent continuation trials (purple) are flat - they cannot. Dotted line = observed PFS.",
       x = expression("Switch-responsive fraction  "*f[sr]),
       y = "Predicted median PFS (months)") +
  theme_qsp
ggsave(file.path(DIR_FIG, "fig05_identifiability.png"), f5, width = 8.5, height = 6, dpi = 150)

# -----------------------------------------------------------------------------
# FIG 6: Competing-hypothesis HR predictions vs observed
# -----------------------------------------------------------------------------
hyp2 <- hyp
hyp2$trial <- factor(hyp2$trial, levels = trial_levels)
hyp2$hypothesis <- factor(hyp2$hypothesis,
  levels = c("drug_dependent","mixed_base","er_driven","agnostic"),
  labels = c("drug-dependent","mixed","ER-driven","agnostic"))
obs_hr_df <- data.frame(trial = factor(names(sapply(TR, `[[`, "hr")), levels = trial_levels),
                        hr = sapply(TR, `[[`, "hr"))
f6 <- ggplot(hyp2, aes(trial, hr_pred, fill = hypothesis)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_point(data = obs_hr_df, aes(trial, hr), inherit.aes = FALSE,
             shape = 23, fill = "black", size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_fill_brewer(palette = "Set2", name = "population\nhypothesis") +
  labs(title = "Competing resistance-composition hypotheses vs observed HRs",
       subtitle = "Black diamonds = observed HR. A predominantly drug-dependent or ER-driven population predicts near-total\nswitch benefit (HR<<observed); a predominantly agnostic population predicts the smallest benefit. No single fixed\ncomposition reproduces all four trials - the aggregate data do not uniquely identify the mechanism.",
       x = NULL, y = "Predicted HR (experimental vs control)") +
  theme_qsp
ggsave(file.path(DIR_FIG, "fig06_hypotheses.png"), f6, width = 8.6, height = 5.2, dpi = 150)

# -----------------------------------------------------------------------------
# FIG 7: Kaplan-Meier curves per trial (control vs experimental)
# -----------------------------------------------------------------------------
km_list <- lapply(trial_levels, function(tn) {
  sim <- simulate_trial_patients(tn, n = 1500, dt = 1, seed = 2024)
  d <- sim$data
  d$event <- as.integer(d$progressed)
  fit <- survfit(Surv(pmax(ttp_days,0.5)/30.4375, event) ~ arm, data = d)
  s <- summary(fit, times = seq(0, 20, by = 0.5))
  data.frame(trial = tn, time = s$time,
             surv = s$surv,
             arm = sub("arm=", "", as.character(s$strata)))
})
km <- do.call(rbind, km_list)
km$trial <- factor(km$trial, levels = trial_levels)
f7 <- ggplot(km, aes(time, surv, colour = arm)) +
  geom_step(linewidth = 0.9, direction = "hv") +
  facet_wrap(~trial) +
  scale_colour_manual(values = arm_cols, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
  labs(title = "Simulated progression-free survival curves by trial",
       subtitle = "Virtual-population Kaplan-Meier curves at the calibrated resistance composition and per-trial aggressiveness.",
       x = "Months", y = "Progression-free (%)") +
  theme_qsp
ggsave(file.path(DIR_FIG, "fig07_km_curves.png"), f7, width = 8.6, height = 6, dpi = 150)

# -----------------------------------------------------------------------------
# FIG 8: Representative signalling trajectory (ODE) - sensitive vs cdk6_remodel
#        vs cdk2_bypass under abemaciclib switch (pRB and tumour volume)
# -----------------------------------------------------------------------------
traj_states <- c("sensitive","cdk6_remodel","cdk2_bypass","rb_loss")
traj <- do.call(rbind, lapply(traj_states, function(s) {
  sim <- simulate_tumour(s, cdk46i = "abemaciclib", et = "fulvestrant",
                         horizon = 360, dt = 2)
  sim$state <- s; sim
}))
traj$state <- factor(traj$state, levels = traj_states,
                     labels = c("sensitive","CDK6 remodel","CDK2 bypass","RB loss"))
pv <- ggplot(traj, aes(time/30.4375, V, colour = state)) +
  geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(x = "Months on abemaciclib + fulvestrant", y = "Subpopulation volume (rel.)",
       title = "Tumour subpopulation volume") + theme_qsp
pr <- ggplot(traj, aes(time/30.4375, pRB, colour = state)) +
  geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(x = "Months on abemaciclib + fulvestrant", y = "pRB (fraction phosphorylated)",
       title = "RB phosphorylation") + theme_qsp
f8 <- (pr | pv) + patchwork::plot_annotation(
  title = "Mechanistic trajectories under an abemaciclib switch",
  subtitle = "Abemaciclib suppresses RB phosphorylation and controls the sensitive and CDK6-remodelling clones, but not the\nCDK2-bypass or RB-loss (agnostic) clones, which continue to proliferate.",
  theme = theme(plot.title = element_text(face="bold", size=11.5),
                plot.subtitle = element_text(size=9.3, colour="grey30")))
ggsave(file.path(DIR_FIG, "fig08_trajectories.png"), f8, width = 9.5, height = 4.6, dpi = 150)

# -----------------------------------------------------------------------------
# FIG 9: Calibrated resistance composition
# -----------------------------------------------------------------------------
cdf <- data.frame(state = factor(names(comp), levels = STATE_NAMES,
                   labels = c("Sensitive /\nER-dependent","CDK6 remodel","CDK2 bypass",
                              "RB loss","ER reactivation")),
                  frac = as.numeric(comp))
f9 <- ggplot(cdf, aes(reorder(state, frac), frac, fill = state)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = percent(frac, accuracy = 1)), hjust = -0.15, size = 3.3) +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(c(0,0.12))) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Best shared resistance composition (calibrated to 4 trials)",
       subtitle = "Latin-hypercube calibration across the 4 randomized trials simultaneously.",
       x = NULL, y = "Fraction of residual tumour") +
  theme_qsp
ggsave(file.path(DIR_FIG, "fig09_composition.png"), f9, width = 7, height = 4.2, dpi = 150)

# -----------------------------------------------------------------------------
# Summary tables
# -----------------------------------------------------------------------------
tab_main <- data.frame(
  Trial = final$trial,
  `Experimental regimen` = sapply(final$trial, function(t) TR[[t]]$exp_label),
  `Control PFS obs` = round(final$ctrl_pfs_obs,2),
  `Control PFS pred` = round(final$ctrl_pfs_pred,2),
  `Exp PFS obs` = round(final$exp_pfs_obs,2),
  `Exp PFS pred` = round(final$exp_pfs_pred,2),
  `HR obs` = round(final$hr_obs,2),
  `HR pred` = round(final$hr_pred,2),
  check.names = FALSE)
write.csv(tab_main, file.path(DIR_TAB, "table_main_results.csv"), row.names = FALSE)

cat("Figures written to results/figures/ :\n")
print(list.files(DIR_FIG, pattern = "\\.png$"))
cat("\nTables written to results/tables/ :\n")
print(list.files(DIR_TAB, pattern = "\\.csv$"))
cat("\n05_figures.R complete.\n")
