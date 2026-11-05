# =============================================================================
# 99_run_all.R
# Master script: runs the entire QSP CDK4/6i resistance & rechallenge pipeline
# end to end and writes all figures and tables to results/.
#
# Usage:  Rscript R/99_run_all.R        (from the project root)
#     or  cd R && Rscript 99_run_all.R
# =============================================================================

t_start <- Sys.time()
here <- if (file.exists("00_config.R")) "." else "R"

message("\n========================================================")
message(" QSP of CDK4/6i resistance & rechallenge - full pipeline")
message("========================================================\n")

source(file.path(here, "00_config.R"))
source(file.path(here, "01_model.R"))
source(file.path(here, "02_calibration.R"))
source(file.path(here, "03_virtual_population.R"))
source(file.path(here, "04_trials.R"))

# ---- 1. Preclinical PD calibration (T47D) -----------------------------------
message("\n[1/5] Preclinical PD calibration to T47D switching data ...")
cal <- run_calibration()
write.csv(cal$calib, file.path(DIR_TAB, "table_calibration_T47D.csv"), row.names = FALSE)
{
  dr <- rbind(transform(cal$dr_ptxr, model_cell = "PTxR -> abemaciclib (drug-dependent)"),
              transform(cal$dr_atxr, model_cell = "ATxR -> palbociclib (CDK4/6-agnostic)"))
  p <- ggplot(dr, aes(theta, pRB, colour = model_cell)) +
    geom_line(linewidth = 1.1) +
    geom_hline(yintercept = 0.30, linetype = "dashed", colour = "grey40") +
    scale_x_log10() +
    scale_colour_manual(values = c("PTxR -> abemaciclib (drug-dependent)" = "#1b7837",
                                   "ATxR -> palbociclib (CDK4/6-agnostic)" = "#b2182b")) +
    labs(title = "Model pRB dose-response reproduces the T47D switching differential",
         x = expression("CDK4/6 inhibitor exposure (potency units "*theta*", log scale)"),
         y = "Quasi-steady-state pRB", colour = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(DIR_FIG, "fig02_calibration_pRB_doseresponse.png"), p,
         width = 8, height = 5.2, dpi = 150)
}

# ---- 2. Trial calibration, f_sr identifiability, hypothesis tests -----------
message("\n[2/5] Clinical calibration + identifiability + hypothesis tests ...")
cal2 <- calibrate_composition(n_search = 160, n_eval = 450, dt = 4, seed = 101)
comp <- cal2$composition
final <- evaluate_composition(comp, n = 1500, dt = 1, seed = 2024, growth = cal2$growth)
write.csv(transform(final, composition = paste(names(comp), round(comp,3), sep="=", collapse="; ")),
          file.path(DIR_TAB, "table_trial_fits.csv"), row.names = FALSE)
fsr <- analyse_fsr(cal2$growth, n = 900, dt = 1, seed = 11)
write.csv(fsr$table, file.path(DIR_TAB, "table_fsr_identifiability.csv"), row.names = FALSE)
hyp <- test_hypotheses(cal2$growth, n = 900, dt = 1, seed = 7)
write.csv(hyp, file.path(DIR_TAB, "table_hypothesis_tests.csv"), row.names = FALSE)
saveRDS(list(composition = comp, final = final, hyp = hyp, fsr = fsr,
             growth = cal2$growth, search_objs = cal2$objs,
             comp_from_fsr = comp_from_fsr, SR_PROFILE = SR_PROFILE, AG_PROFILE = AG_PROFILE,
             TRIALS = TRIALS, ARMS = ARMS, COMPOSITIONS = COMPOSITIONS),
        file.path(DIR_RES, "trial_results.rds"))

# ---- 3. Figures & summary tables --------------------------------------------
message("\n[3/5] Rendering figures and summary tables ...")
source(file.path(here, "05_figures.R"))

# ---- 4. Global sensitivity analysis -----------------------------------------
message("\n[4/5] Global sensitivity analysis (LHS + PRCC) ...")
source(file.path(here, "06_sensitivity.R"))
S <- run_sensitivity(N = 220, n_pat = 350, dt = 3)
Pdf <- as.data.frame(S$prcc); Pdf$parameter <- rownames(Pdf)
write.csv(Pdf, file.path(DIR_TAB, "table_sensitivity_prcc.csv"), row.names = FALSE)
saveRDS(S, file.path(DIR_RES, "sensitivity.rds"))
{
  pretty <- c(g_vol="tumour growth rate", death_bas="baseline cell loss",
              kE_syn="cyclin E synthesis", cdk2_cap="CDK2 capacity",
              kp46="RB phos. by CDK4/6", kp2="RB phos. by CDK2",
              er_dependence="ER-dependence of growth",
              palbo_penalty="palbociclib CDK6 penalty", ribo_penalty="ribociclib CDK6 penalty",
              abema_penalty="abemaciclib CDK6 penalty", prior_residual="prior-drug residual",
              f_sr="switch-responsive fraction")
  long <- do.call(rbind, lapply(c("HR_MAINTAIN","HR_postMONARCH"), function(o)
    data.frame(parameter = rownames(S$prcc), prcc = S$prcc[, o], output = o)))
  long$parameter <- factor(pretty[long$parameter],
                           levels = pretty[order(abs(S$prcc[,"HR_postMONARCH"]))])
  long$output <- factor(long$output, labels = c("MAINTAIN HR","postMONARCH HR"))
  g <- ggplot(long, aes(parameter, prcc, fill = prcc > 0)) +
    geom_col() + coord_flip() + facet_wrap(~output) +
    scale_fill_manual(values = c("TRUE"="#d73027","FALSE"="#4575b4"),
                      labels = c("TRUE"="increases HR (less benefit)","FALSE"="decreases HR (more benefit)"),
                      name = NULL) +
    labs(title = "Global sensitivity of predicted switch-trial HRs (PRCC)",
         x = NULL, y = "PRCC") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(DIR_FIG, "fig10_sensitivity_prcc.png"), g, width = 9, height = 5.2, dpi = 150)
}

# ---- 5. Verification, cross-validation & prediction intervals ---------------
message("\n[5/6] Model verification, leave-one-trial-out CV, bootstrap PIs ...")
source(file.path(here, "07_validation.R"))
Vchk <- verify_model()
write.csv(Vchk, file.path(DIR_TAB, "table_verification.csv"), row.names = FALSE)
message("   verification checks all pass: ", all(Vchk$pass))
BP <- bootstrap_predictions(comp, cal2$growth, B = 500, n = 1500, dt = 1)
write.csv(BP, file.path(DIR_TAB, "table_prediction_intervals.csv"), row.names = FALSE)
LOO <- loo_cv(n_search = 140, n_eval = 450, n_final = 1500, dt = 4, growth = cal2$growth)
write.csv(LOO, file.path(DIR_TAB, "table_loo_cv.csv"), row.names = FALSE)
rmse_pfs <- sqrt(mean((BP$exp_pfs_pred - BP$exp_pfs_obs)^2))
rmse_loghr <- sqrt(mean((log(BP$hr_pred) - log(BP$hr_obs))^2))
summary_tab <- data.frame(
  criterion = c("Experimental-arm median PFS RMSE (months)","log-HR RMSE",
                "Predicted HR within observed CI (of 4)","Correct HR direction (of 4)",
                "LOO-CV: held-out HR within observed CI (of 4)"),
  value = c(sprintf("%.2f", rmse_pfs), sprintf("%.3f", rmse_loghr),
            sprintf("%d / 4", sum(BP$hr_pred>=BP$hr_obs_lo & BP$hr_pred<=BP$hr_obs_hi)),
            sprintf("%d / 4", sum(sign(BP$hr_pred-1)==sign(BP$hr_obs-1))),
            sprintf("%d / 4", sum(LOO$hr_in_CI))))
write.csv(summary_tab, file.path(DIR_TAB, "table_validation_summary.csv"), row.names = FALSE)
saveRDS(list(verification=Vchk, bootstrap=BP, loo=LOO, summary=summary_tab),
        file.path(DIR_RES, "validation.rds"))

# Forest plot with prediction CIs (fig11) and LOO-CV (fig12)
{
  bp <- BP; bp$trial <- factor(bp$trial, levels = rev(c("MAINTAIN","PACE","postMONARCH","PALMIRA")))
  f11 <- ggplot(bp, aes(y = trial)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = hr_obs_lo, xmax = hr_obs_hi), height = 0.0, colour = "#4575b4", linewidth = 3, alpha = 0.25) +
    geom_point(aes(x = hr_obs), colour = "#4575b4", size = 3) +
    geom_errorbarh(aes(xmin = hr_pred_lo, xmax = hr_pred_hi), height = 0.16, colour = "#d73027") +
    geom_point(aes(x = hr_pred), colour = "#d73027", shape = 18, size = 3.6) +
    scale_x_log10(breaks = c(0.4,0.6,0.8,1.0,1.3)) +
    labs(title = "Observed vs QSP-predicted PFS hazard ratios (with prediction intervals)",
         subtitle = "Blue point/band = observed HR and CI; red diamond/bar = model prediction\nwith 95% bootstrap virtual-population (Monte-Carlo) interval.",
         x = "Hazard ratio (experimental vs control, log scale)", y = NULL) +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face="bold", size=11),
         plot.subtitle = element_text(size=9, colour="grey30"), panel.grid.minor = element_blank())
  ggsave(file.path(DIR_FIG, "fig11_forest_prediction_intervals.png"), f11, width = 8.4, height = 4.6, dpi = 150)

  lo <- LOO; lo$trial <- factor(lo$held_out, levels = c("MAINTAIN","PACE","postMONARCH","PALMIRA"))
  f12 <- ggplot(lo, aes(y = trial)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = hr_lo, xmax = hr_hi), height = 0.14, colour = "#4575b4") +
    geom_point(aes(x = hr_obs), colour = "#4575b4", size = 3) +
    geom_point(aes(x = hr_pred), colour = "#1b7837", shape = 17, size = 3.4) +
    scale_x_log10(breaks = c(0.3,0.4,0.6,0.8,1.0,1.3)) +
    labs(title = "Leave-one-trial-out cross-validation: held-out HR prediction",
         subtitle = "Each trial HR predicted from a composition fit ONLY on the other three.\nBlue = observed (CI); green triangle = out-of-sample prediction.",
         x = "Hazard ratio (log scale)", y = NULL) +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face="bold", size=11),
         plot.subtitle = element_text(size=8.8, colour="grey30"), panel.grid.minor = element_blank())
  ggsave(file.path(DIR_FIG, "fig12_loo_cv.png"), f12, width = 8.4, height = 4.4, dpi = 150)
}

# ---- 6. Done ----------------------------------------------------------------
message("\n[6/6] Pipeline complete.")
cat("\nFigures:\n"); print(list.files(DIR_FIG, pattern="\\.png$"))
cat("\nTables:\n");  print(list.files(DIR_TAB, pattern="\\.csv$"))
cat(sprintf("\nTotal elapsed: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
