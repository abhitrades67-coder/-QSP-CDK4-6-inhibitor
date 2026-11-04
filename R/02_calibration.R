# =============================================================================
# 02_calibration.R
# Pharmacodynamic calibration of the QSP model to the preclinical T47D
# CDK4/6i-switching data (Zapatero-Solana et al., Int J Mol Sci 2025;
# DOI 10.3390/ijms26062643), cross-checked against Algethami et al.,
# Sci Rep 2025 (DOI 10.1038/s41598-025-11052-4).
#
# Mapping of experimental resistant models to QSP resistance states
# -----------------------------------------------------------------
#   PTxR (palbociclib + tamoxifen resistant) -> cdk6_remodel
#       Drug-DEPENDENT resistance: resistant to palbociclib, RESENSITISED by a
#       switch to the more potent/broader abemaciclib (or ribociclib partially).
#       Matches: PTxR + abemaciclib+fulvestrant -> decreased %GEM+/Ki67/colony,
#       near-complete pRB depletion at 250 nM, Fu-corrected IC50 ~1.1 uM.
#
#   ATxR (abemaciclib + tamoxifen resistant)  -> cdk2_bypass (CDK4/6-agnostic)
#       Cross-resistant: a switch to palbociclib does NOT restore control
#       (pRB maintained by CDK2); needs non-clinical palbociclib (2.2 uM) to
#       deplete pRB; Fu-corrected IC50 >10.7 uM; higher CDK6 and CDK2.
#
# The key quantitative anchor is the ~9-10x differential in potency required to
# suppress pRB between the drug-dependent (PTxR) and agnostic (ATxR) settings,
# plus the DIRECTION of response to the therapeutic switch.
# =============================================================================

if (!exists("qsp_rhs")) {
  for (cand in c("01_model.R", "R/01_model.R")) if (file.exists(cand)) { source(cand); break }
}

# ---- Explicit calibration contexts for the two T47D resistant lines ----------
# Both are CDK6-driven (drug-dependent) resistance; ATxR has higher CDK6 (and so
# needs a higher opposing-drug dose to re-suppress pRB, matching the 2.2 uM vs
# 250 nM pRB-depletion data and the ~1.8x CDK6 enrichment reported in ATxR).
# (The CDK2-bypass and RB-loss ARCHETYPES remain separate patient-level
# hypotheses; they are not these two cell lines.)
CALIB_STATES <- list(
  PTxR = list(label = "PTxR (palbociclib+tam resistant)",
              mods = list(cdk46_cap = 1.3, cdk2_cap = 1.05, RB_total = 1.0,
                          gf_basal = 0.18, estrogen_independent = 0,
                          cdk6_level = 1.3, kE_syn = 1.05)),
  ATxR = list(label = "ATxR (abemaciclib+tam resistant)",
              mods = list(cdk46_cap = 1.9, cdk2_cap = 1.10, RB_total = 1.0,
                          gf_basal = 0.20, estrogen_independent = 0,
                          cdk6_level = 2.0, kE_syn = 1.10))
)

# ---- Dose-response of pRB vs CDK4/6i exposure (potency units theta) ----------
# For a given resistant state + switch drug, sweep the drug's theta and record
# the quasi-steady-state pRB. Returns the theta at which pRB falls below a
# depletion threshold (the model analogue of the experimental pRB-depletion dose).
prb_dose_response <- function(state_name, drug_name, et = "fulvestrant",
                              theta_grid = 10^seq(-1, 2.6, length.out = 200),
                              core = CORE, states = CALIB_STATES) {
  base_drug <- DRUGS[[drug_name]]
  ER <- untreated_ER(state_name, core, states)
  res <- data.frame(theta = theta_grid, pRB = NA_real_, s_entry = NA_real_)
  for (i in seq_along(theta_grid)) {
    d <- base_drug; d$theta <- theta_grid[i]
    drugs2 <- DRUGS; drugs2[[drug_name]] <- d
    p <- build_context(state_name, drug_name, et, core, drugs2, states = states)
    g <- g1s_steady_state(p, ER)
    res$pRB[i] <- g$pRB
    res$s_entry[i] <- g$s_entry
  }
  res$state <- state_name; res$drug <- drug_name
  res
}

theta_for_prb <- function(dr, target = 0.30) {
  below <- which(dr$pRB <= target)
  if (length(below) == 0) return(NA_real_)
  dr$theta[below[1]]
}

# ---- Build calibration comparison -------------------------------------------
run_calibration <- function() {
  core <- CORE

  # Dose-response: PTxR switched to abemaciclib vs ATxR switched to palbociclib
  # (opposing-drug switch, matching the experimental design).
  dr_ptxr <- prb_dose_response("PTxR", "abemaciclib", states = CALIB_STATES)
  dr_atxr <- prb_dose_response("ATxR", "palbociclib", states = CALIB_STATES)

  th_ptxr <- theta_for_prb(dr_ptxr, 0.30)
  th_atxr <- theta_for_prb(dr_atxr, 0.30)
  model_fold <- if (is.na(th_atxr)) Inf else th_atxr / th_ptxr
  obs_fold   <- 2200 / 250   # ~8.8x from pRB-depletion doses; IC50 gives >9.7x

  # Response direction to the therapeutic switch at CLINICAL exposure
  # (net growth rate sign).
  g_ptxr_abema <- net_growth_rate("PTxR", "abemaciclib", "fulvestrant", states = CALIB_STATES)$r_eff
  g_ptxr_palbo <- net_growth_rate("PTxR", "palbociclib", "fulvestrant", states = CALIB_STATES)$r_eff
  g_atxr_palbo <- net_growth_rate("ATxR", "palbociclib", "fulvestrant", states = CALIB_STATES)$r_eff
  g_atxr_abema <- net_growth_rate("ATxR", "abemaciclib", "fulvestrant", states = CALIB_STATES)$r_eff
  g_atxr_ribo  <- net_growth_rate("ATxR", "ribociclib",  "fulvestrant", states = CALIB_STATES)$r_eff

  resp <- function(r) if (r < 0) "controlled (responds)" else "grows (resistant)"

  calib <- data.frame(
    readout = c(
      "PTxR + abemaciclib switch (clinical exposure)",
      "PTxR + palbociclib (drug that drove resistance)",
      "ATxR + palbociclib switch (clinical exposure)",
      "ATxR + ribociclib switch (clinical exposure)",
      "Potency fold-differential, pRB suppression (ATxR-palbo / PTxR-abema)",
      "CDK6 level, ATxR vs PTxR (model contexts)",
      "CDK2 capacity, ATxR vs PTxR (model contexts)"
    ),
    model = c(
      resp(g_ptxr_abema),
      resp(g_ptxr_palbo),
      resp(g_atxr_palbo),
      resp(g_atxr_ribo),
      sprintf("%.1fx", model_fold),
      sprintf("%.2fx", CALIB_STATES$ATxR$mods$cdk6_level /
                        CALIB_STATES$PTxR$mods$cdk6_level),
      sprintf("%.2fx", CALIB_STATES$ATxR$mods$cdk2_cap /
                        CALIB_STATES$PTxR$mods$cdk2_cap)
    ),
    observed = c(
      "responds (decreased %GEM+/Ki67/colony, apoptosis)",
      "resistant (selected under palbociclib)",
      "cross-resistant (no decrease %GEM+/Ki67/colony)",
      "partial / non-significant trend",
      sprintf("~%.1fx (2200 vs 250 nM pRB; >9.7x by IC50)", obs_fold),
      "~1.8x (CDK6 enriched in ATxR)",
      "~1.5x (CDK2/cyclinD greater in ATxR)"
    ),
    source = "Zapatero-Solana 2025 (10.3390/ijms26062643)",
    stringsAsFactors = FALSE
  )

  list(calib = calib, dr_ptxr = dr_ptxr, dr_atxr = dr_atxr,
       th_ptxr = th_ptxr, th_atxr = th_atxr,
       model_fold = model_fold, obs_fold = obs_fold)
}

# ---- Execute & save ---------------------------------------------------------
if (.is_main("02_calibration.R")) {
  cal <- run_calibration()
  write.csv(cal$calib, file.path(DIR_TAB, "table_calibration_T47D.csv"), row.names = FALSE)
  cat("\n=== Preclinical PD calibration (T47D switching) ===\n")
  print(cal$calib, right = FALSE)
  cat(sprintf("\nModel potency fold-differential (ATxR palbo / PTxR abema) = %.1fx ; observed ~8.8-9.7x\n",
              cal$model_fold))

  # Dose-response figure
  dr <- rbind(
    transform(cal$dr_ptxr, model_cell = "PTxR -> abemaciclib (drug-dependent)"),
    transform(cal$dr_atxr, model_cell = "ATxR -> palbociclib (CDK4/6-agnostic)")
  )
  p <- ggplot(dr, aes(theta, pRB, colour = model_cell)) +
    geom_line(linewidth = 1.1) +
    geom_hline(yintercept = 0.30, linetype = "dashed", colour = "grey40") +
    annotate("text", x = 0.12, y = 0.33, label = "pRB-depletion threshold",
             hjust = 0, size = 3, colour = "grey30") +
    scale_x_log10() +
    scale_colour_manual(values = c("PTxR -> abemaciclib (drug-dependent)" = "#1b7837",
                                   "ATxR -> palbociclib (CDK4/6-agnostic)" = "#b2182b")) +
    labs(title = "Model pRB dose-response reproduces the T47D switching differential",
         subtitle = "Drug-dependent resistance (PTxR) is re-suppressed by a potent switch; agnostic resistance (ATxR) is not",
         x = expression("CDK4/6 inhibitor exposure (potency units, "*theta*", log scale)"),
         y = "Quasi-steady-state pRB (fraction phosphorylated)",
         colour = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 11))
  ggsave(file.path(DIR_FIG, "fig02_calibration_pRB_doseresponse.png"), p,
         width = 8, height = 5.2, dpi = 150)
  cat("Saved: results/figures/fig02_calibration_pRB_doseresponse.png\n")
  cat("Saved: results/tables/table_calibration_T47D.csv\n")
}
