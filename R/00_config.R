# =============================================================================
# 00_config.R
# Quantitative Systems Pharmacology of CDK4/6 inhibitor resistance & rechallenge
# in HR+/HER2- metastatic breast cancer.
#
# Global configuration: paths, physiological/PD parameters, drug library,
# resistance-state archetypes, and treatment-regimen definitions.
#
# All parameters are either (a) literature-informed, (b) calibrated to the
# preclinical T47D switching data (Zapatero-Solana 2025) and clinical PFS
# (MAINTAIN, PACE, postMONARCH, PALMIRA), or (c) reasonable defaults with a
# documented rationale. Provenance is in proposal/research_proposal.md sec 6.
# =============================================================================

suppressWarnings(suppressMessages({
  library(deSolve)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(survival)
  library(patchwork)
  library(lhs)
}))

set.seed(20251104)

# ---- Helper: is this file the one invoked via Rscript? -----------------------
.is_main <- function(name) {
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  length(a) == 1 && basename(a) == name
}

# Fraction of a CDK4/6 inhibitor's state-specific effect that REMAINS when it is
# the SAME agent the patient already progressed on (prior-drug cross-resistance).
# Re-challenge / continuation of the identical agent therefore adds little,
# whereas a SWITCH to a different agent retains its full state-specific activity.
PRIOR_RESIDUAL <- 0.12

# ---- Paths ------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a
# Robust project root detection (works via source() or Rscript)
.get_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) == 1 && nzchar(file_arg)) {
    return(normalizePath(file.path(dirname(file_arg), ".."), mustWork = FALSE))
  }
  if (dir.exists("R") && dir.exists("data")) return(normalizePath("."))
  if (basename(getwd()) == "R") return(normalizePath(".."))
  normalizePath(".")
}
ROOT       <- .get_root()
DIR_DATA   <- file.path(ROOT, "data")
DIR_R      <- file.path(ROOT, "R")
DIR_RES    <- file.path(ROOT, "results")
DIR_FIG    <- file.path(DIR_RES, "figures")
DIR_TAB    <- file.path(DIR_RES, "tables")
for (d in c(DIR_RES, DIR_FIG, DIR_TAB)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Core physiological / signalling parameters -----------------------------
# Time unit = days. Signalling species are dimensionless activity units.
# The G1/S restriction-point module (CDK4/6, CDK2, RB, E2F) is solved at
# quasi-steady state by damped fixed-point iteration each RHS call (see
# 01_model.R :: g1s_steady_state); only ER and tumour volume are integrated.
CORE <- list(
  # --- Estrogen receptor (ER) module ---
  kER_syn   = 1.0,    # ER activity production rate
  kER_deg   = 1.0,    # basal ER turnover
  estrogen  = 1.0,    # ligand drive (postmenopausal baseline, pre-AI)

  # --- Cyclin D module (integrates ER + growth-factor/MAPK bypass) ---
  kD_syn    = 1.2,    # cyclin D synthesis scale
  wER       = 0.70,   # weight of ER drive on cyclin D
  wGF       = 0.30,   # weight of growth-factor/MAPK drive on cyclin D
  gf_basal  = 0.15,   # basal growth-factor drive (RTK/MAPK); reference level
  kD_deg    = 1.0,    # cyclin D degradation

  # --- CDK4/6 complex ---
  cdk46_cap = 1.0,    # CDK4/6 catalytic capacity (scaled by state)
  p16       = 0.10,   # INK4/p16 inhibitory tone

  # --- CDK2 / cyclin E bypass ---
  # Baseline CDK2/cyclinE feedback is kept modest so that in non-bypass clones
  # the restriction-point switch tracks CDK4/6-inhibition level (a potent switch
  # can arrest a CDK6-remodelled clone). Genuine bypass states scale cdk2_cap /
  # kE_syn up so the feedback self-sustains and the clone is CDK4/6-agnostic.
  kE_syn    = 0.70,   # cyclin E produced per unit E2F (CycE = kE_syn * E2F)
  cdk2_cap  = 0.50,   # CDK2 catalytic capacity (scaled up in bypass states)
  cdk2_basal= 0.04,   # E2F-independent basal CDK2 activity

  # --- RB phosphorylation / E2F switch ---
  RB_total  = 1.0,    # total RB (→~0 in RB-loss state)
  kp46      = 1.50,   # RB phosphorylation rate by CDK4/6
  kp2       = 1.20,   # RB phosphorylation rate by CDK2
  kdp_RB    = 1.0,    # RB dephosphorylation (phosphatase)
  E2F_K     = 0.45,   # Hill K for E2F release vs pRB fraction
  E2F_n     = 2.2,    # Hill coefficient (graded restriction-point switch, so
                      # potent inhibition gives PARTIAL control, not full arrest,
                      # matching stable-disease-dominated clinical responses)

  # --- Proliferation / tumour volume kinetics (clinical readout) ---
  # CDK4/6 inhibitors + ET are predominantly CYTOSTATIC, so kill_scale is small
  # and controlled clones become ~stable (shallow nadir), not deeply regressing.
  # g_vol gives realistic HR+ MBC volume-doubling times (~4-6 months).
  g_vol     = 0.0060, # max tumour volume growth rate (per day) at full S-phase entry
  death_bas = 0.0008, # baseline tumour volume loss rate (per day); arrested
                      # clones become ~stable (shallow nadir), matching the low
                      # objective-response / high stable-disease rates clinically
  S_K       = 0.42,   # E2F activity for half-maximal S-phase entry
  S_n       = 2,      # steepness of S-phase entry vs E2F (graded)
  kill_scale= 0.0015, # scaling of (minimal) drug cytotoxic kill on cycling cells
  er_dependence = 0.55 # fraction of proliferative drive that is ER-dependent
                       # (so ET alone slows ER-dependent clones; ligand-
                       #  independent/ER-reactivated clones are unaffected)
)

# ---- Drug library -----------------------------------------------------------
# PD is parameterised in dimensionless "potency units":
#   theta = effective unbound exposure / cellular CDK4/6 IC50 (at CDK6 level 1).
# Fractional CDK4/6 inhibition = i_max * theta_eff^hpd / (1 + theta_eff^hpd),
# where theta_eff = theta / cdk6_factor and
#   cdk6_factor = 1 + (cdk6_penalty - 1) * (cdk6_level - 1).
# cdk6_penalty is DRUG-SPECIFIC: palbociclib loses most potency under CDK6
# over-expression, abemaciclib least (calibrated to the ~9-10x T47D
# pRB-depletion / IC50 differential, Zapatero-Solana 2025).
DRUGS <- list(
  none = list(theta = 0, hpd = 1, i_max = 0,
              breadth = 0, cytotox = 0, cdk6_penalty = 1),

  palbociclib = list(theta = 8.0,  hpd = 1.3, i_max = 0.97,
                     breadth = 0.10, cytotox = 0.12, cdk6_penalty = 9.0),
  ribociclib  = list(theta = 8.0,  hpd = 1.3, i_max = 0.97,
                     breadth = 0.25, cytotox = 0.28, cdk6_penalty = 2.6),
  abemaciclib = list(theta = 10.0, hpd = 1.3, i_max = 0.98,
                     breadth = 0.60, cytotox = 0.45, cdk6_penalty = 1.6)
)

# ---- Endocrine-therapy library ----------------------------------------------
# Effects on the ER module. fulvestrant degrades ER; AI depletes ligand;
# tamoxifen reduces transcriptional efficiency. 'potency' scales effect.
ENDOCRINE <- list(
  none        = list(er_degrade = 0.00, ligand_block = 0.00, txn_block = 0.00),
  AI          = list(er_degrade = 0.00, ligand_block = 0.75, txn_block = 0.00),
  tamoxifen   = list(er_degrade = 0.00, ligand_block = 0.00, txn_block = 0.60),
  fulvestrant = list(er_degrade = 0.80, ligand_block = 0.00, txn_block = 0.20)
)

# =============================================================================
# RESISTANCE-STATE ARCHETYPES
# Each archetype is a multiplicative/additive perturbation of CORE parameters.
# These are the competing MECHANISTIC HYPOTHESES the project tests. The clinical
# populations are modelled as MIXTURES over these states (see 04/05).
# =============================================================================
RESISTANCE_STATES <- list(

  # S0: Sensitive / residual-dependent. Still RB-E2F driven, modest adaptation.
  sensitive = list(
    label = "ER-dependent / sensitive",
    mods  = list(cdk46_cap = 1.0, cdk2_cap = 1.0, RB_total = 1.0,
                 gf_basal = 0.15, estrogen_independent = 0,
                 cdk6_level = 1.0, kE_syn = 1.0)
  ),

  # H1: CDK6 over-expression / CDK4-6-cyclinD complex remodelling.
  # DRUG-DEPENDENT resistance: potent/broad agents (abema) still engage; narrow
  # agent (palbo) does not. This is the state that makes a SWITCH work.
  cdk6_remodel = list(
    label = "CDK6 over-expression / complex remodelling",
    mods  = list(cdk46_cap = 1.8, cdk2_cap = 1.1, RB_total = 1.0,
                 gf_basal = 0.20, estrogen_independent = 0,
                 cdk6_level = 1.8, kE_syn = 1.1)
  ),

  # H2: CDK2 / cyclin E bypass. CDK4/6-AGNOSTIC: no CDK4/6i (same or switched)
  # restores control because RB is phosphorylated by CDK2 instead.
  cdk2_bypass = list(
    label = "CDK2 / cyclin E bypass",
    mods  = list(cdk46_cap = 1.0, cdk2_cap = 2.2, RB_total = 1.0,
                 gf_basal = 0.25, estrogen_independent = 0,
                 cdk6_level = 1.1, kE_syn = 1.8)
  ),

  # H3: RB loss. COMPLETE, agent-independent CDK4/6i resistance.
  rb_loss = list(
    label = "RB loss / functional inactivation",
    mods  = list(cdk46_cap = 1.0, cdk2_cap = 1.2, RB_total = 0.08,
                 gf_basal = 0.25, estrogen_independent = 0,
                 cdk6_level = 1.0, kE_syn = 1.3)
  ),

  # H4: ER reactivation (e.g. ESR1 mutation / ligand-independent ER).
  # ET-SWITCH-dependent: benefit comes from re-engaging ER with a different ET;
  # CDK4/6i still contributes because the RB axis is intact.
  er_reactivation = list(
    label = "ER reactivation (ESR1-mut / ligand-independent)",
    mods  = list(cdk46_cap = 1.2, cdk2_cap = 1.2, RB_total = 1.0,
                 gf_basal = 0.20, estrogen_independent = 1,
                 cdk6_level = 1.2, kE_syn = 1.2)
  )
)

STATE_NAMES <- names(RESISTANCE_STATES)

# ---- RECIST progression definition ------------------------------------------
# A patient tumour is a MIXTURE of subpopulations; total volume V(t) is the
# weighted sum of subpopulation volumes. RECIST diameter D = (V/V0)^(1/3).
# Progressive disease (PD) = 20% increase in DIAMETER over the on-treatment
# nadir (equivalently V > 1.728 * V_nadir), evaluated up to the horizon.
RECIST <- list(
  pd_fraction  = 0.20,   # +20% over nadir (diameter)
  pd_vol_ratio = 1.20^3, # equivalent volume ratio = 1.728
  horizon_days = 720,    # 24-month simulation horizon
  dt = 1                 # daily output grid
)

message("00_config.R loaded: ", length(RESISTANCE_STATES), " resistance states, ",
        length(DRUGS)-1, " CDK4/6 inhibitors, ROOT = ", ROOT)
