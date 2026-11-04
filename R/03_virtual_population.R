# =============================================================================
# 03_virtual_population.R
# Virtual-patient generation and resistance-composition hypotheses.
#
# A second-line patient's residual tumour is modelled as a MIXTURE over the
# resistance-state archetypes (00_config.R). The patient-level volume is the
# weighted sum of subpopulation volumes, each growing/shrinking at the state-
# and-treatment-specific net rate r_eff. Because treatment inputs are constant
# within an arm, each subpopulation grows exponentially and the patient volume
# has the closed form
#
#     V_i(t) = sum_s  w_{i,s} * exp( r_{i,s} * t )
#
# with r_{i,s} = a_i * prolif_{s,arm} - death - kill_{s,arm}, where a_i is a
# per-patient "aggressiveness" multiplier on proliferation (lognormal) and
# w_{i,.} are Dirichlet-distributed subpopulation fractions around a composition
# hypothesis. Time-to-progression uses the RECIST +20%-diameter rule.
#
# COMPETING RESISTANCE-COMPOSITION HYPOTHESES
# -------------------------------------------
# The project does NOT pre-select one mechanism. It defines several population
# compositions and asks which can reproduce the full pattern of four trials:
#   drug_dependent : resistance mostly CDK6-remodelling (switch should work,
#                    continue-same should fail)
#   agnostic       : resistance mostly CDK2-bypass / RB-loss (nothing works)
#   er_driven      : resistance mostly ER-reactivation (ET-switch drives benefit)
#   mixed_base     : biology-informed blend (primary/ calibrated composition)
# =============================================================================

if (!exists("net_growth_rate")) {
  for (cand in c("01_model.R", "R/01_model.R")) if (file.exists(cand)) { source(cand); break }
}

# ---- Resistance-composition hypotheses (fractions over STATE_NAMES) ----------
# order: sensitive, cdk6_remodel, cdk2_bypass, rb_loss, er_reactivation
COMPOSITIONS <- list(
  mixed_base = c(sensitive = 0.10, cdk6_remodel = 0.40, cdk2_bypass = 0.20,
                 rb_loss = 0.08, er_reactivation = 0.22),
  drug_dependent = c(sensitive = 0.10, cdk6_remodel = 0.70, cdk2_bypass = 0.05,
                     rb_loss = 0.03, er_reactivation = 0.12),
  agnostic = c(sensitive = 0.08, cdk6_remodel = 0.07, cdk2_bypass = 0.50,
               rb_loss = 0.25, er_reactivation = 0.10),
  er_driven = c(sensitive = 0.10, cdk6_remodel = 0.12, cdk2_bypass = 0.10,
                rb_loss = 0.05, er_reactivation = 0.63)
)
stopifnot(all(sapply(COMPOSITIONS, function(x) abs(sum(x) - 1) < 1e-8)))
stopifnot(all(sapply(COMPOSITIONS, function(x) identical(names(x), STATE_NAMES))))

# ---- Treatment-arm definitions (CDK4/6i, ET) --------------------------------
# All post-progression arms use a switched/subsequent ET (fulvestrant here as
# the representative SERD, matching postMONARCH/PACE and MAINTAIN's switch).
ARMS <- list(
  ET_only        = list(cdk46i = "none",        et = "fulvestrant", label = "ET alone"),
  continue_palbo = list(cdk46i = "palbociclib", et = "fulvestrant", label = "continue/rechallenge palbociclib + ET"),
  switch_ribo    = list(cdk46i = "ribociclib",  et = "fulvestrant", label = "switch to ribociclib + ET"),
  switch_abema   = list(cdk46i = "abemaciclib", et = "fulvestrant", label = "switch to abemaciclib + ET")
)

# ---- Precompute growth components per (state, arm) ---------------------------
# Returns a list with matrices prolif[state,arm], death (scalar vector),
# kill[state,arm]. Net rate for patient with aggressiveness a:
#   r = a*prolif - death - kill
precompute_growth <- function(core = CORE, states = RESISTANCE_STATES,
                              arms = ARMS, state_names = STATE_NAMES,
                              prior_cdk46i = "palbociclib") {
  arm_names <- names(arms)
  prolif <- matrix(NA_real_, length(state_names), length(arm_names),
                   dimnames = list(state_names, arm_names))
  killm  <- prolif
  deathv <- setNames(rep(NA_real_, length(state_names)), state_names)
  for (s in state_names) for (a in arm_names) {
    g <- net_growth_rate(s, arms[[a]]$cdk46i, arms[[a]]$et, core, states,
                         prior_cdk46i = prior_cdk46i)
    prolif[s, a] <- g$prolif
    killm[s, a]  <- g$kill_drug
    deathv[s]    <- g$death
  }
  list(prolif = prolif, kill = killm, death = deathv, arms = arms,
       state_names = state_names, prior_cdk46i = prior_cdk46i)
}

# ---- Generate a virtual population ------------------------------------------
# Returns weights matrix [n x states] and aggressiveness vector.
make_population <- function(n, composition,
                            dirichlet_conc = 12,   # higher = less heterogeneity
                            aggress_logsd = 0.22,  # lognormal sd of aggressiveness
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  comp <- composition[STATE_NAMES]
  alpha <- comp * dirichlet_conc
  # Dirichlet via normalised Gammas
  W <- matrix(NA_real_, n, length(STATE_NAMES), dimnames = list(NULL, STATE_NAMES))
  for (j in seq_along(STATE_NAMES)) W[, j] <- rgamma(n, shape = max(alpha[j], 1e-3), rate = 1)
  W <- W / rowSums(W)
  aggress <- rlnorm(n, meanlog = -0.5 * aggress_logsd^2, sdlog = aggress_logsd) # mean 1
  list(weights = W, aggress = aggress, composition = comp)
}

# ---- Patient time-to-progression under one arm (vectorised over patients) -----
# aggress_scale multiplies each patient's proliferative aggressiveness, used to
# represent trial-level population differences (tumour burden / pace of disease).
simulate_arm_pfs <- function(pop, arm_name, growth,
                             aggress_scale = 1,
                             horizon = RECIST$horizon_days,
                             dt = RECIST$dt,
                             pd_vol_ratio = RECIST$pd_vol_ratio) {
  W      <- pop$weights                     # n x states
  agg    <- pop$aggress * aggress_scale     # length n
  prolif <- growth$prolif[colnames(W), arm_name]  # states
  kill   <- growth$kill[colnames(W),   arm_name]
  death  <- growth$death[colnames(W)]
  tgrid  <- seq(0, horizon, by = dt)
  n <- nrow(W); nt <- length(tgrid)

  # per-patient, per-state net rate: R[i,s] = agg_i * prolif_s - death_s - kill_s
  R <- outer(agg, prolif) -
       matrix(death + kill, n, length(prolif), byrow = TRUE)   # n x states

  Vmat <- matrix(0, n, nt)
  for (k in seq_len(nt)) {
    Vmat[, k] <- rowSums(W * exp(R * tgrid[k]))
  }
  # nadir (running min) and progression threshold per patient
  nad <- Vmat[, 1]
  ttp <- rep(horizon, n); prog <- rep(FALSE, n)
  for (k in seq_len(nt)) {
    nad <- pmin(nad, Vmat[, k])
    hit <- (!prog) & (Vmat[, k] > nad * pd_vol_ratio)
    if (any(hit)) { ttp[hit] <- tgrid[k]; prog[hit] <- TRUE }
  }
  data.frame(patient = seq_len(n), ttp_days = ttp, progressed = prog,
             pfs_months = ttp / 30.4375)
}

message("03_virtual_population.R loaded: ", length(COMPOSITIONS),
        " composition hypotheses, ", length(ARMS), " arms, ",
        "make_population(), precompute_growth(), simulate_arm_pfs()")
