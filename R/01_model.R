# =============================================================================
# 01_model.R
# The reduced-resistance-state QSP model.
#
# Signalling network (quasi-mechanistic, reduced):
#
#   Estrogen --|ET|-->  ER  --+
#                             |--> CyclinD --> CDK4/6(active) --|CDK4/6i|
#   GF/MAPK bypass -----------+                                   |
#                                                                 v
#                                             pRB <-- (CDK4/6 + CDK2) phosphorylation
#                                               |
#                                               v  (restriction-point switch, Hill)
#                                             E2F --> CyclinE --> CDK2(active)  [bypass feedback]
#                                               |
#                                               v
#                                     S-phase entry --> tumour volume growth
#
# Integrated states y = (ER, V):
#   ER : ER transcriptional activity (ODE)
#   V  : subpopulation tumour volume (ODE; baseline 1 within a patient weighting)
#
# The G1/S module (CycD, CDK4/6, CDK2, CycE, pRB, E2F) is resolved to quasi-
# steady state each RHS evaluation by a damped fixed-point iteration starting
# from a neutral seed (E2F = 0.5). Because the seed is fixed, the converged
# branch is a deterministic function of the parameters (drug, resistance state),
# which is exactly the behaviour we want: a sensitive tumour under an effective
# CDK4/6i collapses to the arrested branch, whereas a CDK2-bypass tumour remains
# on the proliferative branch under the same drug.
# =============================================================================

if (!exists("CORE")) {
  for (cand in c("00_config.R", "R/00_config.R", file.path(dirname(getwd()), "R", "00_config.R")))
    if (file.exists(cand)) { source(cand); break }
}

# ---- Helper: Hill function --------------------------------------------------
hill <- function(x, K, n) {
  x <- pmax(x, 0)
  xn <- x^n
  xn / (xn + K^n)
}

# ---- Build an effective-parameter set for a (state, drug, ET) context -------
build_context <- function(state_name,
                          cdk46i = "none",
                          et = "fulvestrant",
                          core = CORE,
                          drugs = DRUGS,
                          endocrine = ENDOCRINE,
                          states = RESISTANCE_STATES,
                          prior_cdk46i = "none") {
  st <- states[[state_name]]
  if (is.null(st)) stop("Unknown resistance state: ", state_name)
  m <- st$mods
  d <- drugs[[cdk46i]]
  e <- endocrine[[et]]
  if (is.null(d)) stop("Unknown drug: ", cdk46i)
  if (is.null(e)) stop("Unknown endocrine therapy: ", et)

  p <- core
  p$cdk46_cap <- core$cdk46_cap * m$cdk46_cap
  p$cdk2_cap  <- core$cdk2_cap  * m$cdk2_cap
  p$RB_total  <- core$RB_total  * m$RB_total
  p$gf_basal  <- m$gf_basal
  p$kE_syn    <- core$kE_syn    * m$kE_syn
  p$cdk6_level<- m$cdk6_level
  p$estrogen_independent <- m$estrogen_independent

  p$drug        <- d
  p$et          <- e
  p$cdk46i_name <- cdk46i
  p$et_name     <- et
  p$state_name  <- state_name
  p$prior_cdk46i<- prior_cdk46i
  # untreated ER reference (et = none): ligand = estrogen (>= 0.85 floor if
  # ligand-independent), txn_eff = 1, degrade = 0
  ligand0 <- p$estrogen
  if (p$estrogen_independent == 1) ligand0 <- max(ligand0, 0.85)
  p$er_ref <- (p$kER_syn * ligand0) / p$kER_deg
  p
}

# ---- Fractional CDK4/6 inhibition by the current drug ------------------------
# Drug-specific CDK6 penalty raises the apparent potency requirement when CDK6
# is over-expressed / the complex is remodelled. Returns inhibition in [0,1].
cdk46_inhibition <- function(p) {
  d <- p$drug
  if (d$theta <= 0) return(0)
  cdk6_factor <- 1 + (d$cdk6_penalty - 1) * max(p$cdk6_level - 1, 0)
  theta_eff <- d$theta / cdk6_factor
  inh <- d$i_max * (theta_eff^d$hpd) / (1 + theta_eff^d$hpd)
  # prior-drug cross-resistance: re-using the agent the patient already
  # progressed on retains only a small residual effect
  if (!is.null(p$prior_cdk46i) && p$prior_cdk46i != "none" &&
      identical(p$cdk46i_name, p$prior_cdk46i)) {
    inh <- inh * PRIOR_RESIDUAL
  }
  inh
}

# ---- ER steady-state activity under endocrine therapy -----------------------
er_steady <- function(p) {
  e <- p$et
  ligand <- p$estrogen * (1 - e$ligand_block)
  if (p$estrogen_independent == 1) ligand <- max(ligand, 0.85)  # ligand-independent ER
  txn_eff <- (1 - e$txn_block)
  # dER/dt = kER_syn*ligand*txn_eff - kER_deg*ER*(1+degrade) = 0
  (p$kER_syn * ligand * txn_eff) / (p$kER_deg * (1 + e$er_degrade))
}

# ---- Cyclin D from ER + growth-factor/MAPK bypass ---------------------------
cyclinD_level <- function(p, ER) {
  drive <- p$wER * ER + p$wGF * (p$gf_basal / 0.15)   # bypass normalised to ref
  (p$kD_syn * drive) / p$kD_deg
}

# ---- G1/S restriction-point module at quasi-steady state --------------------
# Solves the implicit loop E2F -> CycE -> CDK2 -> pRB -> E2F by damped
# fixed-point iteration. Returns a list of module activities.
g1s_steady_state <- function(p, ER, max_iter = 200, tol = 1e-9, damp = 0.5) {
  inh  <- cdk46_inhibition(p)
  CycD <- cyclinD_level(p, ER)
  a46  <- p$cdk46_cap * CycD * (1 - p$p16) * (1 - inh)
  a46  <- max(a46, 0)

  rb_lost <- (p$RB_total < 0.2)

  E2F <- 0.5
  for (i in seq_len(max_iter)) {
    CycE <- p$kE_syn * E2F
    a2   <- p$cdk2_cap * (CycE + p$cdk2_basal)
    Krb  <- p$kp46 * a46 + p$kp2 * a2
    if (rb_lost) {
      pRB <- 1.0                     # restriction point abolished
    } else {
      pRB <- Krb / (Krb + p$kdp_RB)  # Michaelian RB phosphorylation fraction
    }
    E2F_new <- hill(pRB, p$E2F_K, p$E2F_n)
    E2F_upd <- (1 - damp) * E2F + damp * E2F_new
    if (abs(E2F_upd - E2F) < tol) { E2F <- E2F_upd; break }
    E2F <- E2F_upd
  }
  CycE <- p$kE_syn * E2F
  a2   <- p$cdk2_cap * (CycE + p$cdk2_basal)
  Krb  <- p$kp46 * a46 + p$kp2 * a2
  pRB  <- if (rb_lost) 1.0 else Krb / (Krb + p$kdp_RB)
  list(inh = inh, CycD = CycD, a46 = a46, a2 = a2, CycE = CycE,
       pRB = pRB, E2F = E2F,
       s_entry = hill(E2F, p$S_K, p$S_n))
}

# ---- ER-dependence factor on proliferation ----------------------------------
# A fraction (er_dependence) of proliferative drive requires ER signalling, so
# endocrine therapy that lowers ER slows growth of ER-dependent clones. ER is
# measured relative to the state's untreated ER; ligand-independent / ER-
# reactivated clones keep ER near baseline under ET, so are unaffected.
er_proliferation_factor <- function(p, ER) {
  er_rel <- ER / max(p$er_ref, 1e-6)
  fac <- (1 - p$er_dependence) + p$er_dependence * er_rel
  max(min(fac, 1.2), 0.05)
}

# ---- Growth components for a context at a given ER ---------------------------
# Returns proliferation, baseline death, and drug cytotoxic kill, plus the
# underlying module activities. Cytotoxic kill is gated by the DRUG-INDUCED DROP
# in S-phase entry relative to the same context without the CDK4/6 inhibitor, so
# a CDK4/6-agnostic (e.g. CDK2-bypass, RB-loss) clone - whose cycling the drug
# does not suppress - is neither arrested nor killed by the inhibitor.
growth_components <- function(p, ER) {
  g1s <- g1s_steady_state(p, ER)
  erf <- er_proliferation_factor(p, ER)
  # reference S-phase entry with the SAME ER / ET but no CDK4/6 inhibitor
  if (p$drug$theta > 0) {
    p0 <- p; p0$drug <- DRUGS[["none"]]
    s_ref <- g1s_steady_state(p0, ER)$s_entry
  } else {
    s_ref <- g1s$s_entry
  }
  s_drop    <- max(s_ref - g1s$s_entry, 0)
  prolif    <- p$g_vol * g1s$s_entry * erf
  kill_drug <- p$kill_scale * p$drug$cytotox * s_drop
  list(prolif = prolif, death = p$death_bas, kill_drug = kill_drug,
       r_eff = prolif - p$death_bas - kill_drug,
       s_entry = g1s$s_entry, s_ref = s_ref, s_drop = s_drop,
       pRB = g1s$pRB, E2F = g1s$E2F, inh = g1s$inh,
       a46 = g1s$a46, a2 = g1s$a2, erf = erf)
}

# =============================================================================
# ODE right-hand side: integrate ER and tumour volume V
# =============================================================================
qsp_rhs <- function(t, y, p) {
  ER <- y[["ER"]]
  V  <- y[["V"]]

  e  <- p$et
  ligand <- p$estrogen * (1 - e$ligand_block)
  if (p$estrogen_independent == 1) ligand <- max(ligand, 0.85)
  txn_eff <- (1 - e$txn_block)
  dER <- p$kER_syn * ligand * txn_eff - p$kER_deg * ER * (1 + e$er_degrade)

  gc <- growth_components(p, ER)
  r_eff <- gc$r_eff
  dV <- r_eff * V

  list(c(ER = dER, V = dV),
       a46 = gc$a46, a2 = gc$a2, pRB = gc$pRB, E2F = gc$E2F,
       s_entry = gc$s_entry, r_eff = r_eff, inh = gc$inh)
}

# ---- Net volume growth rate at steady state (closed form convenience) --------
# For constant inputs the subpopulation grows exponentially at this rate.
net_growth_rate <- function(state_name, cdk46i = "none", et = "fulvestrant",
                            core = CORE, states = RESISTANCE_STATES,
                            prior_cdk46i = "none") {
  p   <- build_context(state_name, cdk46i, et, core, states = states,
                       prior_cdk46i = prior_cdk46i)
  ER  <- er_steady(p)
  gc  <- growth_components(p, ER)
  list(r_eff = gc$r_eff,
       prolif = gc$prolif, death = gc$death, kill_drug = gc$kill_drug,
       s_entry = gc$s_entry, s_drop = gc$s_drop, pRB = gc$pRB,
       E2F = gc$E2F, inh = gc$inh, a46 = gc$a46, a2 = gc$a2)
}

# ---- Initial ER (untreated proliferative context) ---------------------------
untreated_ER <- function(state_name, core = CORE, states = RESISTANCE_STATES) {
  p <- build_context(state_name, cdk46i = "none", et = "none", core = core, states = states)
  er_steady(p)
}

# ---- Simulate one subpopulation under a treatment ---------------------------
# Returns a data.frame with time, ER, V, and module activities.
simulate_tumour <- function(state_name,
                            cdk46i = "none",
                            et = "fulvestrant",
                            V0 = 1.0,
                            core = CORE,
                            states = RESISTANCE_STATES,
                            horizon = RECIST$horizon_days,
                            dt = RECIST$dt) {
  p   <- build_context(state_name, cdk46i, et, core, states = states)
  ER0 <- untreated_ER(state_name, core, states)  # start from untreated ER, then treat
  init <- c(ER = ER0, V = V0)
  times <- seq(0, horizon, by = dt)
  out <- tryCatch(
    ode(y = init, times = times, func = qsp_rhs, parms = p,
        method = "lsoda", rtol = 1e-7, atol = 1e-9),
    error = function(e) NULL)
  if (is.null(out)) return(NULL)
  as.data.frame(out)
}

# ---- Time-to-progression from a (patient) volume trajectory -----------------
time_to_progression <- function(time, V,
                                 pd_vol_ratio = RECIST$pd_vol_ratio,
                                 horizon = RECIST$horizon_days) {
  nadir  <- cummin(V)
  thresh <- nadir * pd_vol_ratio
  hit <- which(V > thresh)
  if (length(hit) == 0) return(list(ttp = horizon, progressed = FALSE))
  list(ttp = time[hit[1]], progressed = TRUE)
}

message("01_model.R loaded: QSP (ER+V ODE, algebraic G1/S fixed point), ",
        "simulate_tumour(), net_growth_rate(), time_to_progression()")
