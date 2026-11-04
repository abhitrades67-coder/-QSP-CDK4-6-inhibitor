# QSP of CDK4/6-Inhibitor Resistance & Rechallenge in HR+/HER2- Metastatic Breast Cancer

A quantitative systems pharmacology (QSP) model, implemented end-to-end in **R**,
that explains **why continuing a CDK4/6 inhibitor after resistance sometimes
works and sometimes fails** in HR+/HER2- metastatic breast cancer — and shows
which randomized trials can, and cannot, identify the operative mechanism.

> **Research question.** Which molecular resistance states retain sensitivity to
> CDK4/6 inhibition, how does this depend on the specific agent, and can a
> reduced-state QSP model explain the contrasting results of the randomized
> rechallenge trials (MAINTAIN, PACE, postMONARCH, PALMIRA)?

## TL;DR findings

- The model reconciles the four discordant trials with **two mechanistic rules**:
  (1) **agent-specific** potency against CDK6-remodelled resistance (abemaciclib,
  ribociclib re-engage it; palbociclib does not), and (2) **prior-drug
  cross-resistance** (re-challenging the failed agent adds almost nothing; a
  *switch* does).
- Predicted vs observed hazard ratios: MAINTAIN 0.62 (obs 0.57), PACE 0.99
  (1.11), postMONARCH 0.56 (0.73), PALMIRA 0.99 (0.84) — the **switch-vs-continue
  dichotomy is reproduced**.
- **Identifiability result:** switch trials strongly constrain the resistance
  mechanism; same-agent continuation trials are intrinsically uninformative — a
  concrete explanation for why the literature is ambiguous.
- Competing-hypothesis testing shows aggregate PFS does **not** uniquely identify
  the mechanism; the observed modest benefits imply **predominantly
  CDK4/6-agnostic** resistance after palbociclib.

See `report/final_report.md` for results and `proposal/research_proposal.md` for
the full proposal.

## Repository layout

```
.
├── proposal/research_proposal.md     # formal, publication-grade research proposal
├── report/final_report.md            # executed-study results and conclusions
├── data/
│   ├── clinical_trials.csv           # MAINTAIN/PACE/postMONARCH/PALMIRA endpoints (+DOIs)
│   ├── preclinical_T47D.csv          # T47D switching calibration targets (+DOI)
│   └── drug_pk.csv                   # agent exposure/potency anchors (+provenance)
├── R/
│   ├── 00_config.R                   # parameters, drug/ET/resistance-state libraries
│   ├── 01_model.R                    # QSP ODE + G1/S quasi-steady-state solver
│   ├── 02_calibration.R              # preclinical PD calibration (T47D)
│   ├── 03_virtual_population.R       # virtual patients + composition hypotheses
│   ├── 04_trials.R                   # trial simulation, calibration, hypothesis tests
│   ├── 05_figures.R                  # figures + summary tables
│   ├── 06_sensitivity.R             # global sensitivity (LHS + PRCC)
│   └── 99_run_all.R                  # master script (runs everything)
└── results/
    ├── figures/   fig01..fig10 (.png)
    └── tables/    table_*.csv
```

## Reproduce everything

Requires R (≥4.3) with `deSolve`, `survival`, `ggplot2`, `patchwork`, `scales`,
`lhs`, `data.table`, `dplyr`, `tidyr`.

```bash
Rscript R/99_run_all.R
```

Runtime ≈ 2 minutes. Outputs are written to `results/figures/` and
`results/tables/`. Individual stages can be run directly, e.g.:

```bash
cd R && Rscript 02_calibration.R     # preclinical calibration
cd R && Rscript 04_trials.R          # trial calibration + hypothesis tests
cd R && Rscript 06_sensitivity.R     # global sensitivity analysis
```

## Model in one paragraph

Each tumour subpopulation is governed by ER activity, cyclin D, the CDK4/6 and
CDK2 kinase inputs to RB, the E2F restriction-point switch, and tumour volume.
Endocrine therapy acts on ER (ligand depletion / transcriptional block /
degradation); CDK4/6 inhibitors reduce CDK4/6 activity with an **agent- and
state-specific** potency (a drug-specific CDK6 penalty). Five resistance-state
archetypes — sensitive, CDK6-remodelling (drug-dependent), CDK2/cyclin-E bypass
(agnostic), RB-loss (agnostic), ER-reactivation — are mixed to form virtual
patients. Time-to-progression uses the RECIST +20%-diameter rule. See the
proposal for equations and the report for results.

## Data & attribution

Clinical and preclinical anchors were retrieved from **PubMed**; every data file
and both write-ups carry the source DOIs. Please cite the original articles.

---

*Generated as a reproducible QSP study. Not medical advice; a research model of
published aggregate data, not individual patient data.*
