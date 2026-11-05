# Quantitative Systems Pharmacology of CDK4/6-Inhibitor Resistance and Rechallenge: Integrating Cell-Cycle Signalling, Drug Exposure and Clinical Outcomes in Metastatic Breast Cancer

**A research proposal**

Indication: Hormone-receptor-positive (HR+), HER2-negative metastatic breast cancer (MBC)
Modality: Quantitative Systems Pharmacology (QSP) / translational mathematical oncology
Implementation: R (`deSolve`, `survival`, `ggplot2`, `lhs`)

---

## 1. Lay summary

Nearly all women with the most common form of advanced breast cancer are treated with a class of drugs called CDK4/6 inhibitors (palbociclib, ribociclib, abemaciclib) combined with hormone therapy. These drugs work well at first, but the cancer eventually learns to grow despite them. A pressing clinical question is what to do next: **stop the CDK4/6 inhibitor, keep the same one, or switch to a different one?** Randomized trials have given apparently conflicting answers — some show a real benefit from continuing CDK4/6 blockade, others show none. This project builds a mechanistic computer model of the cancer-cell "engine" that these drugs target, calibrates it to recent laboratory data and to the trials themselves, and uses it to explain **why the same strategy helps some patients and fails others** — and, importantly, to show which of these clinical questions the existing trials can and cannot actually answer.

---

## 2. Background and significance

### 2.1 The target and the drugs
In HR+/HER2- breast cancer, oestrogen-receptor (ER) signalling drives transcription of **cyclin D1**, which activates **CDK4/6**. Active CDK4/6 phosphorylates the retinoblastoma protein (**RB**), relieving its repression of **E2F** transcription factors and permitting passage through the G1/S restriction point into DNA synthesis. CDK4/6 inhibitors (CDK4/6i) hold RB hypophosphorylated and arrest cells in G1. Three agents are approved — palbociclib, ribociclib and abemaciclib — with differing potency, continuous-vs-intermittent schedules, off-target breadth (abemaciclib is the broadest, with appreciable CDK2/CDK9 activity and a cytotoxic component) and unbound exposure.

### 2.2 The clinical dilemma
Combined with endocrine therapy (ET), CDK4/6i are standard first-line therapy, but resistance is nearly universal. After progression, whether to **maintain CDK4/6 blockade** is unresolved, and the randomized evidence is genuinely discordant:

| Trial (randomized) | Post-progression strategy | Median PFS (exp vs control) | HR (95% CI) |
|---|---|---|---|
| **MAINTAIN** | switch ET **+ ribociclib** vs switch ET + placebo | **5.29 vs 2.76 mo** | 0.57 (0.39–0.85) |
| **PACE** | fulvestrant **+ palbociclib** (same agent) vs fulvestrant | 4.60 vs 4.80 mo | 1.11 (0.79–1.55) |
| **postMONARCH** | fulvestrant **+ abemaciclib** vs fulvestrant | **6.00 vs 5.30 mo** | 0.73 (0.57–0.95) |
| **PALMIRA** | **palbociclib** rechallenge + ET (same agent) vs ET | 4.90 vs 3.60 mo | 0.84 (0.66–1.07) |

Comparisons are within-study; populations and prior therapy differ across studies (sources: MAINTAIN [1], PACE [2], postMONARCH [3], PALMIRA [4]; pooled analysis [5]; review [6]). A striking pattern emerges: trials that **switched to a different, more potent/continuously dosed CDK4/6i** (MAINTAIN, postMONARCH) were positive, whereas trials that **continued or re-challenged the same agent** patients had already failed (PACE and PALMIRA, both palbociclib) were null or marginal.

### 2.3 Mechanistic substrate for a model
Resistance to CDK4/6i is heterogeneous and, crucially, partly **drug-dependent**. A 2025 study in HR+/HER2- T47D cells [7] showed that palbociclib-plus-tamoxifen–resistant cells (PTxR) **re-responded** to a switch to abemaciclib (near-complete pRB depletion at 250 nM, Fu-corrected IC50 ~1.1 µM, decreased %GEM+/Ki-67/colony formation and increased apoptosis), whereas abemaciclib-resistant cells (ATxR) were **cross-resistant** to palbociclib (pRB depletion only at a non-clinical 2.2 µM; IC50 >10.7 µM), with higher CDK6 and CDK2. An independent 2025 panel [8] generated ribociclib-, abemaciclib- and palbociclib-resistant MCF7/T47D lines, confirmed cross-resistance, and identified G2/M-bypass routes (PLK1, Aurora B). ctDNA from PACE [9] shows the clinical resistance landscape is a **mixture** — ESR1 (54%), TP53 (36%), PIK3CA (34%), RB1 (10%) — with RB1/cell-cycle and PI3K-pathway mutations marking shorter PFS. Established resistance mechanisms thus span (i) **CDK6 over-expression / CDK4-6-cyclin D complex remodelling** (drug-dependent; a potent switch can still engage), (ii) **CDK2 / cyclin E bypass** (CDK4/6-agnostic), (iii) **RB loss** (agnostic), and (iv) **ER reactivation / ESR1 mutation** (ET-switch-dependent).

### 2.4 Prior modelling and the gap
Mathematical cell-cycle models of CDK4/6i exist — e.g. a 2025 model of palbociclib alternating with fulvestrant + WEE1 inhibition that delays resistance [10] — but they address a *single* mechanism and a *single* drug combination, and are not connected to the discordant rechallenge trials. **No existing model integrates (a) the cyclin D–CDK4/6–RB–E2F axis, ER signalling and the CDK2/cyclin E bypass; (b) agent-specific drug exposure and potency; (c) a heterogeneous mixture of resistance states; and (d) the four randomized rechallenge trials as a joint calibration/validation target.** This proposal fills that gap.

### 2.5 The biological question
> **Which molecular resistance states retain sensitivity to CDK4/6 inhibition, how does this depend on the specific agent, and can a reduced-state QSP model explain the contrasting results of the randomized rechallenge trials — or show that the aggregate clinical data cannot uniquely identify the operative mechanism?**

---

## 3. Objectives, hypotheses and specific aims

**Central hypothesis.** The discordant rechallenge trials are reconciled by two factors acting on a *heterogeneous mixture* of resistance states: (1) **agent-specific** engagement of CDK6-remodelled (drug-dependent) resistance — a potent switch (abemaciclib, ribociclib) re-engages it, re-challenge of the already-failed agent does not; and (2) the **fraction of the residual tumour that is CDK4/6-agnostic** (CDK2-bypass, RB-loss), which no CDK4/6i can control.

**Competing hypotheses to be tested (not pre-selected).** Because aggregate PFS/HR may not identify the mechanism uniquely, we explicitly test several population-composition hypotheses:
- **H1 (drug-dependent-dominant):** residual tumour is predominantly CDK6-remodelling.
- **H2 (agnostic-dominant):** predominantly CDK2-bypass / RB-loss.
- **H3 (ER-driven-dominant):** predominantly ER-reactivation.
- **H4 (mixed):** biologically-informed blend.

- **Aim 1 — Build and verify a reduced-state QSP model** of the cyclin D–CDK4/6–RB–E2F switch with ER input and a CDK2/cyclin E bypass, parameterised by five resistance-state archetypes and three CDK4/6 inhibitors.
- **Aim 2 — Calibrate the pharmacodynamics** to the preclinical T47D switching data [7] (pRB dose-response, response direction, relative CDK6/CDK2 levels), reproducing the drug-dependent vs agnostic distinction.
- **Aim 3 — Simulate the four randomized trials** with a virtual population, calibrate the resistance composition, and **test the competing hypotheses** against the joint PFS/HR pattern.
- **Aim 4 — Quantify identifiability and sensitivity:** determine which trial designs constrain the mechanism and which do not, and rank the parameters driving predicted benefit.

---

## 4. Research design and methods

### 4.1 Model structure (Aim 1)
State variables (per tumour subpopulation): ER activity, cyclin D, cyclin E, phosphorylated-RB fraction (pRB), free E2F, and tumour volume *V*. The G1/S restriction-point module is solved to quasi-steady state each evaluation by damped fixed-point iteration (deterministic branch selection, avoiding integration artefacts of the underlying bistable switch); ER and *V* are integrated with `lsoda`. Governing relations (schematically):

- ER: `dER/dt = kER·ligand·txn − kER·ER·(1+degrade)`, with ET acting through ligand depletion (AI), transcriptional block (tamoxifen) or ER degradation (fulvestrant).
- Cyclin D: `CycD = kD·(wER·ER + wGF·GF)`, integrating ER and a growth-factor/MAPK bypass.
- CDK4/6 activity: `a46 = cap46·CycD·(1−p16)·(1−inh)`, where `inh` is agent- and state-specific (below).
- Restriction point: `pRB = Krb/(Krb+kdp)` with `Krb = kp46·a46 + kp2·a2`; `E2F = Hill(pRB)`; `CycE = kE·E2F`; `a2 = cap2·(CycE+basal)` (the CDK2 bypass feeds back).
- Growth: `dV/dt = [g·S(E2F)·f_ER(ER) − death − kill]·V`, with S-phase entry `S = Hill(E2F)` and a modest, arrest-gated cytotoxic term.

**Agent-specific potency.** Fractional CDK4/6 inhibition is `i_max·θ_eff^h/(1+θ_eff^h)` with `θ_eff = θ / [1 + (penalty−1)(CDK6_level−1)]`. The **CDK6 penalty is drug-specific** (palbociclib ≫ ribociclib > abemaciclib), encoding the observation that abemaciclib overcomes CDK6-driven resistance where palbociclib cannot.

**Prior-drug cross-resistance.** Because every enrolled patient progressed on palbociclib, re-use of the identical agent retains only a small residual activity (`PRIOR_RESIDUAL`), so *continuation/re-challenge* of palbociclib adds little while a *switch* retains full state-specific activity. This single, biologically-grounded rule distinguishes the "continue" from the "switch" trials.

**Resistance-state archetypes:** sensitive/ER-dependent; CDK6-remodelling (drug-dependent); CDK2/cyclin E bypass (agnostic); RB-loss (agnostic); ER-reactivation (ESR1/ligand-independent).

### 4.2 Preclinical calibration (Aim 2)
Map the two experimental resistant lines to model contexts (PTxR → milder, ATxR → higher CDK6 grade of drug-dependent resistance) and calibrate the drug-specific CDK6 penalties and potencies to reproduce: (i) the **direction** of the switch response (PTxR controlled by abemaciclib; ATxR cross-resistant to palbociclib, partial to ribociclib); (ii) the **~9-fold (lower-bound) potency differential** for pRB suppression; (iii) the relative CDK6/CDK2 levels. Calibration anchors are the pRB dose-response and net-growth-rate sign. Data: Zapatero-Solana 2025 [7], cross-checked against Algethami 2025 [8].

### 4.3 Clinical simulation and hypothesis testing (Aim 3)
- **Virtual population.** Each patient is a Dirichlet-distributed mixture over the five states plus a log-normal proliferative-aggressiveness multiplier; the tumour volume is the weighted sum of subpopulation volumes, and time-to-progression uses the RECIST +20%-diameter rule over a 24-month horizon.
- **Trial mapping.** MAINTAIN (switch ribociclib), PACE (continue palbociclib), postMONARCH (switch abemaciclib), PALMIRA (rechallenge palbociclib); all control arms are ET-only.
- **Calibration.** The **per-trial aggressiveness** is fit to each control arm's observed median PFS (so control arms are matched); the **experimental arms are predictions**. A shared resistance composition is calibrated across all four trials by Latin-hypercube search.
- **Competing hypotheses.** Each fixed composition (H1–H4) is scored by its ability to reproduce the four observed HRs jointly.
- **Outcome metrics.** Median PFS, Kaplan–Meier curves and Cox hazard ratios (vs observed with CIs).

### 4.4 Identifiability and global sensitivity (Aim 4)
- **Identifiability.** Re-parameterise the composition by a single interpretable **switch-responsive fraction f_sr** and map how each trial's experimental-arm PFS depends on f_sr. The *range* of this dependence quantifies how strongly a given trial design constrains the mechanism.
- **Global sensitivity.** Latin-hypercube sampling over all key parameters with **partial rank correlation (PRCC)** against the predicted switch-trial HRs and the switch-vs-continue PFS gap.

### 4.5 Drug-exposure (PK) inputs
Effective unbound exposures and cellular CDK4/6 IC50s are expressed as dimensionless potency ratios θ (literature-informed; see `data/drug_pk.csv`), with abemaciclib's active metabolites (M2/M20), continuous dosing and off-target breadth folded into its effective potency and cytotoxic capacity. These are treated as **calibratable PD anchors**, not asserted PK facts, with provenance documented.

### 4.6 Software and reproducibility
All steps run in R (`deSolve`, `survival`, `ggplot2`, `patchwork`, `lhs`, `data.table`, `dplyr`). A single master script (`R/99_run_all.R`) reproduces every figure and table from source; the model, calibration, trial engine, figures and sensitivity analysis are modularised (`R/00`–`R/06`). Random seeds are fixed.

---

## 5. Validation plan and success criteria
1. **Preclinical face validity:** the model reproduces the sign of the T47D switch responses and the order-of-magnitude pRB-potency differential. *(Success: all response directions correct; differential consistent with the experimental lower bound.)*
2. **Clinical reproduction:** predicted experimental-arm HRs lie within, or adjacent to, the observed confidence intervals for ≥3 of 4 trials, with the correct sign for all four. *(Success criterion met if the switch-vs-continue dichotomy is reproduced.)*
3. **Model verification:** automated checks that the closed-form population layer matches the full ODE, the G1/S fixed point is unique/seed-independent, invariants hold (agnostic-clone insensitivity, prior-drug rule), and PFS is monotone in aggressiveness.
4. **Out-of-sample validation:** **leave-one-trial-out cross-validation** — the resistance composition is calibrated on three trials and used to predict the held-out trial's HR, with bootstrap virtual-population prediction intervals. *(Success: ≥3 of 4 held-out HRs within the observed CI.)*
5. **Identifiability:** quantitative demonstration of which trials constrain the mechanism.
6. **Sensitivity:** the drivers of benefit are mechanistically interpretable and agent-specific (PRCC).

---

## 6. Expected outcomes and impact
We expect the model to show that (i) the switch-vs-continue dichotomy is explained by prior-drug cross-resistance plus agent-specific engagement of drug-dependent (CDK6-remodelling) clones; (ii) the magnitude of benefit is governed chiefly by the *switch-responsive fraction* of the residual tumour; and (iii) **same-agent continuation trials are intrinsically uninformative about the resistance mechanism**, whereas switch trials are highly informative — a direct, actionable statement about trial design and about the limits of what aggregate PFS can establish. The framework provides a reusable platform for prioritising biomarker-stratified rechallenge strategies (e.g. selecting patients with a high drug-dependent fraction for a potent switch) and for in-silico design of future trials.

---

## 7. Limitations and mitigations
- **Reduced model.** The G1/S network is deliberately coarse-grained; we mitigate by anchoring to direct pRB/proliferation measurements and by global sensitivity analysis.
- **Non-uniqueness.** Aggregate PFS cannot uniquely identify composition; this is treated as a *result* (competing-hypothesis testing and explicit identifiability analysis) rather than hidden by fitting a single mechanism.
- **Cross-study heterogeneity.** Populations differ across trials; captured by a per-trial aggressiveness scalar justified by enrolment differences, with experimental arms left as predictions.
- **PK abstraction.** Exposures are effective potency ratios, not a full PBPK model; a future extension can couple a mechanistic PK layer and intermittent dosing.
- **No within-patient evolution.** Resistance states are fixed per simulation; a stochastic clonal-evolution extension is a natural next step.

---

## 8. Timeline (indicative, 12 months)
- Months 1–2: model construction and verification (Aim 1).
- Months 3–4: preclinical calibration (Aim 2).
- Months 5–8: virtual population, trial simulation, hypothesis testing (Aim 3).
- Months 9–10: identifiability and global sensitivity (Aim 4).
- Months 11–12: manuscript, open-source release, in-silico trial-design extensions.

*(The accompanying R implementation executes this entire plan; see `report/final_report.md`.)*

---

## 9. References

Evidence for the clinical and preclinical anchors was retrieved from **PubMed**. Key sources:

1. Kalinsky K, et al. *Randomized Phase II Trial of Endocrine Therapy With or Without Ribociclib After Progression on CDK4/6 Inhibition (MAINTAIN).* J Clin Oncol. 2023;41(24):4004–4013. [DOI](https://doi.org/10.1200/JCO.22.02392)
2. Mayer EL, et al. *PACE: A Randomized Phase II Study of Fulvestrant, Palbociclib, and Avelumab After Progression on CDK4/6 Inhibitor and Aromatase Inhibitor.* J Clin Oncol. 2024;42(17):2050–2060. [DOI](https://doi.org/10.1200/JCO.23.01940)
3. Kalinsky K, et al. *Abemaciclib Plus Fulvestrant After Progression on CDK4/6 Inhibition (postMONARCH, Phase III).* J Clin Oncol. 2025;43(9):1101–1112. [DOI](https://doi.org/10.1200/JCO-24-02086)
4. Llombart-Cussac A, et al. *Second-Line Endocrine Therapy With or Without Palbociclib Rechallenge (PALMIRA).* J Clin Oncol. 2025;43(18):2084–2093. [DOI](https://doi.org/10.1200/JCO-24-01865)
5. Ramos-Esquivel A, et al. *Continuation of CDK4/6 Inhibition and Switching of Hormonal Therapy After Progression: Systematic Review and Meta-Analysis.* Cureus. 2024;16(11):e73738. [DOI](https://doi.org/10.7759/cureus.73738)
6. Benvenuti C, et al. *Unveiling the Potential of CDK4/6 Inhibitors Beyond Progression in HR+/HER2- Advanced Breast Cancer.* Curr Treat Options Oncol. 2024;25(12):1517–1537. [DOI](https://doi.org/10.1007/s11864-024-01259-4)
7. Zapatero-Solana E, et al. *Models of Early Resistance to CDK4/6 Inhibitors Unveil Potential Therapeutic Treatment Sequencing.* Int J Mol Sci. 2025;26(6):2643. [DOI](https://doi.org/10.3390/ijms26062643)
8. Algethami M, et al. *Investigation of ribociclib, abemaciclib and palbociclib resistance in ER+ breast cancer cells reveals potential therapeutic opportunities.* Sci Rep. 2025;15(1):28579. [DOI](https://doi.org/10.1038/s41598-025-11052-4)
9. Jeselsohn R, et al. *Circulating tumor DNA mutational landscape after progression on a CDK4/6 inhibitor (PACE).* ESMO Open. 2025;10(8):105506. [DOI](https://doi.org/10.1016/j.esmoop.2025.105506)
10. He W, et al. *WEE1 inhibition delays resistance to CDK4/6 inhibitor and antiestrogen treatment in ER+ MCF7 cells.* NPJ Syst Biol Appl. 2025;11(1):128. [DOI](https://doi.org/10.1038/s41540-025-00604-z)

Additional supporting sources retrieved via PubMed: PACE circulating-tumour-cell dynamics (Gerratana L, et al. Clin Cancer Res. 2025;31(21):4510–4517, [DOI](https://doi.org/10.1158/1078-0432.CCR-25-0327)); abemaciclib rechallenge (Nishimura M, et al. Breast Cancer Res Treat. 2026, [DOI](https://doi.org/10.1007/s10549-026-08031-2)); Models-of-early-resistance treatment sequencing [7].

*Data retrieved from PubMed; please cite the original articles via the DOI links above.*
