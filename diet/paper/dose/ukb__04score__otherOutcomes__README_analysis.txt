Run after foodwas / cohort and score.R. No dependence on global dat or legacy results.
Exposure: dose (compatibility alias diet.dose.sum); scoring method: dose9half; theoretical range 0-9; step 0.5.
Scoring criteria are defined in score.R and stored in dose_score.rds; no alternative scores are analysed here.
Full prepared-cohort quantile cutoffs and SD are fixed across outcomes; Q1 reference; within-outcome median-score trend.
Observed DOSE quartile ranges: 0-2; 2.5-3; 3.5-4.5; 5-9.
Cohort inherits all s0 exclusions (including prevalent overall stroke), then excludes prevalent cases per outcome.
Prostate cancer (C61) is analysed among men only (sex.b = Male); the now-constant sex covariate is omitted and recorded. Other outcomes retain both sexes.
ICD regex and labels are taken unchanged from icd10_other.lst, excluding its duplicate stroke entry.
ICD dates replace same-named cohort dates, especially broad cvd; baseline adjustment covariates remain unchanged.
Matching ICD records without dates are excluded from incident analysis as timing-unknown.
Censoring: earliest loss, death or administrative date 2023-04-01; events must be after baseline and on/before censoring.
Stroke follow-up is rebuilt too; inspect stroke_followup_comparison_QC.csv for differences from the existing helper.
M1 only, matching s3_2 for DOSE. Alcohol in score: TRUE; separate alcohol adjustment is omitted when true. No regional strata or PCs added.
Requested M1 covariates: age, sex.b, ethnic.c3, edu.c, tdi, energy.kJ, smoke.c, days_pa_mod.c, bmi.
Complete cases are selected using only M1 covariates and DOSE within each outcome; continuous, quartile and trend models share that cohort.
Constant covariates after exclusions are omitted and recorded; other singularities/convergence failures yield missing estimates.
FDR families are across all planned outcomes within each model/scale or model/category; tests adjusted separately.
Tables excluding stroke retain the same full outcome-wide FDR family; no recalculation on selected subsets.
Q1 is a reference convention; consult Events, Group_N and model_fit_status before interpreting sparse outcomes.
Outcome-wide associations are exploratory. The current cohort was also used for DOSE development.
MAKE_FIGURES controls circular bars and dot-and-whisker plots in PDF and PNG for both scales. Circular plots use the Fig2 MYTH palette with stroke kept coral red; dotplots retain the reverse reference palette. Use --circular-only to refresh circular plots from saved results. All analysis and plotting code is contained in score.R.
