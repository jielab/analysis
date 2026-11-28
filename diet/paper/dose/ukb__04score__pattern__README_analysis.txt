Primary DOSE and 13 established dietary patterns: incident ischemic stroke.
Oxford WebQ MedDiet (medi24) is included; baseline touchscreen MedDiet-original (medito) is excluded from the comparison and common-cohort eligibility.
Exposure: dose (compatibility alias diet.dose.sum); scoring method: dose9half; theoretical range 0-9; step 0.5.
Primary scoring rules are defined only in score.R and stored in dose_score.rds.
Alternative totals remain in the scoring RDS; only dose enters the current association analyses.
All 14 scores use the same complete-case participants and outcome. Missing established scores are not imputed here.
Quartiles: R quantile type 7; cut right=TRUE, include.lowest=TRUE. Identical scores remain together.
Continuous HRs: P10-P90 increment and one SD from a single linear Cox fit. P-trend uses within-group medians.
BH FDR families: all 14 scores within each model/scale; category contrasts separately by Q2/Q3/Q4.
Primary DOSE RCS: M2, knots P10/P50/P90 type 2, median reference, likelihood-ratio overall/nonlinear tests.
DOSE contains alcohol: TRUE. Alcohol adjustment follows the score metadata; alcohol-containing scores omit separate alcohol adjustment.
Common C-index covariates: age;sex.b;ethnic.c3;edu.c;tdi;energy.kJ;smoke.c;days_pa_mod.c;bmi;fhist_stroke;primary_ahtn;t2dm. Paired IJ covariance accounts for correlated predictions.
C-index, AIC and score comparisons are apparent/in-sample; this is not independent validation of food selection or cutoffs.
Do not select a scoring method solely to force the lowest Q4 HR or P-nonlinear>0.05. No cutoff search is performed.
Primary summary: DOSE_IS_M2_summary.csv. Historical DOSE_versions_* files, if present, belong to the earlier version comparison and are not updated or read by these scripts.
Primary curves/tests: DOSE_RCS_*.csv. All model and plot labels use DOSE.
Full results: pattern_IS_M2_comparison.csv, pattern_categories.csv, pattern_overall_tests.csv, pattern_IS_prediction_common_M2.csv.
Checkpoints require matching input hashes, grouping, model specifications and fitting functions.
