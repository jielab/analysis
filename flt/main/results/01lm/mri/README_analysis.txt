Three FLTs fitted separately; primary cohort from s0, no disease-specific eligibility filter.
Observed instance-2 MRI outcomes/ICV only; no outcome or exposure imputation in this script.
260 planned MRI measures: 66 thickness, 27 FA, 27 MD, 139 volumes and 1 WMH.
All outcomes directly standardized per complete-case model, including WMH (no log transform).
beta/SE/LCI/HCI: outcome SD per one FLT point; *_standardized: per one FLT SD.
HC1 robust SE and residual-df two-sided t inference; minimum N=100.
Adjustment covariates: age, sex, bmi, edu, inc, tdi, smoke_status, alcohol_status, pa_met_tot
Other covariates match s1_3; volumes/WMH additionally adjust observed ICV.
Primary FDR_family: within FLT and MRI family; FDR_trait: within MRI phenotype (3 FLTs).
FDR_FLTs_family pools 3 FLTs per modality; FDR_exposure_all pools 260 MRIs per FLT; FDR_all_tests pools all 780 tests.
Failed/missing tests remain in planned multiplicity counts.
Chronology: not verified: FLT/MRI dates not assessed
Associations only; no causal interpretation. Inspect QC and residual diagnostics before publication.
