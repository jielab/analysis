Observed MRI outcomes and ICV; no KNN and no blanket across-modality missingness filter.
139 volume outcomes retained; 66 local DK thickness fields available in BIG40 (temporal poles absent); two global mean fields excluded.
27 FA + 27 MD outcomes and WMH; all outcomes directly scale()-standardized, without log transformation. Do not confuse absent fields with nonsignificant results.
M2 with age, age squared and centre from: baseline age/center (user requested); NOT imaging visit ; ICV for volume/WMH only.
A single observed centre is explicitly recorded as site_constant=1; its intercept absorbs this constant term.
Primary FDR within exposure and prespecified modality family; broader corrections provided.
Component predictors are scores, not intake amounts. Per-score-point beta also reported.
Chronology: not verified: diet/MRI dates not both supplied
HC3: FALSE
Associations only; no claim of MRI mediation of stroke. Input cohort may already have imputed covariates.
No automated outlier trimming; inspect MRI QC and residual diagnostics before manuscript use.
