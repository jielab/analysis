Flow 1: DOSE -> M linear regression; BH-FDR<0.05.
Flow 2: M -> incident IS Cox regression adjusted for M2 but not DOSE; raw P<0.05.
Candidates are the intersection of flow-1 FDR<0.05 and flow-2 P<0.05.
Formal LHR and MODERN outcome models contain both dose_z and mediator_z plus M2.
LHR uses lm + survreg + mediation::mediate; no dist argument is specified.
MODERN uses lm + logistic regression + mediation::mediate.
The Cox product method uses a*b with the delta/Sobel standard error.
For the Cox product method, proportion=(a*b)/c. Its SE, 95% CI, P value and BH-FDR use a first-order delta approximation including SE(a), SE(b) and SE(c), while assuming cross-model covariances are zero.
Cox product proportions are exploratory and may be unstable when the total-effect coefficient c is near zero or imprecisely estimated.
sims=100 and seed=12345 unless settings in the calling script are changed.
DOSE contains Wine; separate alcohol covariates are excluded from M2.
met.rds was previously cleaned by clean_biom(), which does not transform concentrations.
Only direct absolute-concentration metabolites listed in met.lst are analyzed.
Every analyzed metabolite is transformed by log(x+1) and then Z-standardized within its complete-case sample.
Metabolite sampling may precede WebQ dietary assessment; interpret pathways as exploratory.
