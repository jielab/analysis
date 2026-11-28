Run f/score.R -> score.R -> score.R.
Legacy source: D:/projects/dose_v1/scripts/03_plot_mutiPattern.R.
Included diets, labels and raw score variables follow pattern_spec.csv; outcomes follow analysis_config.rds.
The primary score column and plot label are dose / DOSE. Range and scoring metadata come from pattern_spec.csv. All grouping figures use fitted Q1-Q4.
Primary outcome is ischemic stroke; other stroke outcomes are supporting analyses.
Models forest and consistency scatter each contain separate P10-P90 and per-SD pages.
Quartile PDFs contain one outcome per page; DOSE-DASH combines quartiles and both continuous scales.
Trend figure shows point estimates only; see quartile forest plots for 95% CIs.
Consistency scatter is descriptive, not a statistical heterogeneity test; above y=x means lower IS HR.
Association significance/FDR is read unchanged from s3_2 tables, not recalculated across outcomes.
Correlation uses raw score columns from dat.pattern.rds, pairwise complete observations.
Correlation clustering: Euclidean distance between correlation rows, complete linkage (legacy pheatmap defaults).
Correlation heatmap stars use raw P; exported F matrix is BH-adjusted over unique off-diagonal pairs.
The ten-pattern appendix follows the correlation heatmap, reuses its Spearman matrix, and exports PDF/PNG with Arial text plus 45 pairwise correlations.
Current association ranking and apparent C-index do not establish independently validated superiority.

Fig_Pattern_Stroke_Forest_M2.pdf
Fig_DOSE_Quartiles_RCS_M2.pdf
Fig_Pattern_IS_Comparison_M2.pdf
Fig_Pattern_Stroke_Models_M1_M3.pdf
Fig_Quartiles_dose_M1_M3.pdf
Fig_Quartiles_ahei_M1_M3.pdf
Fig_Quartiles_medi24_M1_M3.pdf
Fig_Quartiles_mind_M1_M3.pdf
Fig_Quartiles_dash_M1_M3.pdf
Fig_Quartiles_hlfd_M1_M3.pdf
Fig_Quartiles_phdi_M1_M3.pdf
Fig_Quartiles_redii_M1_M3.pdf
Fig_Quartiles_digm_M1_M3.pdf
Fig_Quartiles_maha_M1_M3.pdf
Fig_Quartiles_hpdi_M1_M3.pdf
Fig_Quartiles_modern_M1_M3.pdf
Fig_Quartiles_fds_M1_M3.pdf
Fig_Quartiles_hlcd_M1_M3.pdf
Fig_DOSE_DASH_Quartiles_M2.pdf
Fig_Pattern_Quartile_Trends_M2.pdf
Fig_Pattern_Correlation_Cluster.pdf
Fig_Pattern_Correlations_10patterns.pdf
Fig_Pattern_Correlations_10patterns.png
