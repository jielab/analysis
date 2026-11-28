Metascape workflow (not Enrichr)
1. Submit candidates_ensembl.txt to https://metascape.org/; choose Homo sapiens as input and analysis species.
Candidate definition: all FUMA mapped genes union MAGMA bonferroni P < 0.05 genes.
2. Record run date, database versions, background, submitted and successfully recognized gene counts in the job records.
Use the documented Metascape default genome background unless specifying a justified custom assay universe; record your choice.
3. Download the complete result ZIP and extract it under metascape/results.
4. Rerun this script. FLT_METASCAPE_RESULTS can point to Enrichment_GO/_FINAL_GO.csv.
The import uses native Log(q-value), Enrichment and #GeneInGOAndHitList values.
It selects the lowest-q significant term per GroupID (ties broken by raw P), then plots up to 20 groups.
The script does not submit a web job or substitute Enrichr for Metascape.
Do not use Excel group-summary rows as term-level enrichment results.
References: https://metascape.org/blog/?p=122 ; https://metascape.org/blog/?p=252
