# allele3 local analysis

Version: 3.2.0-final-20261022

Source: gnomAD-v4.1-genomes; genome build: GRCh38.

- Qualifying SNV positions: 503,781,385
- Positions with at least two qualifying ALT nucleotides: 37,984,306
- Positions with three qualifying ALT nucleotides: 1,549,318
- Multi-ALT positions with different ALT AN values: 15,767,285
- AC thresholds: [1, 2, 5, 10]
- Per-ALT AN floor: 100
- Callable mask: none; proportions are conditional on observed qualifying SNV positions
- Site AC spectrum coverage: complete (22/22 shards)
- Theoretical candidate storage: deferred

The primary quantity is observed non-reference nucleotide occupancy. REF is not
assumed ancestral, major, present, or reconstructable from split AN values.
ALT order is not a mutation timeline. Failed/missing ALT records are not used as
unlisted counterfactuals. Frequencies are each ALT's AC/AN.

The context background uses other scanned autosomes and reports support and
fallback levels. It estimates a conditional observation frequency, not mutation
rate or evolutionary retention. X/Y/MT and regional/single-chromosome scans are
descriptive. All ancestry comparisons remain sensitive to sample size, coverage,
ancestry composition and genotype filtering.

No AlphaGenome network requests were made by this pipeline. Optional offline AVI
annotation is exploratory: AVI used gnomAD frequency proxy labels in training.
This workflow cannot infer haplotype ancestry, allele age, pleiotropy, or balancing
selection from sites-only counts.

Detailed multi-ALT sites and original per-ALT counts are kept in the per-chromosome
shards. With deferred candidate storage, the redundant three-substitution table
is reconstructed losslessly only if full offline annotate is requested. Sparse
Atlas uses its separate bounded plan. Unlisted candidate AC/AN/AF remain missing.

site_ac_spectrum.tsv and singleton_sensitivity.tsv distinguish multi-ALT sites
containing AC=1 ALTs from sites still having two ALTs after AC>=2. The denominator
is the original eligible multi-ALT cohort. This is a count sensitivity, not loss
over evolutionary time. Old caches lacking this spectrum are explicitly unavailable;
mixed cache summaries are partial and cannot be treated as genome-wide zero values.
