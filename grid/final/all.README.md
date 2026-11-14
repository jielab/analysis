# GRID ABM

Source family: `grid.abm.1`. Only completion-verified results under `tmp/abm/<trait>` enter these tables. The historical blocknet / CountHeads pipeline, benchmark fits, and old `final` tables are not merged.

Predictions average every observed training donor inside a training-frozen genetic radius. Binary predictions are the number of donor cases divided by the number of all donors. The global donor-mean fallback is explicitly labelled and evaluated on the full query population. Support does not guarantee correctness. Read actual coverage, case/control retention, ancestry audit, and common-risk overlap together with error. All individual evidence remains local-only.

The separate count_vs_benchmarks table reads only the new independent benchmark COMPLETE contracts. Comparisons require the exact same test IDs and labels (including missingness); no smaller intersection is taken. It evaluates each method on the same full roster and the same count-supported people, with paired group bootstrap. Missing/incompatible benchmarks are explicitly listed, never replaced by old aggregates. Training protocols and budgets remain separately labelled.

```tsv
trait	status	validation_status	grouping_status	gwas_overlap_status
height	complete	exploratory_previously_inspected_test	verified_from_kinship_components	unverified
ldl	complete	exploratory_previously_inspected_test	verified_from_kinship_components	unverified
t2dm	complete	exploratory_previously_inspected_test	verified_from_kinship_components	unverified
```
