# Hard contamination filter results

Output of [`ANALYSIS/filter_contaminated_genes.R`](../filter_contaminated_genes.R),
run via [`run_decontaminate.sh`](../../run_decontaminate.sh).

## What it does

Removes human genes that pick up appreciable signal in the rat-brain **Control**
samples (`IL64B`, `N168B`, `N269B`; negligible graft). Because those samples
passed through the identical xengsort + Salmon pipeline, any human-gene counts
they carry are residual cross-species contamination — so a gene "expressed" in
them is a contamination sink and is dropped from the count matrix before
differential expression.

**Rule:** drop a gene if its **CPM ≥ cutoff in ≥ `MIN_CONTROLS` of the controls**
(defaults: cutoff = 1, `MIN_CONTROLS` = 2). CPM uses each sample's total
human-assigned counts as the library size.

This is a *hard* filter (the gene is removed entirely), distinct from the
RUVSeq covariate adjustment in `ANALYSIS/results_ruvseq/` — see Phase 2B in the
top-level README. The two are independent takes on the same contamination
problem; this one drops straight into nf-core/differentialabundance.

## Files

| File | Description |
|------|-------------|
| `decontamination_sweep.csv` | Genes flagged at CPM ≥ 0.5 / 1 / 5, for detection in ≥ 1 / 2 / 3 controls. Use it to pick a cutoff. |
| `contaminated_genes.tsv` | The dropped genes (at the chosen primary cutoff) with per-control CPM, mean control CPM, and mean tumor CPM. |
| `salmon.merged.gene_counts.decontaminated.tsv` | Written next to the original in `results_human_final/star_salmon/` (gitignored data dir) — the filtered counts matrix. |
| `salmon.merged.gene_lengths.decontaminated.tsv` | Matching filtered transcript-length matrix (same genes). |

## Reproduce / retune

```bash
# default (CPM>=1 in >=2/3 controls)
sbatch run_decontaminate.sh
# or pick a different cutoff / min-controls after inspecting the sweep:
sbatch run_decontaminate.sh 5 2     # CPM>=5 in >=2 controls (less aggressive)
```

Then run differential abundance on the cleaned matrix — `run_de_pdx_v3.sh`
already points `--matrix` / `--transcript_length_matrix` at the
`.decontaminated.tsv` files:

```bash
sbatch run_de_pdx_v3.sh
```

Everything downstream (GSEA, gProfiler2, the publication figure) then runs
unchanged on the decontaminated gene set.
