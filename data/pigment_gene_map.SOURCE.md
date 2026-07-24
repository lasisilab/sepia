# `pigment_gene_map.tsv` — provenance

Pigmentation **gene** regions on both builds — the substrate for a **gene-region** analysis that
does not depend on a specific ascertained SNP panel (items B6/B7). Complements the SNP-level
[`pigment_snp_map.tsv`](pigment_snp_map.tsv).

## Built by
[`code/build_pigment_gene_map.py`](../code/build_pigment_gene_map.py) — reproducible; needs internet.

## Sources unioned (gene symbols → **1,004 unique genes**)
| Source | file | genes |
|---|---|---|
| Assembled pigmentation network | `gene_network_nodes.csv` | 168 |
| Baxter 2018 curated | `baxter2018_650_pigmentation_genes.csv` | 635 |
| Bajpai 2023 CRISPR hits | `bajpai2023_crispr_hits.csv` | 169 |
| D'Arcy 2023 disease genes | `darcy2023_S1_disease_genes.csv` | 243 |

`sources` / `n_sources` record which sets each gene came from; `in_network` marks the curated
168-gene backbone.

## Coordinate resolution
Ensembl REST `/lookup/symbol/homo_sapiens` (batch POST) — GRCh38 from `rest.ensembl.org`, GRCh37 from
`grch37.rest.ensembl.org`. Gene-body coordinates (start/end/strand) per build. **965/1004 resolved on
both builds**; 11 hg38-only, 18 hg19-only, 10 unresolved (retired symbols — to reconcile). Total gene-body
span ~88.7 Mb (hg38). Tracers: MC1R chr16, SLC24A5/OCA2 chr15, TYR chr11.

## Columns
`gene, ensembl_id, chrom, start_hg38, end_hg38, start_hg19, end_hg19, strand, biotype, length_hg38,
sources, n_sources, in_network, status, notes`. `status` ∈ {ok, hg38_only, hg19_only, unresolved}.

## Intended use
Build hg19 regions (± an optional flank for proximal regulatory variants), extract **all** SGDP +
archaic variants within them, LD-prune, then the same modern-PCA-then-project-archaics recipe as the
SNP analyses — but driven by whole pigmentation-gene regions rather than hand-picked SNPs.
