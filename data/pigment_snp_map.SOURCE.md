# `pigment_snp_map.tsv` — provenance

The **extended** pigmentation SNP map: a superset of the 222-SNP panel
([`panel_hg19_hg38_map.tsv`](panel_hg19_hg38_map.tsv)), unioning every rsID-bearing pigmentation
SNP source in [`tinalasisi/pigmentation-gene-network`](https://github.com/tinalasisi/pigmentation-gene-network)
plus HIrisPlex-S, all re-resolved to verified GRCh38 + GRCh37 coordinates.

## Built by
[`code/build_pigment_snp_map.py`](../code/build_pigment_snp_map.py) — reproducible; re-run with
internet access. Reads the gene-network repo (default `~/GitHub/pigmentation-gene-network`).

## Sources unioned (SNP-level, rsID-keyed)
| Source file | Contributes |
|---|---|
| `nb4_unified_association_base.csv` | GWAS Catalog dedup (1,072) + 105 curated-paper rows; carries `gwas_replicated` / `gwas_n_assoc` |
| `hirisplexs2018_markers.csv` | 36 HIrisPlex-S forensic markers (validated → high confidence) |
| `EXTRACT_Martin2017_loci.csv` | KhoeSan (San/Nama) GWAS lead SNPs (`lead_rsid`) |
| `EXTRACT_Kim2024_loci_v2.csv` | East-Asian skin-color GWAS |
| `discordance_loci.csv` | eye-colour discordance / case-study loci |

## Coordinate resolution
Every rsID re-resolved to **both** builds via Ensembl REST batch POST — GRCh38 from
`rest.ensembl.org`, GRCh37 from `grch37.rest.ensembl.org` — keyed by rsID (source coordinate
columns are mixed-build and not trusted). This is the same method that corrected the +1 off-by-one
in the panel's stored hg38 coordinates. Tracers asserted: SLC24A5 `rs1426654` = hg38 `15:48134287` /
hg19 `15:48426484`; HERC2 `rs12913832` = hg38 `15:28120472` / hg19 `15:28365618`.

## Confidence tiers (`confidence` column)
Best evidence across a SNP's sources:
- **high** — GWAS replicated (`gwas_replicated=True` or `gwas_n_assoc ≥ 2`), OR HIrisPlex-S validated, OR a genome-wide-significant Martin/Kim lead SNP.
- **moderate** — curated literature / case-study support, or a genome-wide-significant single association.
- **low** — a single, unreplicated GWAS-Catalog association.

## Columns
`rsid, rsid_resolved, chrom, pos_hg38, pos_hg19, ref_hg38, alt_hg38, ref_hg19, alt_hg19, gene,
sources, n_sources, confidence, gwas_n_assoc, gwas_replicated, palindromic, multiallelic, merged,
in_222_panel, status, notes`.

`status` ∈ {ok, hg38_only, hg19_only, unresolved}. `in_222_panel=yes` marks the original panel
subset. Full narrative + summary counts: the **`snp_panel.qmd`** site page.
