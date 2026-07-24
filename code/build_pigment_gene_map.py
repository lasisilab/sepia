#!/usr/bin/env python3
# =====================================================================================
# SEPIA B6/B7 — pigmentation GENE map: resolve the gene-network gene sets to genomic
# regions on BOTH builds, so we can analyze the variation across whole pigmentation-gene
# regions (less reliant on a specific ascertained SNP panel).
#
# Unions the gene-symbol sets from tinalasisi/pigmentation-gene-network:
#   assembled network (168) · Baxter (635 curated) · Bajpai (169 CRISPR) · D'Arcy (243 disease)
# and resolves each symbol to its GRCh38 + GRCh37 gene region (chr, start, end, strand,
# Ensembl ID, biotype) via Ensembl REST /lookup/symbol (batch POST).
#
# Output: data/pigment_gene_map.tsv (+ .SOURCE.md). Laptop-side (needs internet).
# Usage: python3 code/build_pigment_gene_map.py [path-to-pigmentation-gene-network]
# =====================================================================================
import csv, json, os, sys, time, urllib.request, collections

GN   = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/GitHub/pigmentation-gene-network")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(REPO, "data", "pigment_gene_map.tsv")
CHROMS = {str(c) for c in range(1, 23)} | {"X", "Y", "MT"}

def load(path):
    with open(path, newline="") as f: return list(csv.DictReader(f))

# ---- 1. union the gene-symbol sets (with source membership) --------------------------
SOURCES = [
    ("network", "data/processed/gene_network_nodes.csv",              "gene"),
    ("baxter",  "data/processed/baxter2018_650_pigmentation_genes.csv","Human gene symbol"),
    ("bajpai",  "data/processed/bajpai2023_crispr_hits.csv",           "Symbol"),
    ("darcy",   "data/processed/darcy2023_S1_disease_genes.csv",       "gene"),
]
genes = collections.defaultdict(set)   # symbol -> set(sources)
for tag, rel, col in SOURCES:
    p = os.path.join(GN, rel)
    if not os.path.exists(p): print("  skip (missing):", rel); continue
    n = 0
    for row in load(p):
        s = (row.get(col) or "").strip()
        if s and not s.startswith(("AC0", "RP1", "CTD-")):   # drop clone/contig-style non-symbols
            genes[s].add(tag); n += 1
    print(f"  {tag}: {n} rows -> {rel}")
symbols = sorted(genes)
print(f"union: {len(symbols)} unique gene symbols")

# ---- 2. resolve each symbol to a gene region on GRCh38 + GRCh37 ----------------------
def lookup(server, syms):
    out = {}
    for i in range(0, len(syms), 900):
        chunk = syms[i:i+900]
        body = json.dumps({"symbols": chunk}).encode()
        req = urllib.request.Request(f"{server}/lookup/symbol/homo_sapiens", data=body, method="POST",
              headers={"Content-Type": "application/json", "Accept": "application/json"})
        for attempt in range(5):
            try:
                with urllib.request.urlopen(req, timeout=120) as r: out.update(json.load(r)); break
            except Exception as e: print(f"  retry {attempt+1} ({server}): {e}"); time.sleep(5*(attempt+1))
        print(f"  {server.split('//')[1].split('.')[0]}: {len(out)}/{len(syms)}"); time.sleep(1)
    return out

print("resolving GRCh38..."); g38 = lookup("https://rest.ensembl.org", symbols)
print("resolving GRCh37..."); g37 = lookup("https://grch37.rest.ensembl.org", symbols)

# ---- 3. write the gene map -----------------------------------------------------------
cols = ["gene","ensembl_id","chrom","start_hg38","end_hg38","start_hg19","end_hg19","strand",
        "biotype","length_hg38","sources","n_sources","in_network","status","notes"]
rows, stats = [], collections.Counter()
for s in symbols:
    a = g38.get(s, {}); b = g37.get(s, {})
    chrom = a.get("seq_region_name") or b.get("seq_region_name") or ""
    ok38 = a.get("start") and (a.get("seq_region_name") in CHROMS)
    ok19 = b.get("start") and (b.get("seq_region_name") in CHROMS)
    status = "ok" if (ok38 and ok19) else ("hg38_only" if ok38 else ("hg19_only" if ok19 else "unresolved"))
    stats[status] += 1
    ln = (a.get("end",0)-a.get("start",0)) if ok38 else ""
    rows.append({"gene": s, "ensembl_id": a.get("id") or b.get("id") or "",
        "chrom": chrom if chrom in CHROMS else "",
        "start_hg38": a.get("start","") if ok38 else "", "end_hg38": a.get("end","") if ok38 else "",
        "start_hg19": b.get("start","") if ok19 else "", "end_hg19": b.get("end","") if ok19 else "",
        "strand": a.get("strand") or b.get("strand") or "", "biotype": a.get("biotype") or b.get("biotype") or "",
        "length_hg38": ln, "sources": ";".join(sorted(genes[s])), "n_sources": len(genes[s]),
        "in_network": "yes" if "network" in genes[s] else "", "status": status,
        "notes": "" if status=="ok" else "not on a primary assembly / symbol unresolved"})

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
print(f"\nwrote {OUT}: {len(rows)} genes"); print("  status:", dict(stats))
tot_mb = sum(r["length_hg38"] for r in rows if isinstance(r["length_hg38"], int))/1e6
print(f"  total gene-body span (hg38): {tot_mb:.1f} Mb")
for tr in ("MC1R","SLC24A5","OCA2","TYR"):
    m = next((r for r in rows if r["gene"]==tr), None)
    print(f"  TRACER {tr}: chr{m['chrom']} hg38 {m['start_hg38']}-{m['end_hg38']} hg19 {m['start_hg19']}-{m['end_hg19']}" if m else f"  {tr}: absent")
