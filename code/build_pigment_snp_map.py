#!/usr/bin/env python3
# =====================================================================================
# SEPIA B1/B3 — build the EXTENDED pigmentation SNP map (superset of the 222-SNP panel).
#
# Unions every rsID-bearing pigmentation SNP source in tinalasisi/pigmentation-gene-network
#   - GWAS Catalog pigmentation pull (nb4 union: 1,072 dedup + 105 curated-paper rows)
#   - HIrisPlex-S forensic markers (36)
#   - Martin 2017 KhoeSan GWAS (lead_rsid)
#   - Kim 2024 East-Asian GWAS (rsid)
#   - discordance / case-study loci (rsid)
# then RE-RESOLVES every rsID to verified GRCh38 + GRCh37 coordinates + alleles via Ensembl
# REST (batch POST) -- because the source coordinate columns are mixed-build and untrustworthy
# (same reason the 222-panel map caught the +1 off-by-one). Tags each SNP high/moderate/low
# confidence from replication signals, and flags palindromic / multiallelic / merged rsIDs.
#
# Output: data/pigment_snp_map.tsv (+ .SOURCE.md).  Laptop-side (needs internet for Ensembl).
# Usage:  python3 code/build_pigment_snp_map.py [path-to-pigmentation-gene-network]
# =====================================================================================
import csv, json, os, sys, time, urllib.request, collections

GN   = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/GitHub/pigmentation-gene-network")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(REPO, "data", "pigment_snp_map.tsv")
PANEL = os.path.join(REPO, "data", "panel_hg19_hg38_map.tsv")
CHROMS = {str(c) for c in range(1, 23)} | {"X", "Y"}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def load(path):
    delim = "\t" if path.endswith((".tsv", ".tab")) else ","   # panel map is TSV; sources are CSV
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter=delim))

def is_rs(v):
    return isinstance(v, str) and v.startswith("rs") and v[2:].isdigit()

# ---- 1. gather rsIDs + per-source confidence signals --------------------------------
# rec[rsid] = {sources:set, genes:set, conf:int(1-3), gwas_n_assoc, gwas_replicated, pvalue, notes}
CONF = {"low": 1, "moderate": 2, "high": 3}
rec = collections.defaultdict(lambda: {"sources": set(), "genes": set(), "conf": 0,
                                       "gwas_n_assoc": "", "gwas_replicated": "", "pvalue": "", "merged_from": set()})

def bump(rsid, src, conf, gene="", **extra):
    r = rec[rsid]
    r["sources"].add(src)
    if gene and is_rs(rsid): r["genes"].add(gene)
    r["conf"] = max(r["conf"], CONF[conf])
    for k, v in extra.items():
        if v not in ("", None) and not r.get(k):
            r[k] = v

# nb4 union (GWAS dedup + curated papers) -- carries gwas_replicated / gwas_n_assoc
for row in load(os.path.join(GN, "data/processed/nb4_unified_association_base.csv")):
    rsid = (row.get("rsid_primary") or row.get("rsid") or "").strip()
    if not is_rs(rsid):
        rsid = (row.get("rsid") or "").strip()
    if not is_rs(rsid): continue
    if row.get("rsid") and row["rsid"] != rsid and is_rs(row["rsid"]):
        rec[rsid]["merged_from"].add(row["rsid"])
    st = row.get("source_type", "")
    repl = row.get("gwas_replicated", ""); na = row.get("gwas_n_assoc", "")
    if st == "gwas_catalog_dedup":
        try: high = (repl == "True") or (int(float(na or 0)) >= 2)
        except Exception: high = (repl == "True")
        bump(rsid, "gwas_catalog", "high" if high else "low", gene=row.get("gene_label", ""),
             gwas_n_assoc=na, gwas_replicated=repl, pvalue=row.get("pvalue", ""))
    else:  # curated_paper (case studies / functional literature)
        bump(rsid, "curated:" + (row.get("paper") or "paper"), "moderate", gene=row.get("gene_label", ""),
             pvalue=row.get("pvalue", ""))

# HIrisPlex-S forensic markers -> validated == high confidence
for row in load(os.path.join(GN, "data/processed/hirisplexs2018_markers.csv")):
    if is_rs(row.get("rsid", "")):
        bump(row["rsid"].strip(), "hirisplex_s", "high", gene=row.get("gene", ""))

# Martin 2017 KhoeSan (lead_rsid; genome_wide_significant flag)
for row in load(os.path.join(GN, "data/processed/EXTRACT_Martin2017_loci.csv")):
    rsid = (row.get("lead_rsid") or "").strip()
    if is_rs(rsid):
        gws = str(row.get("genome_wide_significant", "")).lower() in ("true", "yes", "1")
        bump(rsid, "martin2017", "high" if gws else "moderate", gene=row.get("author_attributed_gene", ""))

# Kim 2024 East-Asian (rsid; genome-wide sig from pvalue columns)
for row in load(os.path.join(GN, "data/processed/EXTRACT_Kim2024_loci_v2.csv")):
    rsid = (row.get("rsid") or "").strip()
    if is_rs(rsid):
        gws = False
        for c in ("pvalue_Lstar", "pvalue_astar", "pvalue_bstar"):
            try:
                if float(row.get(c, "1") or 1) < 5e-8: gws = True
            except Exception: pass
        bump(rsid, "kim2024", "high" if gws else "moderate", gene=row.get("nearest_gene_label", ""))

# discordance / case-study loci
for row in load(os.path.join(GN, "data/processed/discordance_loci.csv")):
    rsid = (row.get("rsid") or "").strip()
    if is_rs(rsid):
        bump(rsid, "discordance", "moderate", gene=row.get("nearest_gene_label", ""),
             pvalue=row.get("pvalue", ""))

# 222-panel membership
panel_rsids = {r["rsid"] for r in load(PANEL) if is_rs(r.get("rsid", ""))} if os.path.exists(PANEL) else set()

rsids = sorted(rec.keys(), key=lambda x: int(x[2:]))
print(f"union: {len(rsids)} unique rsIDs from {len(set().union(*[r['sources'] for r in rec.values()]))} source tags")
print("  in 222 panel:", len(set(rsids) & panel_rsids), "| new:", len(set(rsids) - panel_rsids))

# ---- 2. resolve GRCh38 + GRCh37 coordinates + alleles via Ensembl batch POST --------
def ensembl_batch(ids, server):
    url = f"{server}/variation/homo_sapiens"
    body = json.dumps({"ids": ids}).encode()
    req = urllib.request.Request(url, data=body, method="POST",
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.load(resp)
        except Exception as e:
            print(f"  ensembl retry {attempt+1} ({server}): {e}")
            time.sleep(5 * (attempt + 1))
    return {}

def resolve(server, label):
    out = {}
    for i in range(0, len(rsids), 200):
        chunk = rsids[i:i+200]
        data = ensembl_batch(chunk, server)
        for rid, info in (data or {}).items():
            best = None
            for m in info.get("mappings", []):
                if m.get("seq_region_name") in CHROMS and m.get("coord_system") == "chromosome":
                    best = m; break
            if best:
                al = best.get("allele_string", "")
                out[rid] = {"chrom": best["seq_region_name"], "pos": best.get("start"),
                            "alleles": al, "name": info.get("name", rid)}
        print(f"  {label}: resolved {len(out)}/{len(rsids)}")
        time.sleep(1)
    return out

print("resolving GRCh38 (rest.ensembl.org)...")
g38 = resolve("https://rest.ensembl.org", "hg38")
print("resolving GRCh37 (grch37.rest.ensembl.org)...")
g37 = resolve("https://grch37.rest.ensembl.org", "hg19")

# ---- 3. assemble the map ------------------------------------------------------------
def split_alleles(al):
    parts = [p for p in (al or "").split("/") if p]
    ref = parts[0] if parts else ""
    alt = ",".join(parts[1:]) if len(parts) > 1 else ""
    multi = len(parts) > 2
    palin = (len(parts) >= 2 and len(parts[0]) == 1 and len(parts[1]) == 1
             and COMP.get(parts[0]) == parts[1])
    return ref, alt, multi, palin

cols = ["rsid", "rsid_resolved", "chrom", "pos_hg38", "pos_hg19", "ref_hg38", "alt_hg38",
        "ref_hg19", "alt_hg19", "gene", "sources", "n_sources", "confidence",
        "gwas_n_assoc", "gwas_replicated", "palindromic", "multiallelic", "merged",
        "in_222_panel", "status", "notes"]
CONF_NAME = {1: "low", 2: "moderate", 3: "high"}
rows_out, stats = [], collections.Counter()
for rsid in rsids:
    r = rec[rsid]
    a = g38.get(rsid, {}); b = g37.get(rsid, {})
    r38, a38, m38, p38 = split_alleles(a.get("alleles", ""))
    r19, a19, m19, p19 = split_alleles(b.get("alleles", ""))
    resolved = a.get("name") or b.get("name") or rsid
    merged = (resolved != rsid) or bool(r["merged_from"])
    status = "ok" if (a.get("pos") and b.get("pos")) else ("hg38_only" if a.get("pos") else
             ("hg19_only" if b.get("pos") else "unresolved"))
    stats[status] += 1
    gene = ";".join(sorted(g for g in r["genes"] if g)) or (a.get("gene") or "")
    rows_out.append({
        "rsid": rsid, "rsid_resolved": resolved,
        "chrom": a.get("chrom") or b.get("chrom") or "",
        "pos_hg38": a.get("pos") or "", "pos_hg19": b.get("pos") or "",
        "ref_hg38": r38, "alt_hg38": a38, "ref_hg19": r19, "alt_hg19": a19,
        "gene": gene, "sources": ";".join(sorted(r["sources"])), "n_sources": len(r["sources"]),
        "confidence": CONF_NAME.get(r["conf"], "low"),
        "gwas_n_assoc": r["gwas_n_assoc"], "gwas_replicated": r["gwas_replicated"],
        "palindromic": "yes" if (p38 or p19) else "", "multiallelic": "yes" if (m38 or m19) else "",
        "merged": "yes" if merged else "", "in_222_panel": "yes" if rsid in panel_rsids else "",
        "status": status,
        "notes": ("merged_from=" + ",".join(sorted(r["merged_from"]))) if r["merged_from"] else "",
    })

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows_out)

# ---- 4. verify (tracers) + summarize ------------------------------------------------
conf_counts = collections.Counter(row["confidence"] for row in rows_out)
palin = sum(1 for row in rows_out if row["palindromic"]); multi = sum(1 for row in rows_out if row["multiallelic"])
merged = sum(1 for row in rows_out if row["merged"])
print(f"\nwrote {OUT}: {len(rows_out)} SNPs")
print("  status:", dict(stats)); print("  confidence:", dict(conf_counts))
print(f"  palindromic: {palin} | multiallelic: {multi} | merged: {merged} | in-panel: {len(set(rsids)&panel_rsids)}")
for tr, e38, e19 in [("rs1426654", "48134287", "48426484"), ("rs12913832", "28120472", "28365618")]:
    m = next((row for row in rows_out if row["rsid"] == tr), None)
    ok = m and str(m["pos_hg38"]) == e38 and str(m["pos_hg19"]) == e19
    print(f"  TRACER {tr}: hg38={m['pos_hg38'] if m else '-'} (exp {e38}) hg19={m['pos_hg19'] if m else '-'} (exp {e19}) {'OK' if ok else 'CHECK'}")
