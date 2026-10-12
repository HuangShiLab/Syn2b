#!/usr/bin/env python3
"""Summarize OAPGC-vs-MIBiG MMseqs2 hits and integrate best orthologues."""

import argparse
import array
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


def read_fasta_lengths(path):
    lengths = {}
    name = None
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name = line[1:].rstrip().split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


def percentile(values, p):
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    rank = (len(ordered) - 1) * p / 100.0
    lo = math.floor(rank)
    hi = math.ceil(rank)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (rank - lo)


def feature_key_from_m8_query(query):
    f = query.split("|")
    if len(f) < 4:
        return None
    return (f[0], f[1], f[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mmseqs", required=True)
    ap.add_argument("--query-fasta", required=True)
    ap.add_argument("--target-fasta", required=True)
    ap.add_argument("--antismash-features", required=True)
    ap.add_argument("--mibig-family-map", required=True)
    ap.add_argument("--best-hits", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--integrated-features", required=True)
    args = ap.parse_args()

    print("Reading FASTA lengths", flush=True)
    qlens = read_fasta_lengths(args.query_fasta)
    tlens = read_fasta_lengths(args.target_fasta)

    feature_cluster = {}
    feature_product = {}
    print("Reading antiSMASH feature table", flush=True)
    with open(args.antismash_features, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            locus = row["gene_id"].rsplit("|", 1)[-1]
            key = (row["genome_id"], row["contig_id"], locus)
            feature_cluster[key] = row["cluster_id"]
            feature_product[key] = row.get("product", "")

    families = {}
    with open(args.mibig_family_map, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            families[row["cluster_id"]] = row["family"]

    print("Selecting MMseqs2 best hits", flush=True)
    best = {}
    raw_hits = 0
    with open(args.mmseqs) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            raw_hits += 1
            query, target = f[0], f[1]
            try:
                pident = float(f[2]) * 100.0
                alnlen = int(f[3])
                evalue = float(f[10])
                bits = float(f[11])
                qspan = int(f[7]) - int(f[6]) + 1
                tspan = int(f[9]) - int(f[8]) + 1
            except ValueError:
                continue
            qlen = qlens.get(query, 0)
            tlen = tlens.get(target, 0)
            score = (bits, pident, alnlen)
            old = best.get(query)
            if old is None or score > old[:3]:
                accession = target.split("|", 1)[0]
                best[query] = (
                    bits,
                    pident,
                    alnlen,
                    evalue,
                    qspan / qlen * 100.0 if qlen else 0.0,
                    tspan / tlen * 100.0 if tlen else 0.0,
                    target,
                    accession,
                    families.get(accession, ""),
                )

    fields = [
        "query_protein",
        "query_genome",
        "query_contig",
        "query_locus",
        "query_cluster_id",
        "mibig_accession",
        "mibig_protein",
        "mibig_family",
        "bitscore",
        "identity",
        "alignment_length",
        "query_coverage",
        "target_coverage",
        "evalue",
    ]
    print(f"Writing {len(best):,} best hits", flush=True)
    ids = []
    best_by_feature = {}
    identities = array.array("d")
    qcovs = array.array("d")
    tcovs = array.array("d")
    thresholds = Counter()
    region_set = set()
    accession_set = set()
    family_set = set()
    with open(args.best_hits, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for query in sorted(best):
            hit = best[query]
            parts = query.split("|")
            genome = parts[0] if parts else ""
            contig = parts[1] if len(parts) > 1 else ""
            locus = parts[-1] if len(parts) > 1 else ""
            cluster_id = feature_cluster.get((genome, contig, locus), "")
            row = {
                "query_protein": query,
                "query_genome": genome,
                "query_contig": contig,
                "query_locus": locus,
                "query_cluster_id": cluster_id,
                "mibig_accession": hit[7],
                "mibig_protein": hit[6],
                "mibig_family": hit[8],
                "bitscore": f"{hit[0]:.3f}",
                "identity": f"{hit[1]:.4f}",
                "alignment_length": hit[2],
                "query_coverage": f"{hit[4]:.3f}",
                "target_coverage": f"{hit[5]:.3f}",
                "evalue": f"{hit[3]:.3g}",
            }
            writer.writerow(row)
            ids.append(query)
            identities.append(hit[1])
            qcovs.append(hit[4])
            tcovs.append(hit[5])
            thresholds["identity_ge_50"] += hit[1] >= 50.0
            thresholds["identity_ge_70"] += hit[1] >= 70.0
            thresholds["identity_ge_90"] += hit[1] >= 90.0
            thresholds["strict_orthologue_ge70_qcov50_tcov50"] += (
                hit[1] >= 70.0 and hit[4] >= 50.0 and hit[5] >= 50.0
            )
            thresholds["strict_orthologue_ge90_qcov80_tcov80"] += (
                hit[1] >= 90.0 and hit[4] >= 80.0 and hit[5] >= 80.0
            )
            if cluster_id:
                region_set.add(cluster_id)
            best_by_feature[(genome, contig, locus)] = hit
            accession_set.add(hit[7])
            if hit[8]:
                family_set.add(hit[8])

    summary = {
        "mmseqs_raw_hits": raw_hits,
        "queries_in_fasta": len(qlens),
        "mapped_query_proteins": len(best),
        "unmapped_query_proteins": len(qlens) - len(best),
        "unique_oapgc_regions_with_hits": len(region_set),
        "unique_mibig_accessions_hit": len(accession_set),
        "unique_mibig_families_hit": len(family_set),
        "identity_threshold_counts": dict(thresholds),
        "identity_threshold_fractions": {
            k: v / len(best) if best else 0.0 for k, v in thresholds.items()
        },
        "identity_percentiles": {
            str(p): percentile(identities, p) for p in (5, 25, 50, 75, 90, 95, 99)
        },
        "query_coverage_percentiles": {
            str(p): percentile(qcovs, p) for p in (5, 25, 50, 75, 90, 95, 99)
        },
        "target_coverage_percentiles": {
            str(p): percentile(tcovs, p) for p in (5, 25, 50, 75, 90, 95, 99)
        },
        "bitscore_percentiles": {
            str(p): percentile([best[x][0] for x in ids], p)
            for p in (5, 25, 50, 75, 90, 95, 99)
        },
    }
    with open(args.summary, "w") as out:
        json.dump(summary, out, indent=2)
        out.write("\n")

    print("Writing integrated feature table", flush=True)
    out_fields_extra = [
        "mibig_orthologue",
        "mibig_accession",
        "mibig_family",
        "mibig_identity",
        "mibig_query_coverage",
        "mibig_target_coverage",
        "mibig_bitscore",
    ]
    tmp = str(args.integrated_features) + ".tmp"
    with open(args.antismash_features, newline="") as src, open(
        tmp, "w", newline=""
    ) as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.DictWriter(
            dst,
            fieldnames=reader.fieldnames + out_fields_extra,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        n_integrated = 0
        for row in reader:
            locus = row["gene_id"].rsplit("|", 1)[-1]
            query = best_by_feature.get((row["genome_id"], row["contig_id"], locus))
            if query is not None:
                row.update(
                    {
                        "mibig_orthologue": query[6],
                        "mibig_accession": query[7],
                        "mibig_family": query[8],
                        "mibig_identity": f"{query[1]:.4f}",
                        "mibig_query_coverage": f"{query[4]:.3f}",
                        "mibig_target_coverage": f"{query[5]:.3f}",
                        "mibig_bitscore": f"{query[0]:.3f}",
                    }
                )
                n_integrated += 1
            else:
                row.update({x: "" for x in out_fields_extra})
            writer.writerow(row)
    Path(tmp).replace(args.integrated_features)
    summary["integrated_antismash_feature_rows"] = n_integrated
    with open(args.summary, "w") as out:
        json.dump(summary, out, indent=2)
        out.write("\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
