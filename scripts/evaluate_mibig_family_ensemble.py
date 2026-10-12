#!/usr/bin/env python3
"""Evaluate protein-homology plus functional-profile MIBiG family models."""

import argparse
import collections
import csv
import gzip
import math
from pathlib import Path


def fasta_lengths(path):
    lengths = {}
    name = None
    n = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


def split_terms(value):
    if not value:
        return set()
    return {x.strip() for x in value.replace("|", ",").split(";") if x.strip()}


def percentile(scores, p):
    values = sorted(scores.values())
    if not values:
        return 0.0
    rank = (len(values) - 1) * p / 100.0
    lo, hi = math.floor(rank), math.ceil(rank)
    if lo == hi:
        return values[lo]
    return values[lo] + (values[hi] - values[lo]) * (rank - lo)


def evaluate(name, label_scores_by_query, truth, families):
    stats = collections.Counter()
    details = []
    for query, truth_family in truth.items():
        scores = label_scores_by_query[query]
        ranked = sorted(scores, key=lambda label: (-scores[label], label))[:5]
        predicted = ranked[0] if ranked else "*"
        truth_terms = set(families[query].split(";"))
        ranked_terms = [set(x.split(";")) for x in ranked]
        stats["queries"] += 1
        stats["hit"] += bool(ranked)
        stats["exact1"] += predicted == truth_family
        stats["exact5"] += any(x == truth_family for x in ranked)
        stats["primary1"] += bool(ranked and ranked_terms[0] & truth_terms)
        stats["primary5"] += any(x & truth_terms for x in ranked_terms)
        truth_broad = set(families[query].split(";")[:2])
        stats["broad1"] += bool(ranked and set(ranked[0].split(";")[:2]) & truth_broad)
        stats["broad5"] += any(
            set(x.split(";")[:2]) & truth_broad for x in ranked
        )
        details.append(
            {
                "query_cluster": query,
                "truth_family": truth_family,
                "predicted_family": predicted,
                "top5_families": ";".join(ranked),
                "score": scores.get(predicted, 0.0),
            }
        )
    n = stats["queries"]
    return {
        "method": name,
        "queries": n,
        "hit_fraction": stats["hit"] / n,
        "exact_top1": stats["exact1"] / n,
        "exact_top5": stats["exact5"] / n,
        "primary_top1": stats["primary1"] / n,
        "primary_top5": stats["primary5"] / n,
        "broad_top1": stats["broad1"] / n,
        "broad_top5": stats["broad5"] / n,
    }, details


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mmseqs", required=True)
    ap.add_argument("--query-fasta", required=True)
    ap.add_argument("--target-fasta", required=True)
    ap.add_argument("--query-features", required=True)
    ap.add_argument("--target-features", required=True)
    ap.add_argument("--family-map", required=True)
    ap.add_argument("--output")
    ap.add_argument("--detail-prefix")
    args = ap.parse_args()

    qlens = fasta_lengths(args.query_fasta)
    tlens = fasta_lengths(args.target_fasta)
    families = {}
    with open(args.family_map, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            families[row["cluster_id"]] = row["family"]
    truth = {
        query.split("|", 1)[0]: families.get(query.split("|", 1)[0], query.split("|", 1)[0])
        for query in qlens
    }

    best_protein = {}
    with open(args.mmseqs) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            query, target = f[0], f[1]
            try:
                bits = float(f[11])
                qstart, qend = int(f[6]), int(f[7])
                tstart, tend = int(f[8]), int(f[9])
            except ValueError:
                continue
            qcov = min(1.0, (qend - qstart + 1) / max(1, qlens.get(query, 1)))
            tcov = min(1.0, (tend - tstart + 1) / max(1, tlens.get(target, 1)))
            weighted = bits * math.sqrt(qcov * tcov)
            old = best_protein.get((query, target))
            if old is None or weighted > old[0]:
                best_protein[(query, target)] = (weighted, bits, qcov, tcov)

    raw_best = collections.defaultdict(float)
    raw_sum = collections.defaultdict(float)
    q_with_hits = collections.defaultdict(set)
    t_with_hits = collections.defaultdict(set)
    for (query, target), (weighted, *_rest) in best_protein.items():
        qcluster = query.split("|", 1)[0]
        tcluster = target.split("|", 1)[0]
        key = (qcluster, tcluster)
        raw_best[key] = max(raw_best[key], weighted)
        raw_sum[key] += weighted
        q_with_hits[key[0]].add(query)
        t_with_hits[key[1]].add(target)

    protein_methods = {
        "protein_weighted_best": dict,
        "protein_gene_normalized": dict,
    }
    protein_scores = {
        "protein_weighted_best": collections.defaultdict(dict),
        "protein_gene_normalized": collections.defaultdict(dict),
    }
    for (qcluster, tcluster), value in raw_best.items():
        protein_scores["protein_weighted_best"][qcluster][tcluster] = value
        normalized = raw_sum[(qcluster, tcluster)] / math.sqrt(
            len(q_with_hits[qcluster]) * len(t_with_hits[tcluster])
        )
        protein_scores["protein_gene_normalized"][qcluster][tcluster] = normalized

    profiles = collections.defaultdict(set)
    for path in (args.query_features, args.target_features):
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                cluster = row["cluster_id"].removeprefix("query__")
                profiles[cluster].update(split_terms(row.get("ko", "")))
                profiles[cluster].update(split_terms(row.get("pfam", "")))
    document_frequency = collections.Counter()
    for terms in profiles.values():
        document_frequency.update(terms)
    total_docs = max(1, len(profiles))
    idf = {term: math.log((total_docs + 1) / (df + 1)) + 1.0 for term, df in document_frequency.items()}

    def weighted_vector(cluster):
        return {term: idf.get(term, 1.0) for term in profiles.get(cluster, set())}

    target_vectors = {cluster: weighted_vector(cluster) for cluster in {
        tcluster for _q, tcluster in raw_best
    } | set(families)}
    profile_scores = collections.defaultdict(dict)
    for query in truth:
        qvec = weighted_vector(query)
        qnorm = math.sqrt(sum(x * x for x in qvec.values()))
        for target, tvec in target_vectors.items():
            if target == query:
                continue
            overlap = sum(weight * tvec.get(term, 0.0) for term, weight in qvec.items())
            tnorm = math.sqrt(sum(x * x for x in tvec.values()))
            score = overlap / (qnorm * tnorm) if qnorm and tnorm else 0.0
            profile_scores[query][target] = score

    methods = dict(protein_scores)
    methods["profile_idf_cosine"] = {
        query: dict(scores) for query, scores in profile_scores.items()
    }
    for protein_weight in (0.8, 0.7, 0.6, 0.5):
        name = f"ensemble_protein{protein_weight:.1f}_profile{1-protein_weight:.1f}"
        methods[name] = collections.defaultdict(dict)
        for query in truth:
            pmax = max(protein_scores["protein_gene_normalized"].get(query, {0: 0}).values())
            fmax = max(profile_scores.get(query, {0: 0}).values())
            targets = set(protein_scores["protein_gene_normalized"].get(query, {})) | set(
                profile_scores.get(query, {})
            )
            for target in targets:
                pscore = protein_scores["protein_gene_normalized"].get(query, {}).get(target, 0.0)
                fscore = profile_scores.get(query, {}).get(target, 0.0)
                pscore = pscore / pmax if pmax else 0.0
                fscore = fscore / fmax if fmax else 0.0
                methods[name][query][target] = protein_weight * pscore + (1 - protein_weight) * fscore

    results = []
    all_details = {}
    for name, scores in methods.items():
        label_scores = collections.defaultdict(dict)
        for query, targets in scores.items():
            for target, score in targets.items():
                label = families.get(target, target)
                current = label_scores[query].get(label, 0.0)
                label_scores[query][label] = max(current, score)
        result, details = evaluate(name, label_scores, truth, families)
        results.append(result)
        all_details[name] = details

    fields = list(results[0])
    out = Path(args.output) if args.output else Path("mibig_family_ensemble_summary.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)
    if args.detail_prefix:
        for name, details in all_details.items():
            detail_path = Path(f"{args.detail_prefix}.{name}.tsv")
            detail_path.parent.mkdir(parents=True, exist_ok=True)
            with open(detail_path, "w", newline="") as fh:
                writer = csv.DictWriter(
                    fh,
                    fieldnames=list(details[0]),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(details)
    for row in results:
        print("\t".join(f"{row[x]:.4f}" if isinstance(row[x], float) else str(row[x]) for x in fields))


if __name__ == "__main__":
    main()
