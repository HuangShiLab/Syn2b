#!/usr/bin/env python3
"""Fragmentation benchmark for Syn2b-Contig MVP.

The reference may be an existing TGT directory. Query complete genomes are
fragmented, digested one fragment at a time (because one FASTA is one TGT
genome in `syn2b digest`), concatenated into a multi-record text TGT, placed,
and evaluated against the genome identity encoded in the fragment header.

Example:
  scripts/benchmark_contig_fragmentation.py \
    --reference-tgt-dir tgt/refdb \
    --query-fasta complete_genomes/*.fna \
    --syn2b target/release/syn2b \
    --output benchmark_contig_fragmentation.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

LEVELS = [500_000, 250_000, 100_000, 50_000, 25_000, 10_000]


def run(cmd: list[str], capture: bool = False) -> subprocess.CompletedProcess:
    print("+", " ".join(map(str, cmd)), file=sys.stderr)
    return subprocess.run(list(map(str, cmd)), check=True, capture_output=capture, text=capture)


def fasta_records(path: Path):
    genome_id = path.name
    for suffix in (".fna.gz", ".fna", ".fa.gz", ".fa"):
        if genome_id.endswith(suffix):
            genome_id = genome_id[: -len(suffix)]
            break
    name = None
    seq: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield genome_id, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            elif line:
                seq.append(line)
    if name is not None:
        yield genome_id, "".join(seq)


def fragment_fastas(paths: list[Path], sizes: list[int], destination: Path) -> list[dict]:
    jobs = []
    for source in paths:
        for genome_id, sequence in fasta_records(source):
            for size in sizes:
                if len(sequence) < size:
                    continue
                count = (len(sequence) + size - 1) // size
                for index in range(count):
                    fragment = sequence[index * size:(index + 1) * size]
                    if len(fragment) < 100:
                        continue
                    query_id = f"frag__{genome_id}__{size}__{index:06d}"
                    fragment_path = destination / f"{query_id}.fna"
                    fragment_path.write_text(f">{query_id}\n{fragment}\n")
                    jobs.append({
                        "query_id": query_id,
                        "truth": genome_id,
                        "level": size,
                        "fasta": fragment_path,
                    })
    return jobs


def digest_queries(jobs: list[dict], syn2b: Path, destination: Path, mode: str,
                   enzymes: str) -> Path:
    query_tgt = destination / "fragments.tgt"
    with query_tgt.open("wb") as combined:
        for job in jobs:
            fragment_tgt = destination / f"{job['query_id']}.tgt"
            run([syn2b, "digest", "--input", job["fasta"], "--output", fragment_tgt,
                 "--format", "text", "--mode", mode, "--enzymes", enzymes], capture=True)
            combined.write(fragment_tgt.read_bytes())
    return query_tgt


def build_reference(args: argparse.Namespace, work: Path) -> Path:
    if args.reference_tgt_dir:
        return Path(args.reference_tgt_dir)
    reference_dir = work / "reference_tgt"
    reference_dir.mkdir(parents=True)
    for source in map(Path, args.reference_fasta):
        output = reference_dir / f"{source.stem}.tgt"
        run([args.syn2b, "digest", "--input", source, "--output", output,
             "--format", "text", "--mode", args.mode, "--enzymes", args.enzymes], capture=True)
    return reference_dir


def read_clusters(path: Path | None) -> dict[str, str]:
    clusters = {}
    if not path:
        return clusters
    with path.open() as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        for row in rows:
            clusters[row["genome_id"]] = row.get("species_cluster") or row["genome_id"]
    return clusters


def evaluate(placements_path: Path, hits_path: Path, clusters: dict[str, str],
             query_count: int) -> tuple[dict, list[dict]]:
    hits = defaultdict(list)
    with hits_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hits[row["query_id"]].append(row)

    with placements_path.open() as handle:
        placements = list(csv.DictReader(handle, delimiter="\t"))
    if len(placements) != query_count:
        raise RuntimeError(f"expected {query_count} placements, got {len(placements)}")

    by_level = defaultdict(lambda: {"queries": 0, "top1": 0, "top5": 0,
                                    "species": 0, "orientation": 0,
                                    "coverage_error": 0.0, "confident_correct": 0})
    rows = []
    for placement in placements:
        query_id = placement["query_id"]
        try:
            _, truth, level_text, _index = query_id.split("__", 3)
        except ValueError as exc:
            raise RuntimeError(f"fragment {query_id} does not have a benchmark truth header") from exc
        level = int(level_text)
        ranked = hits.get(query_id, [])
        predicted = ranked[0]["genome"] if ranked else "*"
        predicted_cluster = ranked[0]["species_cluster"] if ranked else "*"
        truth_cluster = clusters.get(truth, truth)
        stats = by_level[level]
        stats["queries"] += 1
        stats["top1"] += int(predicted == truth)
        stats["top5"] += int(any(hit["genome"] == truth for hit in ranked[:5]))
        stats["species"] += int(predicted_cluster == truth_cluster)
        # Fragmentation preserves forward coordinate order; reverse placements are errors.
        stats["orientation"] += int(bool(ranked) and ranked[0]["orientation"] == "+")
        try:
            stats["coverage_error"] += abs(float(placement["query_coverage"]) - 1.0)
        except ValueError:
            pass
        stats["confident_correct"] += int(
            placement["confidence"] in ("high", "medium") and predicted == truth
        )
        rows.append({
            "level": level,
            "query_id": query_id,
            "truth": truth,
            "predicted": predicted,
            "correct": predicted == truth,
            "rank": next((i + 1 for i, hit in enumerate(ranked)
                          if hit["genome"] == truth), "NA"),
            "score": ranked[0]["score"] if ranked else 0.0,
            "confidence": placement["confidence"],
            "query_coverage": placement["query_coverage"],
        })

    summary = []
    for level in sorted(by_level):
        stats = by_level[level]
        n = stats["queries"]
        summary.append({
            "level_bp": level,
            "queries": n,
            "top1_accuracy": stats["top1"] / n,
            "top5_accuracy": stats["top5"] / n,
            "species_cluster_accuracy": stats["species"] / n,
            "orientation_accuracy": stats["orientation"] / n,
            "mean_coverage_error": stats["coverage_error"] / n,
            "confident_placement_precision": (stats["confident_correct"] / n),
        })
    return summary, rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-tgt-dir")
    ap.add_argument("--reference-fasta", nargs="*", default=[])
    ap.add_argument("--query-fasta", nargs="*", default=[])
    ap.add_argument("--query-list")
    ap.add_argument("--metadata")
    ap.add_argument("--syn2b", default="target/release/syn2b")
    ap.add_argument("--output", required=True)
    ap.add_argument("--summary-output")
    ap.add_argument("--workdir")
    ap.add_argument("--mode", default="2brad")
    ap.add_argument("--enzymes", default="BcgI")
    ap.add_argument("--min-tags", type=int, default=5)
    ap.add_argument("--min-adjacencies", type=int, default=3)
    ap.add_argument("--max-targets", type=int, default=25)
    args = ap.parse_args()

    syn2b = Path(args.syn2b).resolve()
    if not args.reference_tgt_dir and not args.reference_fasta:
        ap.error("provide --reference-tgt-dir or --reference-fasta")

    owned_temp = args.workdir is None
    work = Path(tempfile.mkdtemp(prefix="syn2b_contig_bench.")) if owned_temp else Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    try:
        fragments = work / "fragments"
        fragments.mkdir(exist_ok=True)
        query_fastas = list(map(Path, args.query_fasta))
        if args.query_list:
            with open(args.query_list) as fh:
                query_fastas.extend(
                    Path(line.rstrip("\n")) for line in fh if line.strip()
                )
        jobs = fragment_fastas(query_fastas, LEVELS, fragments)
        if not jobs:
            raise RuntimeError("no fragments met the minimum 100 bp length")

        reference_tgt = build_reference(args, work)
        query_tgt = digest_queries(jobs, syn2b, fragments, args.mode, args.enzymes)
        index = work / "reference.s2cidx"
        metadata = Path(args.metadata) if args.metadata else None
        index_cmd = [syn2b, "contig-index", "--genome-tgt-dir", reference_tgt, "--output", index]
        if metadata:
            index_cmd += ["--genome-metadata", metadata]
        build_start = time.perf_counter()
        run(index_cmd, capture=True)
        build_seconds = time.perf_counter() - build_start

        placements = work / "placements.tsv"
        details = work / "details.tsv"
        hits = placements.with_name(f"{placements.stem}.hits.tsv")
        search_start = time.perf_counter()
        run([syn2b, "place-contig", "--index", index, "--query", query_tgt,
             "--output", placements, "--details", details,
             "--min-tags", args.min_tags, "--min-adjacencies", args.min_adjacencies,
             "--max-targets", args.max_targets], capture=True)
        search_seconds = time.perf_counter() - search_start

        clusters = read_clusters(metadata)
        summary, detail_rows = evaluate(placements, hits, clusters, len(jobs))
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(detail_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(detail_rows)
        summary_path = Path(args.summary_output) if args.summary_output else output.with_suffix(".summary.tsv")
        with summary_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(summary)

        print(f"index_build_seconds\t{build_seconds:.6f}", file=sys.stderr)
        print(f"search_seconds\t{search_seconds:.6f}", file=sys.stderr)
        print(f"queries_per_second\t{len(jobs) / search_seconds:.3f}", file=sys.stderr)
        print(f"wrote {output} and {summary_path}")
        return 0
    finally:
        if owned_temp:
            shutil.rmtree(work)


if __name__ == "__main__":
    raise SystemExit(main())
