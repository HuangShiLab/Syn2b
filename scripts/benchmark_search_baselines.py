#!/usr/bin/env python3
"""Small wrapper for contig-search baseline comparison.

Supported tools: minimap2 (-x asm20), skani search, and mash dist. The script
fragments complete genomes with the same headers as
`benchmark_contig_fragmentation.py`, builds a single concatenated reference
FASTA, records wall time, and writes one ranked-prediction row per query/tool.

This wrapper intentionally uses one thread so serial runtime is comparable
across tools. Production benchmarks may vary threads explicitly.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

from benchmark_contig_fragmentation import LEVELS, fasta_records, fragment_fastas


def require(tool: str) -> bool:
    return shutil.which(tool) is not None


def run(cmd: list[str], output: Path | None = None, stdout=None):
    print("+", " ".join(map(str, cmd)), file=sys.stderr)
    with output.open("wb") if output else open("/dev/null", "wb") as sink:
        subprocess.run(list(map(str, cmd)), stdout=stdout or sink, check=True)


def base_rows(jobs):
    return [{"query_id": job["query_id"], "truth": job["truth"], "level": job["level"]}
            for job in jobs]


def minimap2(reference: Path, query: Path, threads: int) -> dict[str, list[tuple[str, float]]]:
    raw = query.with_suffix(".minimap2.paf")
    run(["minimap2", "-t", threads, "-x", "asm20", reference, query], output=raw)
    scores = defaultdict(float)
    with raw.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 10:
                continue
            query_id, _, _, reference_id = fields[0], fields[1], fields[2], fields[5]
            block_length = int(fields[10])
            scores[(query_id, reference_id)] += block_length
    hits = defaultdict(list)
    for (query_id, reference_id), score in scores.items():
        hits[query_id].append((reference_id, score))
    for values in hits.values():
        values.sort(key=lambda pair: (-pair[1], pair[0]))
    return hits


def skani(reference_dir: Path, query: Path, threads: int) -> dict[str, list[tuple[str, float]]]:
    sketch_dir = reference_dir.parent / "reference.skani"
    raw = query.with_suffix(".skani.tsv")
    fasta_list = reference_dir.parent / "reference fasta list.txt"
    fasta_list.write_text("\n".join(str(path) for path in sorted(reference_dir.iterdir())) + "\n")
    run(["skani", "sketch", "-l", fasta_list, "-o", sketch_dir, "-t", threads])
    run(["skani", "search", "-q", query, "--qi", "-d", sketch_dir, "-t", threads], output=raw)
    hits = defaultdict(list)
    with raw.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 2:
                continue
            # skani v0.2 writes Ref_file, Query_file, ANI, ..., Ref_name, Query_name.
            reference_id, query_id = fields[5], fields[6]
            try:
                score = float(fields[2])
            except (ValueError, IndexError):
                score = 0.0
            hits[query_id].append((reference_id, score))
    for values in hits.values():
        values.sort(key=lambda pair: (-pair[1], pair[0]))
    return hits


def mash(reference: Path, query: Path) -> dict[str, list[tuple[str, float]]]:
    sketch = query.parent / "reference.msh"
    raw = query.with_suffix(".mash.tsv")
    run(["mash", "sketch", "-o", query.parent / "reference", "-s", "1000", reference])
    # `mash sketch -o` appends .msh.
    run(["mash", "dist", sketch, query, "-p", "1"], output=raw)
    scores = defaultdict(dict)
    with raw.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                continue
            reference_id, query_id, distance = fields[0], fields[1], float(fields[2])
            scores[query_id][reference_id] = 1.0 - distance
    hits = defaultdict(list)
    for query_id, values in scores.items():
        hits[query_id] = sorted(values.items(), key=lambda pair: (-pair[1], pair[0]))
    return hits


def write(output: Path, tool: str, jobs, hits, runtime_seconds):
    with output.open("w", newline="") as handle:
        fields = ["tool", "query_id", "truth", "level", "predicted", "correct", "rank",
                  "runtime_seconds"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for job in jobs:
            ranked = hits.get(job["query_id"], [])
            predicted = ranked[0][0] if ranked else "*"
            rank = next((i + 1 for i, pair in enumerate(ranked)
                         if pair[0] == job["truth"]), "NA")
            writer.writerow({
                "tool": tool,
                "query_id": job["query_id"],
                "truth": job["truth"],
                "level": job["level"],
                "predicted": predicted,
                "correct": predicted == job["truth"],
                "rank": rank,
                "runtime_seconds": runtime_seconds,
            })


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-fasta", nargs="+", required=True)
    ap.add_argument("--query-fasta", nargs="+", required=True)
    ap.add_argument("--output-prefix", required=True)
    ap.add_argument("--tools", nargs="*", default=["minimap2", "skani", "mash"])
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--workdir")
    args = ap.parse_args()

    requested = set(args.tools)
    unavailable = {tool for tool in requested if not require(tool)}
    if unavailable:
        print("skipped unavailable tools: " + ", ".join(sorted(unavailable)), file=sys.stderr)
    requested -= unavailable
    if not requested:
        print("no requested baseline tools are installed", file=sys.stderr)
        return 2

    owned = args.workdir is None
    work = Path(tempfile.mkdtemp(prefix="syn2b_search_baselines.")) if owned else Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    try:
        fragments = work / "fragments"
        reference_dir = work / "reference_fasta"
        fragments.mkdir()
        reference_dir.mkdir()
        jobs = fragment_fastas(list(map(Path, args.query_fasta)), LEVELS, fragments)
        if not jobs:
            raise RuntimeError("no fragments met the minimum 100 bp length")
        for source in map(Path, args.reference_fasta):
            shutil.copy2(source, reference_dir / source.name)
        reference = work / "reference.concat.fna"
        with reference.open("w") as out:
            for source in reference_dir.iterdir():
                for name, sequence in fasta_records(source):
                    out.write(f">{name}\n{sequence}\n")
        query = work / "fragments.fna"
        with query.open("w") as out:
            for path in sorted(fragments.glob("*.fna")):
                out.write(path.read_text())

        runtimes = {}
        prefix = Path(args.output_prefix)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        if "minimap2" in requested:
            started = time.perf_counter()
            hits = minimap2(reference, query, args.threads)
            runtimes["minimap2"] = time.perf_counter() - started
            write(prefix.with_suffix(".minimap2.tsv"), "minimap2", jobs, hits, runtimes["minimap2"])
        if "skani" in requested:
            started = time.perf_counter()
            hits = skani(reference_dir, query, args.threads)
            runtimes["skani"] = time.perf_counter() - started
            write(prefix.with_suffix(".skani.tsv"), "skani", jobs, hits, runtimes["skani"])
        if "mash" in requested:
            started = time.perf_counter()
            hits = mash(reference, query)
            runtimes["mash"] = time.perf_counter() - started
            write(prefix.with_suffix(".mash.tsv"), "mash", jobs, hits, runtimes["mash"])
        print(f"wrote baseline tables with prefix {prefix}")
        return 0
    finally:
        if owned:
            shutil.rmtree(work)


if __name__ == "__main__":
    raise SystemExit(main())
