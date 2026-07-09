#!/usr/bin/env python3
"""Evaluate min_tags threshold impact on ABHQ (real draft genome)."""

import subprocess
import re


def run_scaffold(min_tags, ref, draft):
    cmd = [
        "./target/release/syn2b", "scaffold",
        "-r", ref, "-d", draft, "-o", f"/tmp/scaffold_abhq_{min_tags}.agp",
        "--min-tags", str(min_tags)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd="/Users/shihuang/Downloads/iibSyn")
    output = result.stdout + result.stderr

    anchored = 0
    for line in output.split('\n'):
        if "Anchored" in line and "contig(s)" in line:
            match = re.search(r'Anchored (\d+) contig\(s\)', line)
            if match:
                anchored = int(match.group(1))
    return anchored


ref = "/tmp/syn2b_benchmark/tgts/E_coli_K-12_MG1655.tgt"
draft = "/tmp/syn2b_benchmark/tgts/ASM17195v1_multicontig_v2.tgt"
total_contigs = 135  # Total contigs in ABHQ

print(f"{'min_tags':>8} | {'anchored':>8} | {'sensitivity':>11}")
print("-" * 35)

for min_tags in [1, 2, 3, 5, 10, 15, 20, 25, 30, 40, 50]:
    anchored = run_scaffold(min_tags, ref, draft)
    sensitivity = anchored / total_contigs
    print(f"{min_tags:>8} | {anchored:>8} | {sensitivity:>10.2%}")
