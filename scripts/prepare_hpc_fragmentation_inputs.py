#!/usr/bin/env python3
"""Create deterministic metadata/query lists for HPC Syn2b fragmentation runs."""

import argparse
import gzip
import json
import os
import random
from pathlib import Path


def genome_id_from_path(path):
    name = Path(path).name
    if name.endswith(".fna.gz"):
        return name[:-7]
    if name.endswith(".fna"):
        return name[:-4]
    if name.endswith(".fa.gz"):
        return name[:-6]
    if name.endswith(".fa"):
        return name[:-3]
    return name


def prepare_gtdb(path_list, output, query_count, seed):
    output.mkdir(parents=True, exist_ok=True)
    query_dir = output / "query_fastas"
    query_dir.mkdir(exist_ok=True)
    rows = []
    with open(path_list) as fh:
        for line in fh:
            species, fasta = line.rstrip("\n").split("\t", 1)
            rows.append(
                {
                    "genome_id": genome_id_from_path(fasta),
                    "species_cluster": species,
                    "taxonomy": species,
                    "fasta_path": fasta,
                    "is_representative": 1,
                }
            )
    metadata = output / "gtdb_metadata.tsv"
    with open(metadata, "w") as out:
        out.write("genome_id\tspecies_cluster\ttaxonomy\tfasta_path\tis_representative\n")
        for row in rows:
            out.write(
                "{genome_id}\t{species_cluster}\t{taxonomy}\t{fasta_path}\t{is_representative}\n".format(
                    **row
                )
            )
    selected = sorted(
        random.Random(seed).sample(rows, query_count), key=lambda row: row["genome_id"]
    )
    query_list = output / "query_list.tsv"
    with open(query_list, "w") as out:
        for row in selected:
            source = Path(row["fasta_path"])
            destination = query_dir / f"{row['genome_id']}.fna"
            if destination.exists() or destination.is_symlink():
                destination.unlink()
            if source.suffix == ".gz":
                destination = query_dir / f"{row['genome_id']}.fna.gz"
                if destination.exists() or destination.is_symlink():
                    destination.unlink()
                os.symlink(source.resolve(), destination)
            else:
                os.symlink(source.resolve(), destination)
            out.write(f"{destination.resolve()}\n")
    return {
        "reference_count": len(rows),
        "query_count": len(selected),
        "metadata": str(metadata),
        "query_list": str(query_list),
    }


def prepare_hrom(metadata_path, taxonomy_path, tgt_dir, representative_fasta_dir,
                 output, query_count, seed):
    output.mkdir(parents=True, exist_ok=True)
    query_dir = output / "query_fastas"
    query_dir.mkdir(exist_ok=True)
    taxonomy = {}
    with open(taxonomy_path) as fh:
        for line in fh:
            if line.startswith("accession\t"):
                continue
            accession, value = line.rstrip("\n").split("\t", 1)
            taxonomy[accession] = value
    fasta_by_id = {
        path.stem: path for path in Path(representative_fasta_dir).glob("*.fna")
    }
    rows = []
    query_candidates = []
    with open(metadata_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        required = {"genome", "cluster", "length_bp", "quality"}
        if not required.issubset(header):
            raise RuntimeError(f"HROM metadata missing columns: {required - set(header)}")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            genome = row["genome"]
            tgt = Path(tgt_dir) / f"{genome}.tgt"
            if not tgt.exists():
                continue
            fasta = fasta_by_id.get(genome)
            parsed = {
                "genome_id": genome,
                "species_cluster": row["cluster"],
                "taxonomy": taxonomy.get(genome, row["cluster"]),
                "fasta_path": str(fasta) if fasta else "",
                "is_representative": int(fasta is not None),
            }
            rows.append(parsed)
            if fasta is not None:
                query_candidates.append(parsed)
    metadata_out = output / "hrom_metadata.tsv"
    with open(metadata_out, "w") as out:
        out.write("genome_id\tspecies_cluster\ttaxonomy\tfasta_path\tis_representative\n")
        for row in rows:
            out.write(
                "{genome_id}\t{species_cluster}\t{taxonomy}\t{fasta_path}\t{is_representative}\n".format(
                    **row
                )
            )
    selected = sorted(
        random.Random(seed).sample(query_candidates, query_count),
        key=lambda row: row["genome_id"],
    )
    query_list = output / "query_list.tsv"
    with open(query_list, "w") as out:
        for row in selected:
            destination = query_dir / f"{row['genome_id']}.fna"
            if destination.exists() or destination.is_symlink():
                destination.unlink()
            os.symlink(Path(row["fasta_path"]).resolve(), destination)
            out.write(f"{destination.resolve()}\n")
    return {
        "reference_count": len(rows),
        "query_count": len(selected),
        "metadata": str(metadata_out),
        "query_list": str(query_list),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtdb-list", required=True)
    ap.add_argument("--hrom-metadata", required=True)
    ap.add_argument("--hrom-taxonomy", required=True)
    ap.add_argument("--hrom-tgt-dir", required=True)
    ap.add_argument("--hrom-fasta-dir", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--gtdb-queries", type=int, default=1000)
    ap.add_argument("--hrom-queries", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20261012)
    args = ap.parse_args()
    root = Path(args.output)
    gtdb_info = prepare_gtdb(
        Path(args.gtdb_list), root / "gtdb", args.gtdb_queries, args.seed
    )
    hrom_info = prepare_hrom(
        Path(args.hrom_metadata),
        Path(args.hrom_taxonomy),
        Path(args.hrom_tgt_dir),
        Path(args.hrom_fasta_dir),
        root / "hrom",
        args.hrom_queries,
        args.seed,
    )
    manifest = {"gtdb": gtdb_info, "hrom": hrom_info, "seed": args.seed}
    with open(root / "prepare_manifest.json", "w") as out:
        json.dump(manifest, out, indent=2)
        out.write("\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
