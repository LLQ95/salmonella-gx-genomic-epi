#!/usr/bin/env python3
"""Build pyseer phenotype tables for the fluoroquinolone GWAS.

Chromosomal QRDR mutations are read from PointFinder and plasmid-mediated
qnrS and oqxAB carriage is read from ResFinder. The output is a TSV with one
row per sample and binary phenotype columns, suitable for pyseer.
"""
import argparse
import json
import os
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
import pandas as pd


def strain_from_pointfinder(path):
    parts = path.split("/")
    return parts[1] if len(parts) > 1 else path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pointfinder", required=True)
    ap.add_argument("--resfinder", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.metadata, encoding="utf-8") as handle:
        metadata = json.load(handle)
    samples = sorted(metadata.keys())
    table = pd.DataFrame({"samples": samples})
    for column in ("qrdr", "gyrA", "parC", "qnrS1", "oqxAB"):
        table[column] = 0

    point = pd.read_csv(args.pointfinder, sep="\t", header=None, dtype=str)
    point["strain"] = point[0].map(strain_from_pointfinder)
    for _, row in point.iterrows():
        mutation = row[1]
        strain = row["strain"]
        index = table.index[table["samples"] == strain]
        if len(index) == 0:
            continue
        i = index[0]
        if mutation.startswith("gyrA") or mutation.startswith("gyrB"):
            table.at[i, "gyrA"] = 1
            table.at[i, "qrdr"] = 1
        if mutation.startswith("parC") or mutation.startswith("parE"):
            table.at[i, "parC"] = 1
            table.at[i, "qrdr"] = 1

    res = pd.read_csv(args.resfinder, sep="\t", dtype=str)
    res["strain"] = res["#FILE"].map(
        lambda p: os.path.basename(p).replace(".fasta", ""))
    for _, row in res.iterrows():
        gene = row["GENE"]
        strain = row["strain"]
        index = table.index[table["samples"] == strain]
        if len(index) == 0:
            continue
        i = index[0]
        if gene.startswith("qnrS"):
            table.at[i, "qnrS1"] = 1
        if gene.startswith("OqxA") or gene.startswith("OqxB"):
            table.at[i, "oqxAB"] = 1

    table.to_csv(args.out, sep="\t", index=False)
    for column in ("qrdr", "gyrA", "parC", "qnrS1", "oqxAB"):
        print(column, int(table[column].sum()))
    print("samples", len(table))


if __name__ == "__main__":
    main()
