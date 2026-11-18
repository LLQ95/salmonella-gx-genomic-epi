#!/usr/bin/env python3
"""Aggregate binary resistance determinants into annual prevalence series.

The GWAS phenotype table is joined with the sample year from the metadata and
summarized by year. For each determinant the output reports the number of
positive samples, the total number of samples and the prevalence.
"""
import argparse
import json
import sys
sys.stdout.reconfigure(encoding="utf-8")
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phenotypes", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--start-year", type=int, default=2014)
    ap.add_argument("--end-year", type=int, default=2024)
    args = ap.parse_args()

    phenotypes = pd.read_csv(args.phenotypes, sep="\t", dtype=str)
    with open(args.metadata, encoding="utf-8") as handle:
        metadata = json.load(handle)
    years = {key: record.get("date") for key, record in metadata.items()}
    phenotypes["year"] = phenotypes["samples"].map(years)
    phenotypes = phenotypes[phenotypes["year"].notna()]
    phenotypes["year"] = phenotypes["year"].astype(int)
    phenotypes = phenotypes[
        (phenotypes["year"] >= args.start_year) &
        (phenotypes["year"] <= args.end_year)]

    determinants = ["qrdr", "gyrA", "parC", "qnrS1", "oqxAB"]
    for column in determinants:
        phenotypes[column] = phenotypes[column].astype(int)

    rows = []
    for year, group in phenotypes.groupby("year"):
        row = {"year": int(year), "total": len(group)}
        for determinant in determinants:
            positive = int(group[determinant].sum())
            row[determinant + "_positive"] = positive
            row[determinant + "_prevalence"] = positive / len(group)
        rows.append(row)

    result = pd.DataFrame(rows).sort_values("year")
    result.to_csv(args.out, sep="\t", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
