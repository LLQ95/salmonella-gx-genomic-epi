#!/usr/bin/env python3
"""Pairwise core-SNP distances between Guangxi and reference isolates.

The script compares every Guangxi isolate against the Southeast Asia,
adjacent-province and anchor isolates in the joint core alignment. It counts
close pairs at several SNP thresholds and writes both a JSON summary and a
CSV table of close pairs with reference metadata.
"""
import argparse
import json
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd


def read_alignment(path):
    names = []
    sequences = []
    current = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                    current = []
                names.append(line[1:].strip().split()[0])
            else:
                current.append(line.strip())
    if current:
        sequences.append("".join(current))
    matrix = np.array([[ord(c) for c in sequence] for sequence in sequences],
                      dtype=np.int16)
    return names, matrix


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--outprefix", required=True)
    ap.add_argument("--thresholds", default="5,12,25")
    args = ap.parse_args()

    names, matrix = read_alignment(args.alignment)
    print("alignment:", matrix.shape)
    with open(args.metadata, encoding="utf-8") as handle:
        metadata = json.load(handle)
    thresholds = [int(value) for value in args.thresholds.split(",")]

    gx_index = [i for i, name in enumerate(names)
                if metadata.get(name, {}).get("dataset") == "Guangxi"]
    ref_index = [i for i in range(len(names)) if i not in gx_index]
    ref_matrix = matrix[ref_index]
    ref_names = [names[i] for i in ref_index]

    pair_counts = {threshold: 0 for threshold in thresholds}
    gx_with_close = {threshold: set() for threshold in thresholds}
    rows = []

    for i in gx_index:
        query = matrix[i]
        valid = (query > 64) & (ref_matrix > 64)
        distances = np.count_nonzero(
            (ref_matrix != query) & valid, axis=1)
        for ref_position, distance in enumerate(distances):
            distance = int(distance)
            for threshold in thresholds:
                if distance <= threshold:
                    pair_counts[threshold] += 1
                    gx_with_close[threshold].add(names[i])
            if distance <= max(thresholds):
                ref_id = ref_names[ref_position]
                record = metadata.get(ref_id, {})
                rows.append({
                    "gx_id": names[i],
                    "ref_id": ref_id,
                    "snp_distance": distance,
                    "ref_dataset": record.get("dataset", ""),
                    "ref_country": record.get("country", ""),
                    "ref_region": record.get("region", ""),
                    "ref_serovar": record.get("serovar", ""),
                    "ref_st": record.get("st", ""),
                    "ref_date": record.get("date", ""),
                })

    close_pairs = pd.DataFrame(rows)
    close_pairs.to_csv(args.outprefix + "_close_pairs.csv", index=False)

    summary = {
        "total_guangxi": len(gx_index),
        "total_reference": len(ref_index),
        "close_pairs": pair_counts,
        "guangxi_isolates_with_close_ref": {
            str(threshold): len(gx_with_close[threshold])
            for threshold in thresholds},
    }
    with open(args.outprefix + "_summary.json", "w",
              encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))
    print("close pair rows:", len(close_pairs))


if __name__ == "__main__":
    main()
