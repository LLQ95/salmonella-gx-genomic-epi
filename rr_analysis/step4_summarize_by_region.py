"""Summarize RR signatures by region and match novel RR patterns.

Long-read calls provide the structural evidence. Short-read strains that share
the same cargo pattern are reported as candidates rather than confirmed novel
RRs, because fragmented assemblies cannot establish the deletion boundaries.
"""
import argparse
import os
import sys
import pandas as pd

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--longread-calls", required=True)
    ap.add_argument("--shortread-scan", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    longread = pd.read_csv(args.longread_calls,
                           dtype={"label": str, "cargo_pattern": str})
    shortread = pd.read_csv(args.shortread_scan,
                            dtype={"label": str, "cargo_pattern": str})

    confirmed_novel = longread[longread["rr_call"] == "RR_variant"]
    novel_patterns = sorted(confirmed_novel["cargo_pattern"].unique())
    print("confirmed novel RR cargo patterns from long reads:", novel_patterns)

    # Region-level distribution of cargo patterns across all short reads.
    region_table = (shortread.groupby(["dataset", "cargo_pattern"])
                    .size().reset_index(name="n")
                    .sort_values(["dataset", "n"], ascending=[True, False]))
    region_path = os.path.join(args.out_dir, "region_cargo_patterns.csv")
    region_table.to_csv(region_path, index=False)

    # Short-read strains sharing a long-read-confirmed novel pattern.
    matches = shortread[shortread["cargo_pattern"].isin(novel_patterns)].copy()
    matches["evidence"] = "shortread_candidate_needs_longread"
    match_path = os.path.join(args.out_dir, "novel_rr_shortread_candidates.csv")
    matches.to_csv(match_path, index=False)

    print("\nshort-read candidates by dataset:")
    if len(matches):
        print(matches.groupby(["dataset", "cargo_pattern"]).size()
              .to_string())
    else:
        print("none")

    print("\nregion cargo pattern counts (top entries per dataset):")
    for dataset in region_table["dataset"].unique():
        sub = region_table[region_table["dataset"] == dataset].head(8)
        print("\n", dataset)
        print(sub.to_string(index=False))


if __name__ == "__main__":
    main()
