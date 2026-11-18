"""Filter a strain-level table against an exclusion list.

Labels listed in the exclusion file are removed. Lines starting with # and
blank lines are ignored. The filtered table is written to a new path so the
original output is left unchanged.
"""
import argparse
import sys
import pandas as pd

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def read_exclusion(path):
    labels = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                labels.append(line)
    return set(labels)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--exclude", required=True)
    ap.add_argument("--label-column", default="label")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    excluded = read_exclusion(args.exclude)
    df = pd.read_csv(args.input, dtype=str)
    before = len(df)
    df = df[~df[args.label_column].isin(excluded)]
    df.to_csv(args.out, index=False)
    print("before", before, "after", len(df), "removed",
          before - len(df))


if __name__ == "__main__":
    main()
