#!/usr/bin/env python3
"""Compute a pairwise core-genome similarity matrix for pyseer.

Similarity is the fraction of aligned positions that are identical between
two samples. Gaps and ambiguous bases are excluded from the comparison. The
matrix is written with a header row and row labels in the format expected by
pyseer --similarity.
"""
import argparse
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np


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
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    names, matrix = read_alignment(args.alignment)
    n = len(names)
    valid = matrix > 64
    similarity = np.eye(n, dtype=np.float64)

    block = 128
    for start in range(0, n, block):
        end = min(start + block, n)
        block_valid = valid[start:end][:, None, :] & valid[None, :, :]
        block_equal = (matrix[start:end][:, None, :] ==
                       matrix[None, :, :]) & block_valid
        shared = block_valid.sum(axis=2)
        identical = block_equal.sum(axis=2)
        with np.errstate(invalid="ignore"):
            value = np.where(shared > 0, identical / shared, 0.0)
        similarity[start:end, :] = value

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("\t" + "\t".join(names) + "\n")
        for i, name in enumerate(names):
            row = "\t".join("{:.6f}".format(v) for v in similarity[i])
            handle.write(name + "\t" + row + "\n")
    print("wrote similarity matrix", n, "by", n)


if __name__ == "__main__":
    main()
