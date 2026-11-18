"""Extract a single continuous fljAB locus from each long-read genome.

Anchor genes are searched as a left chain and a right chain. The first
detected anchor on each side defines the locus span. The span is built with
min/max coordinates because the linearization origin of a circular chromosome
can place either anchor at the lower coordinate.
"""
import argparse
import glob
import os
import sys
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rr_common

LEFT_CHAIN = ["STM2744", "STM2743", "STM2742", "STM2741", "STM2769"]
RIGHT_CHAIN = ["iroC", "iroD", "iroE", "iroN", "STM2778", "STM2779",
               "STM2780", "STM2786", "STM2787", "STM2788", "STM2789",
               "STM2797"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--longread-dir", required=True)
    ap.add_argument("--anchor-probes", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    summary = []
    for fasta in sorted(glob.glob(os.path.join(args.longread_dir, "*.fasta"))):
        label = os.path.basename(fasta).replace(".fasta", "")
        hits = rr_common.blastn(args.anchor_probes, fasta, 85, 70,
                               ["qseqid", "sstart", "send"])
        coords = {}
        for line in hits:
            qseqid, sstart, send = line.split("\t")
            coords[qseqid] = (int(sstart), int(send))
        left = next((x for x in LEFT_CHAIN if x in coords), None)
        right = next((x for x in RIGHT_CHAIN if x in coords), None)
        if not left or not right:
            summary.append((label, "no_anchor", 0))
            continue
        lo = min(min(coords[left]), min(coords[right]))
        hi = max(max(coords[left]), max(coords[right]))
        record = next(SeqIO.parse(fasta, "fasta"))
        sub = record.seq[lo - 1:hi]
        out = os.path.join(args.out_dir, label + "_locus.fasta")
        SeqIO.write(type(record)(sub, id=label, description="fljAB_locus"),
                    out, "fasta")
        summary.append((label, "ok", len(sub)))
    for row in summary:
        print(row)
    ok = [r for r in summary if r[1] == "ok"]
    print("extracted", len(ok), "of", len(summary), "long-read loci")


if __name__ == "__main__":
    main()
