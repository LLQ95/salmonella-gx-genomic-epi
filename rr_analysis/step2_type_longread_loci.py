"""Type each resolved long-read locus and flag novel RR candidates.

A locus above the size limit uses distant fallback anchors and is not treated
as one resolved locus. The RR5 pattern is recognized directly so that strains
carrying it are not mislabeled as novel.
"""
import argparse
import glob
import os
import sys
import pandas as pd
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rr_common

SIZE_LIMIT = 100000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--locus-dir", required=True)
    ap.add_argument("--fine-probes", required=True)
    ap.add_argument("--resfinder", required=True)
    ap.add_argument("--mer-probe", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    metadata = pd.read_csv(args.metadata, sep="\t", dtype={"node": str})
    meta = metadata.set_index("node").to_dict("index")
    rows = []
    for fasta in sorted(glob.glob(os.path.join(args.locus_dir, "*_locus.fasta"))):
        label = os.path.basename(fasta).replace("_locus.fasta", "")
        record = next(SeqIO.parse(fasta, "fasta"))
        locus_len = len(record.seq)
        genes = rr_common.present_genes(args.fine_probes, fasta)
        intact = "fljB" in genes and "hin" in genes
        cargo = rr_common.cargo_families(args.resfinder, fasta)
        mer = rr_common.mer_present(args.mer_probe, fasta)
        if locus_len > SIZE_LIMIT:
            call = "locus_oversized"
        else:
            call = rr_common.classify(cargo, mer, intact)
        info = meta.get(label, {})
        rows.append({
            "label": label,
            "dataset": info.get("dataset", "External"),
            "region": info.get("region", "Unknown"),
            "rr_call": call,
            "cargo_pattern": rr_common.cargo_string(mer, cargo),
            "locus_length": locus_len,
            "flj_intact": int(intact),
            "mer": mer, **cargo})
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    print(df["rr_call"].value_counts().to_string())
    print("\nRR_variant (novel candidates):")
    print(df[df["rr_call"] == "RR_variant"][["label", "region",
          "cargo_pattern"]].to_string(index=False))


if __name__ == "__main__":
    main()
