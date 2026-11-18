"""Precise RR typing of resolved long-read loci.

The cargo pattern alone cannot separate RR4 from RR7, and the left boundary
genes distinguish the deletion endpoints. This script adds structural markers
from the RR3 reference and reports the boundary genes alongside the final
RR call. Oversized loci are kept separate because they are not one resolved
module.
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
BOUNDARY_GENES = ["STM2757", "STM2758", "STM2759", "STM2760", "STM2761"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--locus-dir", required=True)
    ap.add_argument("--fine-probes", required=True)
    ap.add_argument("--structure-probes", required=True)
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
        markers = rr_common.structure_markers(args.structure_probes, fasta)
        boundary = {gene: int(gene in genes) for gene in BOUNDARY_GENES}
        if locus_len > SIZE_LIMIT:
            call = "locus_oversized"
        else:
            call = rr_common.precise_classify(cargo, mer, intact, boundary,
                                              markers)
        info = meta.get(label, {})
        rows.append({
            "label": label,
            "dataset": info.get("dataset", "External"),
            "region": info.get("region", "Unknown"),
            "rr_call": call,
            "cargo_pattern": rr_common.cargo_string(mer, cargo),
            "locus_length": locus_len,
            "flj_intact": int(intact),
            "left_boundary": ",".join(g for g in BOUNDARY_GENES if boundary[g]),
            "IS26_count": markers["IS26_count"],
            "MAK1": int(markers["MAK1_78"] or markers["MAK1_75"]
                        or markers["MAK1_73"]),
            "meth": markers["meth"],
            "mer": mer, **cargo})
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    print(df["rr_call"].value_counts().to_string())
    print("\nRR_variant (novel candidates):")
    print(df[df["rr_call"] == "RR_variant"][["label", "region",
          "cargo_pattern", "left_boundary"]].to_string(index=False))


if __name__ == "__main__":
    main()
