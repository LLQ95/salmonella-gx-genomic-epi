"""Scan every short-read assembly for boundary genes and resistance cargo.

Short-read assemblies are fragmented and do not yield one continuous locus,
so the analysis is restricted to gene-content signatures. Each fasta may hold
several contigs and is passed directly to BLASTn, which searches all records.
Strains are processed in parallel and results are joined to the trait table.
"""
import argparse
import glob
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rr_common

BOUNDARY_GENES = ["STM2744", "STM2745", "STM2757", "STM2758", "STM2759",
                  "STM2760", "STM2761", "STM2765", "STM2768", "STM2769",
                  "fljA", "fljB", "hin", "iroB", "iroC"]


def scan_one(arguments):
    fasta, fine_probes, resfinder, mer_probe, meta_lookup = arguments
    label = os.path.basename(fasta).replace(".fasta", "")
    genes = rr_common.present_genes(fine_probes, fasta)
    cargo = rr_common.cargo_families(resfinder, fasta)
    mer = rr_common.mer_present(mer_probe, fasta)
    info = meta_lookup.get(label, {})
    row = {"label": label,
           "dataset": info.get("dataset", "Unknown"),
           "region": info.get("region", "Unknown"),
           "cargo_pattern": rr_common.cargo_string(mer, cargo),
           "mer": mer, **cargo}
    for gene in BOUNDARY_GENES:
        row[gene] = int(gene in genes)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta-dir", required=True)
    ap.add_argument("--fine-probes", required=True)
    ap.add_argument("--resfinder", required=True)
    ap.add_argument("--mer-probe", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    metadata = pd.read_csv(args.metadata, sep="\t", dtype={"node": str})
    meta_lookup = metadata.set_index("node").to_dict("index")
    files = sorted(glob.glob(os.path.join(args.fasta_dir, "*.fasta")))
    tasks = [(f, args.fine_probes, args.resfinder, args.mer_probe,
              meta_lookup) for f in files]

    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(scan_one, task) for task in tasks]
        for done, future in enumerate(as_completed(futures), start=1):
            rows.append(future.result())
            if done % 250 == 0:
                print("processed", done, "of", len(files))
    df = pd.DataFrame(rows).sort_values("label")
    df.to_csv(args.out, index=False)
    print("wrote", args.out, len(df), "rows")
    print(df["dataset"].value_counts().to_string())


if __name__ == "__main__":
    main()
