"""Build structural probes from the annotated RR3 reference KR856283.

The reference carries the genes used to separate RR types that share the same
resistance cargo, including IS26 (counted by copy number), the MAK1 cluster,
the meth gene, and the IS1 and IS10 transposase genes.
"""
import argparse
import os
import sys
from Bio import SeqIO

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# One representative gene is kept for each structural marker. Repeated IS26
# copies are identical, so a single IS26 probe is used and its copy number is
# counted from the number of BLAST hits.
MARKERS = {
    "IS26": "tnpAIS26",
    "IS1": "insB_IS1",
    "IS10": "tnpAIS10",
    "MAK1_78": "MAK1.78c",
    "MAK1_75": "MAK1.75",
    "MAK1_73": "MAK1.73",
    "meth": "methdelta",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    record = next(SeqIO.parse(args.reference, "genbank"))
    found = {}
    for feature in record.features:
        if feature.type != "gene":
            continue
        label = (feature.qualifiers.get("gene", [""])[0]
                 or feature.qualifiers.get("locus_tag", [""])[0])
        for marker, source in MARKERS.items():
            if label == source and marker not in found:
                found[marker] = feature.location.extract(record).seq
    missing = [m for m in MARKERS if m not in found]
    if missing:
        raise RuntimeError("missing markers in reference: " + ",".join(missing))
    with open(args.out, "w", encoding="utf-8") as handle:
        for marker, sequence in found.items():
            handle.write(f">{marker}\n")
            for index in range(0, len(sequence), 70):
                handle.write(str(sequence[index:index + 70]) + "\n")
    print("wrote", args.out, len(found), "markers")


if __name__ == "__main__":
    main()
