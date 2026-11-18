"""Shared helpers for resistance region (RR) analysis.

The module wraps BLASTn calls and maps ResFinder hits to a small set of
resistance gene families so that different allele numbers of the same family
are counted only once.
"""
import os
import sys
import subprocess

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

CARGO_FAMILIES = ["blaTEM", "strA", "strB", "sul2", "tetB"]
FAMILY_PREFIX = {
    "blaTEM": "blaTEM-1B",
    "strA": "aph(3'')-Ib",
    "strB": "aph(6)-Id",
    "sul2": "sul2",
    "tetB": "tet(B)",
}
# Order used for the compact binary cargo string.
CARGO_ORDER = ["mer", "blaTEM", "strA", "strB", "sul2", "tetB"]


def blastn(query, subject, identity, qcov, out_fields, extra=None):
    cmd = ["blastn", "-query", query, "-subject", subject, "-task", "blastn",
           "-outfmt", "6 " + " ".join(out_fields),
           "-perc_identity", str(identity),
           "-qcov_hsp_perc", str(qcov)]
    if extra:
        cmd.extend(extra)
    run = subprocess.run(cmd, capture_output=True, text=True)
    if run.returncode != 0:
        raise RuntimeError(run.stderr)
    return run.stdout.splitlines()


def present_genes(probe_fasta, subject_fasta, identity=85, qcov=70):
    """Return the set of probe genes detected in a subject fasta."""
    lines = blastn(probe_fasta, subject_fasta, identity, qcov, ["qseqid"])
    return {line.split("\t")[0] for line in lines}


def cargo_families(resfinder_fasta, subject_fasta, identity=90, qcov=80):
    """Return a 0/1 dict for each resistance gene family.

    ResFinder genes are used as the query. A long structural region is the
    subject, which avoids the very low query coverage that occurs when a
    multi-gene region is treated as the query.
    """
    lines = blastn(resfinder_fasta, subject_fasta, identity, qcov,
                   ["qseqid", "pident", "length", "qlen"])
    best = {}
    for line in lines:
        qseqid, pident, length, qlen = line.split("\t")
        coverage = int(length) / int(qlen)
        for family, prefix in FAMILY_PREFIX.items():
            if qseqid.startswith(prefix):
                if family not in best or coverage > best[family]:
                    best[family] = coverage
    return {family: int(family in best) for family in CARGO_FAMILIES}


def mer_present(mer_probe, subject_fasta, identity=80, qcov=70):
    lines = blastn(mer_probe, subject_fasta, identity, qcov,
                   ["pident", "length", "qlen"])
    for line in lines:
        pident, length, qlen = line.split("\t")
        if int(length) / int(qlen) >= 0.7 and float(pident) >= 80:
            return 1
    return 0


def cargo_string(mer, cargo):
    values = [str(mer)] + [str(cargo[f]) for f in CARGO_FAMILIES]
    return "".join(values)


def structure_markers(structure_probe, subject_fasta, identity=85, qcov=70):
    """Return structural marker presence and the IS26 copy number.

    Each BLAST line is one alignment, so the number of IS26 lines gives the
    copy number rather than a simple presence flag.
    """
    lines = blastn(structure_probe, subject_fasta, identity, qcov, ["qseqid"])
    counts = {}
    for line in lines:
        marker = line.split("\t")[0]
        counts[marker] = counts.get(marker, 0) + 1
    result = {marker: int(marker in counts) for marker in
              ("IS1", "IS10", "MAK1_78", "MAK1_75", "MAK1_73", "meth")}
    result["IS26_count"] = counts.get("IS26", 0)
    return result


def classify(cargo, mer, flj_intact):
    has_any_assu = any(cargo[k] for k in ("blaTEM", "strA", "strB", "sul2"))
    if flj_intact:
        return "intact_fljAB_with_cargo" if has_any_assu else "intact_fljAB"
    pattern = cargo_string(mer, cargo)
    # Known cargo patterns read from Figure 1A of Wang et al. 2023, in the
    # order mer, blaTEM, strA, strB, sul2, tetB.
    if pattern == "111111":
        return "RR3"
    if pattern == "101111":
        return "RR5"
    if pattern == "001110":
        return "RR8"
    if pattern == "100001":
        return "RR6"
    if pattern == "000001":
        # RR4 and RR7 share this cargo and differ only in structure.
        return "RR4_RR7_like"
    if pattern == "000000":
        return "deletion_no_cargo"
    return "RR_variant"


def precise_classify(cargo, mer, flj_intact, boundary, markers):
    """Type a resolved locus using cargo and structural evidence together."""
    has_any_assu = any(cargo[k] for k in ("blaTEM", "strA", "strB", "sul2"))
    if flj_intact:
        return "intact_fljAB_with_cargo" if has_any_assu else "intact_fljAB"
    pattern = cargo_string(mer, cargo)
    has_mak = markers["MAK1_78"] or markers["MAK1_75"] or markers["MAK1_73"]
    if pattern == "111111":
        return "RR3"
    if pattern == "101111":
        return "RR5"
    if pattern == "001110":
        return "RR8"
    if pattern == "100001":
        return "RR6"
    if pattern == "000001":
        # RR4 keeps the MAK1/meth region; RR7 has lost it.
        return "RR4" if (has_mak or markers["meth"]) else "RR7"
    if pattern == "000000":
        return "deletion_no_cargo"
    return "RR_variant"
