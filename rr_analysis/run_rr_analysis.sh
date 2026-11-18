#!/usr/bin/env bash
# Resistance region (RR) and novel RR analysis across Guangxi, adjacent
# provinces, and Southeast Asia. Override any path through the environment.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

WORKDIR="${WORKDIR:-/data/lilanqi/salmonella/joint_gx_sea_tree}"
LONGREAD_DIR="${LONGREAD_DIR:-${WORKDIR}/longreads}"
FASTA_DIR="${FASTA_DIR:-${WORKDIR}/input_fasta}"
TRAIT_TABLE="${TRAIT_TABLE:-${WORKDIR}/trait_table.tsv}"

PROBES="${SCRIPT_DIR}/probes"
OUT_BASE="${OUT_BASE:-${WORKDIR}/rr_novel_analysis}"
LOCUS_DIR="${OUT_BASE}/longread_loci"
mkdir -p "${OUT_BASE}"

echo "=== dependency check ==="
command -v blastn >/dev/null || { echo "blastn not found in PATH"; exit 1; }
python - <<'PY'
import importlib.util
for name in ("Bio", "pandas"):
    if not importlib.util.find_spec(name):
        raise SystemExit(f"missing python module: {name}")
print("python dependencies OK")
PY

if [ ! -f "${TRAIT_TABLE}" ]; then
    echo "trait table not found: ${TRAIT_TABLE}"
    exit 1
fi
echo "WORKDIR=${WORKDIR}"
echo "long reads: ${LONGREAD_DIR}"
echo "short-read assemblies: ${FASTA_DIR}"
echo "trait table: ${TRAIT_TABLE}"

echo
echo "=== Reference: RR3 GenBank and structural probes ==="
if [ ! -s "${PROBES}/KR856283.gb" ]; then
  python "${SCRIPT_DIR}/download_references.py" --out-dir "${PROBES}"
fi
if [ ! -s "${PROBES}/structure_probes.fasta" ]; then
  python "${SCRIPT_DIR}/prepare_structure_probes.py" \
    --reference "${PROBES}/KR856283.gb" \
    --out "${PROBES}/structure_probes.fasta"
fi

echo
echo "=== Step 1: extract continuous long-read loci ==="
python "${SCRIPT_DIR}/step1_extract_longread_loci.py" \
    --longread-dir "${LONGREAD_DIR}" \
    --anchor-probes "${PROBES}/robust_anchor_probes.fasta" \
    --out-dir "${LOCUS_DIR}"

echo
echo "=== Step 2: type long-read loci and confirm novel RR ==="
python "${SCRIPT_DIR}/step2_type_longread_loci.py" \
    --locus-dir "${LOCUS_DIR}" \
    --fine-probes "${PROBES}/fine_probes.fasta" \
    --resfinder "${PROBES}/all.fsa" \
    --mer-probe "${PROBES}/mer_probe.fasta" \
    --metadata "${TRAIT_TABLE}" \
    --out "${OUT_BASE}/longread_rr_calls.csv"

echo
echo "=== Step 5: precise structural typing of long-read loci ==="
python "${SCRIPT_DIR}/step5_precise_typing.py" \
    --locus-dir "${LOCUS_DIR}" \
    --fine-probes "${PROBES}/fine_probes.fasta" \
    --structure-probes "${PROBES}/structure_probes.fasta" \
    --resfinder "${PROBES}/all.fsa" \
    --mer-probe "${PROBES}/mer_probe.fasta" \
    --metadata "${TRAIT_TABLE}" \
    --out "${OUT_BASE}/longread_precise_calls.csv"

python "${SCRIPT_DIR}/apply_exclusion.py" \
    --input "${OUT_BASE}/longread_precise_calls.csv" \
    --exclude "${SCRIPT_DIR}/exclude_external_ndm.txt" \
    --out "${OUT_BASE}/longread_precise_calls_filtered.csv"

echo
echo "=== Step 3: scan all short-read assemblies ==="
python "${SCRIPT_DIR}/step3_scan_all_assemblies.py" \
    --fasta-dir "${FASTA_DIR}" \
    --fine-probes "${PROBES}/fine_probes.fasta" \
    --resfinder "${PROBES}/all.fsa" \
    --mer-probe "${PROBES}/mer_probe.fasta" \
    --metadata "${TRAIT_TABLE}" \
    --workers "${RR_WORKERS:-8}" \
    --out "${OUT_BASE}/shortread_rr_scan.csv"

echo
echo "=== Step 4: summarize by region and match novel RR ==="
python "${SCRIPT_DIR}/step4_summarize_by_region.py" \
    --longread-calls "${OUT_BASE}/longread_rr_calls.csv" \
    --shortread-scan "${OUT_BASE}/shortread_rr_scan.csv" \
    --out-dir "${OUT_BASE}/region_summary"

echo
echo "all outputs written under ${OUT_BASE}"
