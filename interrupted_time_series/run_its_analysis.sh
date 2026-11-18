#!/usr/bin/env bash
# Interrupted time series analysis of resistance determinants.
# Paths can be overridden through environment variables.
set -euo pipefail

WORKDIR="${WORKDIR:-/data/lilanqi/salmonella/joint_gx_sea_tree}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${OUTDIR:-${WORKDIR}/its_output}"
METADATA="${METADATA:-${WORKDIR}/unified_metadata.json}"
PHENOTYPES="${PHENOTYPES:-${WORKDIR}/gwas_output/gwas_phenotypes.tsv}"
START_YEAR="${START_YEAR:-2014}"
END_YEAR="${END_YEAR:-2024}"
INTERVENTION_YEAR="${INTERVENTION_YEAR:-2017}"
RSCRIPT="${RSCRIPT:-Rscript}"

mkdir -p "${OUTDIR}"
ANNUAL="${OUTDIR}/annual_series.tsv"
COEFFICIENTS="${OUTDIR}/its_coefficients.tsv"
FIGURE="${OUTDIR}/its_figure.png"

echo "=== dependency check ==="
command -v python >/dev/null
command -v "${RSCRIPT}" >/dev/null

echo "=== prepare annual series ==="
python "${SCRIPT_DIR}/prepare_annual_series.py" \
    --phenotypes "${PHENOTYPES}" \
    --metadata "${METADATA}" \
    --start-year "${START_YEAR}" \
    --end-year "${END_YEAR}" \
    --out "${ANNUAL}"

echo "=== fit interrupted time series ==="
"${RSCRIPT}" "${SCRIPT_DIR}/fit_interrupted_series.R" \
    "${ANNUAL}" "${COEFFICIENTS}" "${FIGURE}" \
    "${INTERVENTION_YEAR}"

echo "=== ITS completed ==="
echo "outputs are under ${OUTDIR}"
