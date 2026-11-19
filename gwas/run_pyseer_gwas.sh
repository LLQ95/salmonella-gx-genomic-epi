#!/usr/bin/env bash
# Fluoroquinolone GWAS with pyseer, corrected for population structure.
# All paths can be overridden through environment variables.
set -euo pipefail

WORKDIR="${WORKDIR:-/data/lilanqi/salmonella/joint_gx_sea_tree}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${OUTDIR:-${WORKDIR}/gwas_output}"
VCF="${VCF:-${WORKDIR}/joint.core.vcf}"
ALIGNMENT="${ALIGNMENT:-${WORKDIR}/joint.core.aln}"
ANNOTATION_DIR="${ANNOTATION_DIR:-/data/lilanqi/salmonella/gx_mono_typhi_1874/seq_mono_typhi}"
POINTFINDER="${POINTFINDER:-${ANNOTATION_DIR}/pointfinder.tab}"
RESFINDER="${RESFINDER:-${ANNOTATION_DIR}/resfinder.tab}"
METADATA="${METADATA:-${WORKDIR}/unified_metadata.json}"
TRAIT_TABLE="${TRAIT_TABLE:-${WORKDIR}/trait_table.tsv}"
MDS_TYPE="${MDS_TYPE:-classic}"
MAX_DIMENSIONS="${MAX_DIMENSIONS:-10}"
PHENOTYPES="${OUTDIR}/gwas_phenotypes.tsv"
SIMILARITY="${OUTDIR}/core_similarity.tsv"

mkdir -p "${OUTDIR}"

echo "=== dependency check ==="
command -v python >/dev/null
command -v pyseer >/dev/null || {
    echo "pyseer not found. Install it with: pip install pyseer"
    exit 1
}

echo "=== prepare phenotypes ==="
python "${SCRIPT_DIR}/prepare_gwas_phenotypes.py" \
    --pointfinder "${POINTFINDER}" \
    --resfinder "${RESFINDER}" \
    --metadata "${METADATA}" \
    --trait-table "${TRAIT_TABLE}" \
    --out "${PHENOTYPES}"

echo "=== prepare similarity matrix ==="
if [ ! -f "${SIMILARITY}" ]; then
    python "${SCRIPT_DIR}/prepare_similarity.py" \
        --alignment "${ALIGNMENT}" \
        --out "${SIMILARITY}"
fi

echo "=== run pyseer for each phenotype ==="
for phenotype in qrdr gyrA parC qnrS1 oqxAB; do
    echo "--- ${phenotype} ---"
    pyseer \
        --vcf "${VCF}" \
        --phenotypes "${PHENOTYPES}" \
        --phenotype-column "${phenotype}" \
        --similarity "${SIMILARITY}" \
        --mds "${MDS_TYPE}" \
        --max-dimensions "${MAX_DIMENSIONS}" \
        --output-patterns "${OUTDIR}/${phenotype}_patterns.txt" \
        > "${OUTDIR}/${phenotype}_gwas.tsv"
done

echo "=== GWAS completed ==="
echo "outputs are under ${OUTDIR}"
