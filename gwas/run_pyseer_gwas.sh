#!/usr/bin/env bash
# Fluoroquinolone GWAS with pyseer, corrected for population structure.
# All paths can be overridden through environment variables.
set -euo pipefail

WORKDIR="${WORKDIR:-/data/lilanqi/salmonella/joint_gx_sea_tree}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${OUTDIR:-${WORKDIR}/gwas_output}"
VCF="${VCF:-${WORKDIR}/joint.core.vcf}"
ALIGNMENT="${ALIGNMENT:-${WORKDIR}/joint.core.aln}"
POINTFINDER="${POINTFINDER:-${WORKDIR}/guangxi_1874/pointfinder.tab}"
RESFINDER="${RESFINDER:-${WORKDIR}/guangxi_1874/resfinder.tab}"
METADATA="${METADATA:-${WORKDIR}/unified_metadata.json}"
TRAIT_TABLE="${TRAIT_TABLE:-${WORKDIR}/trait_table.tsv}"
MDS_COMPONENTS="${MDS_COMPONENTS:-10}"
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
        --mds "${MDS_COMPONENTS}" \
        --output-patterns "${OUTDIR}/${phenotype}_patterns.txt" \
        > "${OUTDIR}/${phenotype}_gwas.tsv"
done

echo "=== GWAS completed ==="
echo "outputs are under ${OUTDIR}"
