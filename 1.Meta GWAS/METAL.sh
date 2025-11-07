#=============================#
# METAL.sh
#=============================#
#！/bin/bash
set -euo pipefail

##############################################################
#  Flexible METAL runner (Enhanced)
#  - Supports unlimited GWAS inputs (--gwas, --gwas1, --gwas2 ...)
#  - Supports --output-prefix <prefix>
#  - Supports --maf <float>  (default 0.01)
#
#  Example:
#    bash METAL.sh \
#        --output-prefix meta_test \
#        --maf 0.05 \
#        --gwas1 EUR.gz --gwas2 AFR.gz --gwas3 AMR.gz
##############################################################
OUT_PREFIX="meta_output"
GWAS_FILES=()
MAF=0.01

# === getopt ===
while [[ $# -gt 0 ]]; do
    case "$1" in 
        --output-prefix|--outname)
            OUT_PREFIX="$2"
            shift 2
            ;;
        --maf)
            MAF="$2"
            shift 2
            ;;
        --gwas|--gwas[0-9]*)
            GWAS_FILES+=("$2")
            shift 2
            ;;
        -h|--help)
            echo "Usage:"
            echo " bash $0 --outname output_prefix --maf 0.01 --gwas1 file1 --gwas2 file2 [--gwas3 file3 ...]"
            echo ""
            echo "example:"
            echo " bash $0 --outname ALL_MVP \\"
            echo "    --maf 0.01 \\"
            echo "    --gwas1 /path/EUR.gz --gwas2 /path/AFR.gz --gwas3 /path/EAS.gz"
            exit 0
            ;;
        *)
            echo "❌ Unknown parameter: $1"
            exit 1
            ;;
    esac
done

# === params check ===
if [[ ${#GWAS_FILES[@]} -eq 0 ]]; then 
    echo "❌ Error: No GWAS files provided. Use --gwas or --gwas1, --gwas2, etc. to specify input files."
    exit 1
fi

MAXFREQ=$(awk -v maf="$MAF" 'BEGIN { printf "%.6f", 1 - maf }')
echo "🧩 Outfile prefix: $OUT_PREFIX"
echo "🧬 Input ${#GWAS_FILES[#]} GWAS files: "
for f in "${GWAS_FILES[@]}"; do
    echo "  - $f"
done
echo "📈 Filter: ${MAF} < freq < ${MAXFREQ}"

# === check gwas files ===
for f in "${GWAS_FILES[@]}"; do
    if [[  ! -f "$f" ]]; then
        echo "⚠️  Warning: GWAS file not exist -> $f"
    fi
done

# === Generate metal config file ===
CONF_FILE="${OUT_PREFIX}_metal.conf"


cat > "$CONF_FILE" << EOF
# ==================================
# METAL Configuration Auto-Generated
# ==================================
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
VERBOSE OFF

# Input columns:
CHROMOSOME CHR
MARKER SNP
ALLELE A1 A2
EFFECT beta
STDERR SE
FREQ freq
PVALUE p
WEIGHT N
GENOMICCONTROL OFF

# Filters:
ADDFILTER freq > ${MAF}
ADDFILTER freq < ${MAXFREQ}
TRACKPOSITIONS ON 
CHROMOSOME CHR
POSITION POS

# additional options
EFFECT_PRINT_PRECISION 12
STDERR_PRINT_PRECISION 12
CUSTOMVARIABLE TotalSampleSize  
LABEL TotalSampleSize as N
EOF

# === add GWAS files to metal config ===
for f in "${GWAS_FILES[@]}"; do
    echo "PROCESS ${f}" >> "$CONF_FILE"
done


cat >> "$CONF_FILE" << EOF

OUTFILE ${OUT_PREFIX} .tbl
ANALYZE HETEROGENEITY
QUIT
EOF

# === run metal ===
metal="/public/home/shilulu/software/METAL/build/bin/metal"
echo "🚀 Config file generated: $CONF_FILE"
echo "🔧 Submitting METAL job ..."
qsubshcom "$metal $CONF_FILE" 1 100G METAL 2:00:00 ""

echo "✅ METAL analysis submitted."
echo "📂 Output prefix: ${OUT_PREFIX}"