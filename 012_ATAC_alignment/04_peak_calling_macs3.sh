#!/usr/bin/env bash
set -euo pipefail

BAM=/hwmaster/hujingling/SNAREseq2/012_ATAC_alignment/minimap2_output/pre_exp2_AC_contami/output_ATACseq_bo_trimmed_aligned_sorted_mapq20_UBCB.bam
NAME="output_ATACseq_bo_trimmed_aligned_sorted_mapq20_UBCB_q0.1"
OUTDIR=/hwmaster/hujingling/SNAREseq2/012_ATAC_alignment/called_peaks/pre_exp2_AC_contami

[[ -f "$BAM" ]] || { echo "ERROR: BAM not found: $BAM"; exit 1; }

mkdir -p "$OUTDIR"

echo "MACS3 version:"
macs3 --version

macs3 callpeak \
  -t "$BAM" \
  -f BAM \
  -g hs \
  -n "$NAME" \
  --outdir "$OUTDIR" \
  --nomodel \
  -q 0.1
  -B \
  --shift -100 \
  --extsize 200 \

echo "Done. Peak file:"
echo "$OUTDIR/${NAME}_peaks.narrowPeak"
