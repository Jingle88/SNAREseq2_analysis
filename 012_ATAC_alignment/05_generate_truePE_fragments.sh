#!/bin/bash
set -euo pipefail

BAM="/hwmaster/hujingling/SNAREseq2/012_ATAC_alignment/minimap2_output/pre_exp1/output_ATAC_umi_trimmed_fixed_aligned.CBUB_dup_mapq20.bam"
OUT="PE_fragments.tsv.gz"

echo "[1] Filtering valid paired reads + MAPQ>=20"

samtools view -b \
  -f 2 \
  -q 20 \
  "$BAM" > filtered.bam

echo "[2] Sorting by name (required for bedtools bamtobed)"

samtools sort -n -o filtered.nameSorted.bam filtered.bam

echo "[3] Converting BAMPE -> BEDPE"

bedtools bamtobed -bedpe -i filtered.nameSorted.bam > fragments.bedpe

echo "[4] Converting BEDPE -> fragment format"

awk 'BEGIN{OFS="\t"} {
    if ($1==$4) {
        print $1, $2, $6, $7, 1
    }
}' fragments.bedpe | bgzip > PE_fragments.tsv.gz

echo "[DONE] fragments.tsv.gz generated"
