#!/bin/bash

FASTQ_DIR=$1
FIXED_DIR=$2
THREADS=${3:-4}

mkdir -p "$FIXED_DIR"

echo "Rewriting FASTQ headers using multi-threading..."

for f in ${FASTQ_DIR}/*.fastq.gz; do
  (
    base=$(basename "$f" .fastq.gz)
    echo "  Processing $base ..."

    pigz -dc "$f" | \
    mawk 'NR%4==1 {
            readid=$1; sub(/^@/, "", readid);
            cb=""; ub="";
            for(i=2;i<=NF;i++){
                if($i ~ /^CB:Z:/) cb=substr($i,6);
                else if($i ~ /^UB:Z:/) ub=substr($i,6);
            }
            print "@" readid "_CB:Z:" cb "_UB:Z:" ub;
            next
         }
         {print}' \
    | pigz -p $THREADS > "${FIXED_DIR}/${base}_fixed.fastq.gz"

    echo "  ✔ Done $base"
  ) &
done

wait
echo "All FASTQ files processed."
