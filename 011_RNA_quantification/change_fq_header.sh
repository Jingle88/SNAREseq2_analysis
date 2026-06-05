#!/bin/bash

FASTQ_DIR=$1

# --- Step 1. Fix headers in all fastq.gz files --------------------
echo "Rewriting FASTQ headers..."
mkdir -p "${FASTQ_DIR}/fixed"

for f in ${FASTQ_DIR}/*.fastq.gz; do
  base=$(basename "$f" .fastq.gz)
  echo "  Processing $base..."
  gunzip -c "$f" | \
  awk 'NR%4==1 {
          # split header into tokens, remove leading "@"
          readid=$1; sub(/^@/,"",readid);
          cb=""; ub="";
          for(i=2;i<=NF;i++){
              if($i ~ /^CB:Z:/) cb=substr($i,6);
              else if($i ~ /^UB:Z:/) ub=substr($i,6);
          }
          # rebuild header as one token: READID_CB:Z:..._UB:Z:...
          print "@" readid "_CB:Z:" cb "_UB:Z:" ub;
          next
       }
       {print}' | gzip > "${FASTQ_DIR}/fixed/${base}_fixed.fastq.gz"
done

echo "Header rewriting complete. Fixed files in ${FASTQ_DIR}/fixed" 
