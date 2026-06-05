#!/usr/bin/env python3
import pysam
from collections import defaultdict, Counter

in_bam  = "/hwmaster/hujingling/SNAREseq2/011_RNA_quantification/STARsolo_output/pre_exp1/pre_exp1Aligned.out.CBUB.sorted.bam"
out_bam = "/hwmaster/hujingling/SNAREseq2/011_RNA_quantification/STARsolo_output/pre_exp1/pre_exp1Aligned.out.CBUB_dedup.bam"

bam_in  = pysam.AlignmentFile(in_bam, "rb")
bam_out = pysam.AlignmentFile(out_bam, "wb", template=bam_in)

def five_prime(read):
    if read.is_reverse:
        return read.reference_end
    return read.reference_start

print("PASS 1 scanning...")

# -------- PASS 1 ---------
# We only store ONE representative qname per molecule key → low memory
cbub_rep = {}  # key = (cb, ub, strand, ref, pos_bin)  → representative qname

total_all  = 0
total_cbub = 0
primary    = 0

ub_counts       = Counter()
cbub_counts     = Counter()  # (cb, ub, strand)
cbub_pos_counts = Counter()  # (cb, ub, strand, ref, pos_bin)

for read in bam_in:
    total_all += 1

    # require CB & UB
    if not (read.has_tag("CB") and read.has_tag("UB")):
        continue
    if read.is_unmapped:
        continue

    total_cbub += 1

    # only dedup/use primary alignments
    if read.is_secondary or read.is_supplementary:
        continue

    primary += 1

    cb = read.get_tag("CB")
    ub = read.get_tag("UB")
    strand = "-" if read.is_reverse else "+"
    ref = bam_in.get_reference_name(read.reference_id)
    pos = five_prime(read)
    pos_bin = (pos // 5) * 5
    key = (cb, ub, strand, ref, pos_bin)

    # statistics
    ub_counts[ub] += 1
    cbub_counts[(cb, ub, strand)] += 1
    cbub_pos_counts[key] += 1

    # keep only the FIRST encountered primary read as the representative
    if key not in cbub_rep:
        cbub_rep[key] = read.query_name

bam_in.close()

# -------- STATISTICS ----------
dup_umi_groups      = sum(1 for v in ub_counts.values() if v > 1)
dup_cbub_groups     = sum(1 for v in cbub_counts.values() if v > 1)
dup_cbubpos_groups  = sum(1 for v in cbub_pos_counts.values() if v > 1)

print("\n=== Deduplication Stats ===")
print(f"Total alignments in BAM:                 {total_all}")
print(f"Alignments with CB+UB tags:             {total_cbub}")
print(f"Primary reads used for dedup:           {primary}")
print(f"Total unique UMIs:                      {len(ub_counts)}")
print(f"Total unique CB+UB+strand:              {len(cbub_counts)}")
print(f"Total unique CB+UB+strand+pos groups:   {len(cbub_rep)}")
print(f"Total unique CB+UB+strand+pos groups:   {len(cbub_pos_counts)}")
print(f"Duplicate UMI groups:                   {dup_umi_groups}")
print(f"Duplicate CB+UB+strand groups:          {dup_cbub_groups}")
print(f"Duplicate molecule position groups:     {dup_cbubpos_groups}")

# -------- PASS 2: write one representative read per molecule key ---------
print("PASS 2 writing output...")

bam_in2 = pysam.AlignmentFile(in_bam, "rb")
written = 0

# Extract a set of representative qnames (very small memory)
representative_qnames = set(cbub_rep.values())

for read in bam_in2:

    # We write alignments inc. (primary, secondary, supplementary) -> one primary alignment per molecule
    if read.query_name in representative_qnames:
        bam_out.write(read)
        written += 1


bam_in2.close()
bam_out.close()

print(f"\nWrote {written} deduplicated molecules (reads) to {out_bam}")
