#!/usr/bin/env python3
import sys, re, pysam

in_bam  = sys.argv[1]
out_bam = sys.argv[2]

# allow underscore OR whitespace between CB and UB in the read name
pattern = re.compile(r".*CB:Z:([ACGTN]+)[_\s]UB:Z:([ACGTN]+)")

# open input (try BAM first, fallback to SAM)
try:
    bam_in = pysam.AlignmentFile(in_bam, "rb")
    in_mode = "bam"
except ValueError:
    bam_in = pysam.AlignmentFile(in_bam, "r")  # SAM text
    in_mode = "sam"

# always write BAM (binary) with the same header
bam_out = pysam.AlignmentFile(out_bam, "wb", template=bam_in)

n_total = n_tagged = 0
for read in bam_in:
    n_total += 1
    # if CB/UB already present, keep them; otherwise try to extract from read name
    has_cb = read.has_tag("CB")
    has_ub = read.has_tag("UB")
    if not (has_cb and has_ub):
        m = pattern.match(read.query_name)
        if m:
            cb, ub = m.groups()
            read.set_tag("CB", cb, value_type="Z")
            read.set_tag("UB", ub, value_type="Z")
            n_tagged += 1
    bam_out.write(read)

bam_in.close()
bam_out.close()

sys.stderr.write(
    f"[INFO] Input mode: {in_mode}; reads processed: {n_total}; newly tagged: {n_tagged}\n"
)

