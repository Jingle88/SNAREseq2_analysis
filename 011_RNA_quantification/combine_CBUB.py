#!/usr/bin/env python3
import pysam, sys

in_bam = sys.argv[1]
out_bam = sys.argv[2]

bam_in  = pysam.AlignmentFile(in_bam,  "rb")
bam_out = pysam.AlignmentFile(out_bam, "wb", template=bam_in)

for read in bam_in:
    try:
        cb = read.get_tag("CB")
        ub = read.get_tag("UB")
        combined = cb + ub  # concatenate
        read.set_tag("UB", combined, value_type="Z")  # overwrite UB
    except KeyError:
        # skip reads missing CB or UB
        pass
    bam_out.write(read)

bam_in.close()
bam_out.close()
print(" Done! Combined CB+UB into one 34-bp UB tag.")
