samtools view -F 4 /hwmaster/hujingling/SNAREseq2/012_ATAC_alignment/minimap2_output/pre_exp2/output_ATACseq_bo_trimmed_aligned_sorted_mapq20_UBCB.bam | \
awk '
BEGIN{OFS="\t"}

{
    chr=$3
    pos=$4

    cb=""

    for(i=12;i<=NF;i++)
    {
        if($i ~ /^CB:Z:/)
        {
            split($i,a,":")
            cb=a[3]
        }
    }

    start=pos-100
    if(start<0) start=0

    end=pos+100

    print chr,start,end,cb,1
}
' | \
sort -k1,1 -k2,2n > pseudo_fragments_mapq20.tsv
