#! /bin/bash

#SBATCH -J umi_dedup
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=32
#SBATCH --nodelist=node[5]
#SBATCH --time=1-00:00:00
#SBATCH --output=umi_dedup.o
#SBATCH --error=umi_dedup.e


#samtools sort -@ 8 -o /hwmaster/hujingling/SNAREseq2/011_RNA_quantification/STARsolo_output/pre_exp1/pre_exp1Aligned.out.CBUB.sorted.bam /hwmaster/hujingling/SNAREseq2/011_RNA_quantification/STARsolo_output/pre_exp1/pre_exp1Aligned.out.CBUB.bam
#samtools index /hwmaster/hujingling/SNAREseq2/011_RNA_quantification/STARsolo_output/pre_exp1/pre_exp1Aligned.out.CBUB.sorted.bam

#umi_tools dedup \
#  -I /hwmaster/hujingling/SNAREseq2/011_RNA_quantification/STARsolo_output/pre_exp1/pre_exp1Aligned.out.CBUB.sorted.bam \
#  --extract-umi-method=tag \
#  --umi-tag=UB --cell-tag=CB \
#  --paired \
#  --output-stats=dedup_sum.txt \
#  -S /hwmaster/hujingling/SNAREseq2/011_RNA_quantification/STARsolo_output/pre_exp1/pre_exp1Aligned.out.umitools.dedup.bam \
#  --log=dedup.log

rumi /hwmaster/hujingling/SNAREseq2/011_RNA_quantification/STARsolo_output/pre_exp1/pre_exp1Aligned.out.combined_CBUB.sorted.bam \
  --output /hwmaster/hujingling/SNAREseq2/011_RNA_quantification/STARsolo_output/pre_exp1/pre_exp1Aligned.out.rumi.dedup.bam \
  --umi_tag UB
