#! /bin/bash

#SBATCH -J align_matrix
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=32
#SBATCH --nodelist=node[6]
#SBATCH --time=1-00:00:00
#SBATCH --output=align_matrix.o
#SBATCH --error=align_matrix.e

word_dir="/hwmaster/hujingling/SNAREseq2/scripts/011_RNA_quantification"
cd $workdir

MANIFEST_TSV="/hwmaster/hujingling/SNAREseq2/011_RNA_quantification/manifest_pre_exp2_random1k.tsv"
OUTPUT_NAME="/hwmaster/hujingling/SNAREseq2/011_RNA_quantification/STARsolo_output/pre_exp2/random1k"

singularity exec -B /hwmaster:/hwmaster /hwmaster/hujingling/singularity_images/starsolo_2.7.11b.sif \
  STAR \
  --readFilesManifest "${MANIFEST_TSV}" \
  --genomeDir /hwmaster/hujingling/SNAREseq2/reference_geno/STAR_index/genome_index \
  --outFileNamePrefix "${OUTPUT_NAME}" \
  --soloFeatures GeneFull_Ex50pAS   \
  --runThreadN 16 \
  --soloType SmartSeq \
  --outSAMtype BAM Unsorted \
  --outBAMcompression -1 \
  --soloUMIdedup NoDedup \
  --soloStrand Unstranded \
  --limitOutSJcollapsed 500000 \
  --soloCellFilter None \
  --outReadsUnmapped Fastx \
  --readFilesCommand zcat
 # --soloUMIdedup Exact \

#singularity exec --bind /hwmaster:/hwmaster /hwmaster/wenjichen/HumanAS_meta/singularity_container/starsolo:2.7.11b.sif \
#  STAR \
#  --readFilesIn /hwmaster/share/ZN_heart_plate/data/test_RNA/output_test_R1.fastq /hwmaster/share/ZN_heart_plate/data/test_RNA/output_test_R2.fastq \
#  --genomeDir /hwmaster/share/ZN_heart_plate/data/reference/refdata-gex-GRCh38-2024-A-STATsolo \
#  --outFileNamePrefix /hwmaster/share/ZN_heart_plate/output/01_RNA_quantification/PE \
#  --soloFeatures GeneFull_Ex50pAS \
#  --runThreadN 16 \
#  --soloType CB_samTagOut \
#  --soloCBposition 0_0_0_-1 \
#  --soloUMIposition 0_0_0_-1 \
#  --outSAMtype BAM Unsorted \
#  --outBAMcompression -1 \
#  --soloStrand Unstranded \
#  --limitOutSJcollapsed 500000 \
#  --outReadsUnmapped Fastx \
#  --outSAMattributes NH HI nM AS CR UR CB GX GN \
#  --soloCBmatchWLtype Exact \
#  --soloCBwhitelist None \
#  --soloBarcodeReadLength 0
