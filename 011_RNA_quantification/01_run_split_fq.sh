#! /bin/bash

#SBATCH -J split_fq
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=32
#SBATCH --nodelist=node[6]
#SBATCH --time=1-00:00:00
#SBATCH --output=split_fq.o
#SBATCH --error=split_fq.e

word_dir="/hwmaster/hujingling/SNAREseq2/scripts/011_RNA_quantification"
cd $workdir

pre_filtered_fq="/hwmaster/hujingling/SNAREseq2/011_RNA_quantification/pre_filtered_fq/pre_exp2/random1k"
splited_fq="/hwmaster/hujingling/SNAREseq2/011_RNA_quantification/splited_fq/pre_exp2/random1k"


bash /hwmaster/hujingling/SNAREseq2/scripts/011_RNA_quantification/01_split_fq.sh "${pre_filtered_fq}" "${splited_fq}"

echo "Finishing the spliting of fastq file accoridng to cell barcode combination"
