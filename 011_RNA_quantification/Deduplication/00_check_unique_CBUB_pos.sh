#! /bin/bash

#SBATCH -J umi_pos_dedup
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=32
#SBATCH --nodelist=node[2]
#SBATCH --time=1-00:00:00
#SBATCH --output=UMI_pos.o
#SBATCH --error=UMI_pos.e

python3 /hwmaster/hujingling/SNAREseq2/scripts/011_RNA_quantification/Deduplication/00_check_unique_CBUB_pos.py
