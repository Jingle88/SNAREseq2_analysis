#! /bin/bash

#SBATCH -J RNA_BC_filter
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=32
#SBATCH --nodelist=node[6]
#SBATCH --time=1-00:00:00
#SBATCH --output=RNA_BC_filtering.o
#SBATCH --error=RNA_BC_filtering.e

# Load the environment
source /hwmaster/hujingling/miniforge3/etc/profile.d/conda.sh
conda activate jupyter_env

python3 /hwmaster/hujingling/SNAREseq2/scripts/00_demultiplexing/RNAseq_hjl_redesBC.py
