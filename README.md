Prefiltering + demultiplexing
1. download the reference genome （.fasta + .gtf) from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/. index the genome using STARsolo
2. Check the barcode mismatch, AC R1 linker match, and poly T match with RNA_mismatch_check.py
3. pre-filter the raw fastq file with ATACseq_hjl_redesBC.py and RNAseq_hjl_redesBC.py seperately, removing all the 3-round barcode mismatch reads + umi duplicates
4. linke the UB and CB information in header with '_' by change_fq_header


scRNAseq analysis:
1. split the fq files according to 3-round barcode combination (single nucleus) - 01_split_fq.sh
2. generate manifest information - 02_generate_manifest.sh
3. Alignment to reference genome by STARsolo and output bam file - 03_align_matrix.sh
4. attribute the UB and CB information in bam - bam_CB_UB.py 
5. Generate cell x expression matrix
6. perform downstream QC with Seurat, see tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial



scATACseq analysis:
1. trim adapters in reads - 01_trim_adapter.sh (for single-end mode, using 01_trim_adapter_R1.sh)
2. Alignment to reference genome by minimap2 and output bam file - 02_snaptools_minimap2_alignment.sh
3. attribute the UB and CB information in bam - 03_rename_bam.py
4. perform QC (alignment filtering according to MAPQ..) ( depends on the data quality) on bam file, index the final bam file
5. perform peak calling using macs3 - 04_peak_calling_macs3.sh
6. generate the fragments file -05_generate_truePE_fragments.sh (if single-end mode, using 05_generate _pseudo_fragments.sh)
7. perform downstream QC with signac - ATAC_QC_visualization.R, see tutorial https://stuartlab.org/signac/articles/pbmc_vignette.html
