# for single-end ATACseq analysis
install.packages('tidyverse')
BiocManager::install(version = "3.22")
install.packages("ggplot2")
install.packages("data.table")
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
BiocManager::install(c(
  "GenomicRanges",
  "IRanges",
  "S4Vectors",
  "Rsamtools",
  "rtracklayer",
  "EnsDb.Hsapiens.v98",
  "biovizBase"
), update = TRUE, ask = FALSE)
# BiocManager::install("rtracklayer") # may not be a good idea to install it directly
library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(S4Vectors)
library(GenomeInfoDb)

install.packages("Signac")
library(Signac)
library(Seurat)

BiocManager::install("AnnotationHub")
library(AnnotationHub)



pre_exp2_AC = ""
pre_exp2_AC_contami = ""
pre_exp1_AC_PE = ""
pre_exp1_AC_SE = ""



# standard QC for the alignment files (for those mapped reads > 10,000)
# distribution of aligned reads across genome
idx <- fread(file.path(pre_exp2_AC_contami, "ACcontami_idxstat_mapq20.txt"), header = FALSE)
colnames(idx) <- c("chr", "length", "mapped", "unmapped")
total_mapped <- sum(idx[chr != "*", mapped])

idx_plot <- idx[chr != "*" & mapped > 10000]

ggplot(idx_plot, aes(x = reorder(chr, mapped), y = mapped)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(x = "Chromosome", y = "Mapped reads",
       title = "Mapped reads per chromosome") 
# chrM percentage = 0.002709326 (20755/7660577) for ATACseq
# chrM percentage = 0.000947465 (12285/12966172) for AC contami (mapq>20)


# standard QC for peaks called from macs3
# check the peak width distribution (insert size distribution is meaningless here because we don't have paired end reads anymore -> no true fragment -> no distribution of NFR, mononucleosome...)
peaks <- fread(file.path(pre_exp1_AC_PE, "output_ATACseq_bo_trimmed_aligned_sorted_mapq20_UBCB_q0.1_peaks.narrowPeak"), header = FALSE)
colnames(peaks)[1:3] <- c("chr", "start", "end")

peaks[, width := end - start]

ggplot(peaks, aes(x = width)) +
  geom_histogram(bins = 50) +
  theme_bw() +
  labs(x = "Peak width", y = "Number of peaks",
       title = "Peak width distribution")

# check the peak score distribution
colnames(peaks)[5] <- "score"

ggplot(peaks, aes(x = score)) +
  geom_histogram(bins = 50) +
  theme_bw() +
  labs(x = "MACS2 peak score", y = "Number of peaks",
       title = "Peak score distribution")


'''
01. Standard Pre-processing workflow 
'''

# load called peaks file
peaks <- import(file.path(pre_exp2_AC, "output_ATACseq_bo_trimmed_aligned_sorted_mapq20_UBCB_peaks.narrowPeak"))

# load the pseudo-fragment file
fragpath <- file.path(pre_exp2_AC, "pseudo_fragments_mapq20.tsv.gz")
total_counts <- CountFragments(fragpath) #get cell barcode combinations -> here 512
cutoff <- 100
barcodes <- total_counts[total_counts$frequency_count > cutoff,]$CB
frags <- CreateFragmentObject(path = fragpath, cells = barcodes)


# quantify fragment counts in each peak -> create peak x cell matrix
counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)

# create the chromatin assays
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  fragments = frags
)

chrom_assay_mincell10_minfea200 <- CreateChromatinAssay(
  counts = counts,
  fragments = frags,
  min.cells = 10,
  min.features = 200
)

# create the seurat object
obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

obj_mincell10_minfea200 <- CreateSeuratObject(
  counts = chrom_assay_mincell10_minfea200,
  assay = "peaks"
)

'''
> obj
An object of class Seurat 
3247 features across 512 samples within 1 assay 
Active assay: peaks (3247 features, 0 variable features)
 2 layers present: counts, data
 
 > obj_mincell10_minfea200
An object of class Seurat 
2425 features across 227 samples within 1 assay 
Active assay: peaks (2425 features, 0 variable features)
 2 layers present: counts, data
'''
# check distribution of the genomic ranges of called peaks
granges(obj)

# keep those peaks in standard chromosome (chr1-22+X+Y)
peaks.keep <- seqnames(granges(obj)) %in% standardChromosomes(granges(obj))
obj_standardChr <- obj[as.vector(peaks.keep), ]

# Add gene annotation
Sys.setenv(BIOCONDUCTOR_ONLINE_MIRROR="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
library(AnnotationHub)
ah <- AnnotationHub()

# Search for the Ensembl 98 EnsDb for Homo sapiens on AnnotationHub
query(ah, "EnsDb.Hsapiens.v98")
ensdb_v98 <- ah[["AH75011"]]

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(obj) <- annotations



'''
02. QC Metrics
'''

# compute nucleosome signal score per cell -> meaningless because here we use single-end read + pseudo-fragments
obj <- NucleosomeSignal(object = obj)
FragmentHistogram(object = obj)
                  
# compute TSS enrichment score per cell
obj <- TSSEnrichment(object =obj)
DensityScatter(obj, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

# add fraction of reads in peaks
library(readr)
frag_df <- read_tsv(fragpath, col_names = FALSE)
names(frag_df) <- c("chrom", "start", "end", "barcode", "readCount")
peaks <- granges(obj)
frag_gr <- GRanges(
  seqnames = frag_df$chrom,
  ranges = IRanges(frag_df$start, frag_df$end),
  barcode = frag_df$barcode
)
hits <- findOverlaps(frag_gr, peaks)
frag_df$in_peak <- FALSE
frag_df$in_peak[queryHits(hits)] <- TRUE
frip_per_cell <- tapply(
  frag_df$in_peak,
  frag_df$barcode,
  mean
)
obj$FRiP <- frip_per_cell[colnames(obj)]
FRiP_global <- length(unique(queryHits(hits))) / nrow(frag_df)
FRiP_global
obj$pct_reads_in_peaks <- obj$FRiP * 100

# add blacklist ratio
blacklist_regions <- ah[['AH107305']] # blacklist regions for hg38
obj$blacklist_ratio <- FractionCountsInRegion(
  object = obj, 
  assay = 'peaks',
  regions = blacklist_regions
)

VlnPlot(
  object = obj,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

VlnPlot(
  obj,
  features = c(
    "nCount_peaks",
    "nFeature_peaks"
  )
)
ggplot(
  data = obj@meta.data, # <=== THE FIX: Pass the metadata data frame, not the whole object
  aes(
    x = nCount_peaks,
    y = nFeature_peaks
  )
) +
  geom_point(size = 0.5, alpha = 0.5, color = "midnightblue") +
  theme_classic() +
  labs(
    title = "Library Complexity Profile",
    x = "Sequencing Depth (nCount_peaks)",
    y = "Unique Peaks Detected (nFeature_peaks)"
  )

# downstream dimension reduction and clutering ? -> after getting reliable QC metrics for those peak calling result!!!
