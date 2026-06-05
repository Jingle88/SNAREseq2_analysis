#!/usr/bin/env python3

import gzip
from collections import defaultdict
from tqdm import tqdm


R1fadir = '/hwmaster/hujingling/SNAREseq2/fastq_raw/MEIJI_seq/20260518/AC.R1.raw_random1k.fastq.gz'
R2fqdir = '/hwmaster/hujingling/SNAREseq2/fastq_raw/MEIJI_seq/20260518/AC.R2.raw_random1k.fastq.gz'

bc1_list = ['AACGTGAC','AACTGTGA','AAGCCTAG','AATGCATG','ACATCACG','ACCTGGTG','ACGTACTA','ACGTTGAC']
bc2_list = ['AACGTGAC','AACTGTGA','AAGCCTAG','AATGCATG','ACATCACG','ACCTGGTG','ACGTACTA','ACGTTGAC']
bc3_list = ['AACGTGAC','AACTGTGA','AAGCCTAG','AATGCATG','ACATCACG','ACCTGGTG','ACGTACTA','ACGTTGAC']

# Convert lists to sets for faster lookup
bc1_set = set(bc1_list)
bc2_set = set(bc2_list)
bc3_set = set(bc3_list)
bias = 0
max_mismatch = 1   # allow 0, 1, or more mismatches

def hamming_dist(s1, s2):
    return sum(a != b for a, b in zip(s1, s2))

def match_barcode(obs_bc, bc_list, max_mismatch=0):
    """
    Return matched whitelist barcode if obs_bc is within max_mismatch.
    If multiple barcodes match equally well, return None to avoid ambiguity.
    """
    best_bc = None
    best_dist = max_mismatch + 1
    n_best = 0

    for ref_bc in bc_list:
        dist = hamming_dist(obs_bc, ref_bc)

        if dist < best_dist:
            best_dist = dist
            best_bc = ref_bc
            n_best = 1
        elif dist == best_dist:
            n_best += 1

    if best_dist <= max_mismatch and n_best == 1:
        return best_bc
    else:
        return None

def open_maybe_gz(path, mode="rt"):
    """
    Open plain text or gzipped FASTQ.
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def fastq_iter(path):
    """
    Streaming FASTQ parser.
    Yields: header_without_at, sequence, quality
    """
    with open_maybe_gz(path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break

            seq = f.readline()
            plus = f.readline()
            qual = f.readline()

            if not qual:
                raise ValueError(f"Incomplete FASTQ record detected in {path}")

            # Similar to pyfastx .name: remove @ and keep the read name
            name = header.rstrip("\n\r")[1:].split()[0]

            yield (
                name,
                seq.rstrip("\n\r"),
                qual.rstrip("\n\r")
            )


# For each CB, record UMIs already seen
seen_umis_by_bc = defaultdict(set)

# Same meaning as the old final outputs
n_pos_index = 0          # number of reads passing barcode filtering before UMI deduplication
bc_set = set()           # unique CBs among filtered reads
n_remain = 0             # number of reads retained after UMI deduplication

r1_iter = fastq_iter(R1fadir)
r2_iter = fastq_iter(R2fqdir)

with open('output_RNAseq_bo_R2.fastq', 'w', buffering=1024 * 1024) as out_r2, \
     open('output_RNAseq_bo_R1.fastq', 'w', buffering=1024 * 1024) as out_r1:

    for (r1_name, r1_seq, r1_qual), (r2_name, r2_seq, r2_qual) in tqdm(
        zip(r1_iter, r2_iter),
        desc="Filtering + deduplicating + writing"
    ):

        # Extract observed barcodes
        obs_bc1 = r2_seq[86 + bias:94 + bias]
        obs_bc2 = r2_seq[48 + bias:56 + bias]
        obs_bc3 = r2_seq[10 + bias:18 + bias]

        # Match to whitelist with mismatch allowance
        bc1 = match_barcode(obs_bc1, bc1_list, max_mismatch)
        if bc1 is None:
            continue

        bc2 = match_barcode(obs_bc2, bc2_list, max_mismatch)
        if bc2 is None:
            continue

        bc3 = match_barcode(obs_bc3, bc3_list, max_mismatch)
        if bc3 is None:
            continue
        umi = r2_seq[0:10]
        bc = bc1 + bc2 + bc3

        n_pos_index += 1
        bc_set.add(bc)

        # Same logic as:
        # for each CB, keep first occurrence of each UMI
        if umi in seen_umis_by_bc[bc]:
            continue

        seen_umis_by_bc[bc].add(umi)
        n_remain += 1

        # Write R2 output, same as original:
        # fq[pos_index[i]].seq[109:150]
        out_r2.write(f"@{r2_name + ' CB:Z:' + bc + ' UB:Z:' + umi}\n")
        out_r2.write(f"{r2_seq[109:150]}\n")
        out_r2.write("+\n")
        out_r2.write(f"{r2_qual[109:150]}\n")

        # Write R1 output, same as original:
        # fq[pos_index[i]].seq
        out_r1.write(f"@{r1_name + ' CB:Z:' + bc + ' UB:Z:' + umi}\n")
        out_r1.write(f"{r1_seq}\n")
        out_r1.write("+\n")
        out_r1.write(f"{r1_qual}\n")


print('pos_index')
print(n_pos_index)

print('bc')
print(len(bc_set))

print('remain')
print(n_remain)
