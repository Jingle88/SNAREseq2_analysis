import gzip
import io
import sys

# Define configuration constants
ACR1_LINKER = "ACGTACTGCAGACTATGTCTACAG"
REDESIGNED_ADAPTER2 = "CGGTAGCCGTCGCTCAGATGTGTATAAGAGACA"

# Expected coordinates (0-indexed for Python slicing)
# 94:118 captures length 24 | 118:151 captures length 33
ACR1_START, ACR1_END = 94, 118
RA2_START, RA2_END = 118, 151


def mismatch_le1(a, b):
    """Returns True if string 'a' and 'b' have <= 1 mismatch."""
    # Quick length check safety guard
    if len(a) != len(b):
        return False

    mismatches = 0
    for x, y in zip(a, b):
        if x != y:
            mismatches += 1
            if mismatches > 1:
                return False
    return True


def filter_paired_fastq(r1_in_path, r2_in_path, r1_out_path, r2_out_path):
    """Filters paired-end FASTQ.gz files based on R2 sequence specific regions."""
    print("Starting filtering process... This may take a while depending on file size.")

    count_total = 0
    count_passed = 0

    # Open all 4 files simultaneously using gzip with text mode wrappers
    with gzip.open(r1_in_path, "rb") as f1_in, gzip.open(r2_in_path, "rb") as f2_in, \
         gzip.open(r1_out_path, "wb") as f1_out, gzip.open(r2_out_path, "wb") as f2_out:

        # Wrap in BufferedReader/TextIOWrapper for fast line-by-line string reading
        r1_reader = io.TextIOWrapper(io.BufferedReader(f1_in), encoding="utf-8")
        r2_reader = io.TextIOWrapper(io.BufferedReader(f2_in), encoding="utf-8")
        
        r1_writer = io.TextIOWrapper(io.BufferedWriter(f1_out), encoding="utf-8")
        r2_writer = io.TextIOWrapper(io.BufferedWriter(f2_out), encoding="utf-8")

        while True:
            # Read 4 lines for R1
            r1_id = r1_reader.readline()
            if not r1_id:
                break  # End of file reached
            r1_seq = r1_reader.readline()
            r1_plus = r1_reader.readline()
            r1_qual = r1_reader.readline()

            # Read 4 lines for R2
            r2_id = r2_reader.readline()
            r2_seq = r2_reader.readline()
            r2_plus = r2_reader.readline()
            r2_qual = r2_reader.readline()

            count_total += 1

            # Extract the exact coordinates from R2 sequence
            # .strip() removes the trailing newline character '\n'
            clean_r2_seq = r2_seq.strip()
            
            # Extract target substrings
            read_acr1 = clean_r2_seq[ACR1_START:ACR1_END]
            read_ra2 = clean_r2_seq[RA2_START:RA2_END]

            # Check if BOTH regions satisfy the maximum 1-mismatch criteria
            if mismatch_le1(read_acr1, ACR1_LINKER) and mismatch_le1(read_ra2, REDESIGNED_ADAPTER2):
                count_passed += 1
                
                # Write matching records to output files
                r1_writer.write(f"{r1_id}{r1_seq}{r1_plus}{r1_qual}")
                r2_writer.write(f"{r2_id}{r2_seq}{r2_plus}{r2_qual}")

            if count_total % 1000000 == 0:
                print(f"Processed {count_total} reads... Found {count_passed} valid pairs.")

    print("\n--- Processing complete ---")
    print(f"Total reads evaluated: {count_total}")
    # Display percentage with up to 2 decimal places
    print(f"Reads passing filter:  {count_passed} ({count_passed / count_total * 100:.2f}%)")


if __name__ == "__main__":
    INPUT_R1 = "/hwmaster/hujingling/SNAREseq2/fastq_raw/MEIJI_seq/20260518/RNA.R1.raw.fastq.gz"
    INPUT_R2 = "/hwmaster/hujingling/SNAREseq2/fastq_raw/MEIJI_seq/20260518/RNA.R2.raw.fastq.gz"
    
    OUTPUT_R1 = "/hwmaster/hujingling/SNAREseq2/fastq_raw/MEIJI_seq/20260518/AC_contami/AC_contami_R1_raw.fastq.gz"
    OUTPUT_R2 = "/hwmaster/hujingling/SNAREseq2/fastq_raw/MEIJI_seq/20260518/AC_contami/AC_contami_R2_raw.fastq.gz"

    filter_paired_fastq(INPUT_R1, INPUT_R2, OUTPUT_R1, OUTPUT_R2)
