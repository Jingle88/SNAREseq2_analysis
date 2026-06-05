#!/usr/bin/env python3

import gzip
import sys
from tqdm import tqdm

fq = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else "fastq_structure_check.tsv"

bc_list = [
    "AACGTGAC", "AACTGTGA", "AAGCCTAG", "AATGCATG",
    "ACATCACG", "ACCTGGTG", "ACGTACTA", "ACGTTGAC"
]

polyT = "T" * 15
acr1_linker = "ACGTACTGCAGACTATGTCTACAG" 
redesignedAdapter2 = "CGGTAGCCGTCGCTCAGATGTGTATAAGAGACA"

max_mismatch = 1


def mismatch_le1(a, b):
    mismatches = 0
    for x, y in zip(a, b):
        if x != y:
            mismatches += 1
            if mismatches > 1:
                return False
    return True


def match_any(seq, targets):
    for t in targets:
        if mismatch_le1(seq, t):
            return True
    return False


total = 0
bc1_n = bc2_n = bc3_n = all_bc_n = 0
polyT_n = acr1_linker_n = RA2_n = all_bc_polyT_n = all_bc_ac_n = 0

with gzip.open(fq, "rt") as f:
    pbar = tqdm(desc="Checking reads", unit=" reads")

    while True:
        name = f.readline()
        if not name:
            break

        seq = f.readline().strip()
        f.readline()
        f.readline()

        total += 1
        pbar.update(1)

        bc3 = seq[10:18]
        bc2 = seq[48:56]
        bc1 = seq[86:94]
        dT = seq[94:109]
        acr1l = seq[94:118]
        ra2 = seq[118:151]

        bc1_ok = match_any(bc1, bc_list)
        bc2_ok = match_any(bc2, bc_list)
        bc3_ok = match_any(bc3, bc_list)

        all_bc_ok = bc1_ok and bc2_ok and bc3_ok
        polyT_ok = mismatch_le1(dT, polyT)
        acr1_link_ok = mismatch_le1(acr1l, acr1_linker)
        ra2_ok = mismatch_le1(ra2, redesignedAdapter2)

        bc1_n += bc1_ok
        bc2_n += bc2_ok
        bc3_n += bc3_ok
        all_bc_n += all_bc_ok
        polyT_n += polyT_ok
        acr1_linker_n += acr1_link_ok
        RA2_n += ra2_ok
        all_bc_polyT_n += all_bc_ok and polyT_ok
        all_bc_ac_n += all_bc_ok and acr1_link_ok and ra2_ok

    pbar.close()


def percent(n):
    return 0 if total == 0 else n / total * 100


rows = [
    ("total_reads_checked", total),
    ("BC1_match_86_94", bc1_n),
    ("BC2_match_48_56", bc2_n),
    ("BC3_match_10_18", bc3_n),
    ("all_BC_match", all_bc_n),
    ("polyT_match_94_109", polyT_n),
    ("acr1_linker_match_94_118", acr1_linker_n),
    ("ra2_match_118_151", RA2_n),
    ("all_BC_and_polyT", all_bc_polyT_n),
    ("all_BC_and_AC_DNA", all_bc_ac_n),
]

print("category\tcount\tpercent")
for name, count in rows:
    print(f"{name}\t{count}\t{percent(count):.2f}")

with open(out, "w") as o:
    o.write("category\tcount\tpercent\n")
    for name, count in rows:
        o.write(f"{name}\t{count}\t{percent(count):.2f}\n")

print(f"\nSaved to: {out}")
