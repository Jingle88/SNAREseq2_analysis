#!/bin/bash
[[ $# -lt 2 ]] && { echo "用法: $0 输入目录 输出目录 [线程数]"; exit 1; }

for r1 in "$1"/*R1*fixed.{fq,fastq}{,.gz}; do
    [[ -f "$r1" ]] || continue
    r2="${r1/R1/R2}"
    [[ -f "$r2" ]] || continue
    seqkit split -j ${3:-8} --by-id --id-regexp ".*CB:Z:([ATCGN]+).*" --out-dir "$2" "$r1"
    seqkit split -j ${3:-8} --by-id --id-regexp ".*CB:Z:([ATCGN]+).*" --out-dir "$2" "$r2"
done

# 重命名输出文件，去除不必要的前缀
for f in "$2"/*.fastq{,.gz}; do
    [[ -f "$f" ]] || continue
    base=$(basename "$f")
    
    # 匹配: prefix_R1.part_BARCODE.fastq(.gz)
    if [[ "$base" =~ ^(.+)_(R[12])\.part_([ATCGN]+)\.(fastq\.gz|fastq)$ ]]; then
        barcode="${BASH_REMATCH[3]}"
        read_type="${BASH_REMATCH[2]}"
        ext="${BASH_REMATCH[4]}"
        mv "$f" "$2/${barcode}_${read_type}.${ext}"
    fi
done
