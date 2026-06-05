#!/bin/bash

# 用法检查
if [[ $# -ne 2 ]]; then
    echo "用法: $0 输入目录 输出manifest文件"
    echo "示例: $0 /path/to/fastq manifest_smartseq.tsv"
    exit 1
fi

input_dir="$1"
output_manifest="$2"

# 检查输入目录是否存在
if [[ ! -d "$input_dir" ]]; then
    echo "错误：输入目录不存在：$input_dir"
    exit 1
fi

# 获取输出文件的绝对路径
output_manifest=$(realpath "$output_manifest")

# 创建输出目录（如果不存在）
output_dir=$(dirname "$output_manifest")
mkdir -p "$output_dir"

# 清空或创建manifest文件
> "$output_manifest"

# 计数器
count=0

# 进入输入目录
cd "$input_dir" || exit 1

# 遍历所有R1文件
for r1_file in *_R1_*fixed.part_*.{fastq,fq}{,.gz}; do
    if [[ -f "$r1_file" ]]; then
        # 提取CB（去掉R1后缀）
        if [[ "$r1_file" =~ fixed\.part_([ACGTN]+)\.(fastq|fq)(\.gz)?$ ]]; then
	#if [[ "$r1_file" =~ (.*)_R1_fixed\.(fastq|fq)(\.gz)?$ ]]; then
            cb="${BASH_REMATCH[1]}"
        else
            continue
        fi

        # 对应的R2文件
        r2_file="${r1_file/R1/R2}"

        # 检查R2文件是否存在
        if [[ -f "$r2_file" ]]; then
            # 获取完整路径
            full_path_r1="$(pwd)/$r1_file"
            full_path_r2="$(pwd)/$r2_file"

            # 写入manifest文件（格式：R2路径 R1路径 样本名）
            echo -e "$full_path_r2\t$full_path_r1\t$cb" >> "$output_manifest"
            ((count++))
        else
            echo "警告：找不到对应的R2文件：$r2_file"
        fi
    fi
done

echo "完成：生成了 $count 个样本的manifest文件：$output_manifest"
