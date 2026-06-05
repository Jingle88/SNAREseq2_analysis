#!/bin/bash

# 使用说明函数
usage() {
    echo "Usage: $0 --input-dir <dir> --sample <sample_name> --output-dir <dir> [--threads <n>]"
    echo "Example: $0 --input-dir ./data --sample ATAC_SNARE2 --output-dir ./trimmed_data --threads 8"
    exit 1
}

# 默认线程数
THREADS=4

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --sample)
            SAMPLE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown parameter: $1"
            usage
            ;;
    esac
done

# 检查必需参数
if [[ -z "$INPUT_DIR" || -z "$SAMPLE" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required parameters"
    usage
fi

# 检查输入目录是否存在
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 定义可能的文件扩展名模式
EXTENSIONS=("fastq.gz" "fq.gz" "fastq" "fq")

# 查找 R1 文件
R1_FILE=""
for ext in "${EXTENSIONS[@]}"; do
    # 尝试多种命名模式
    for pattern in "${SAMPLE}_R1_fixed.${ext}"; do
        found=$(find "$INPUT_DIR" -maxdepth 1 -name "$pattern" | head -n 1)
        if [[ -n "$found" ]]; then
            R1_FILE="$found"
            break 2
        fi
    done
done

# 查找 R2 文件
R2_FILE=""
for ext in "${EXTENSIONS[@]}"; do
    # 尝试多种命名模式
    for pattern in "${SAMPLE}_R2_fixed.${ext}"; do
        found=$(find "$INPUT_DIR" -maxdepth 1 -name "$pattern" | head -n 1)
        if [[ -n "$found" ]]; then
            R2_FILE="$found"
            break 2
        fi
    done
done

# 检查是否找到文件
if [[ -z "$R1_FILE" ]]; then
    echo "Error: R1 file not found for sample: $SAMPLE"
    echo "Searched in: $INPUT_DIR"
    exit 1
fi

if [[ -z "$R2_FILE" ]]; then
    echo "Error: R2 file not found for sample: $SAMPLE"
    echo "Searched in: $INPUT_DIR"
    exit 1
fi

# 定义输出文件
OUTPUT_R1="${OUTPUT_DIR}/${SAMPLE}_trimmed_R1.fastq.gz"
OUTPUT_R2="${OUTPUT_DIR}/${SAMPLE}_trimmed_R2.fastq.gz"

# 显示找到的文件
echo "Found R1 file: $R1_FILE"
echo "Found R2 file: $R2_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo ""

# 运行 fastp 进行接头去除
echo "Running fastp for adapter trimming..."
fastp \
  -i "$R1_FILE" \
  -I "$R2_FILE" \
  -o "$OUTPUT_R1" \
  -O "$OUTPUT_R2" \
  --thread "$THREADS" \
  --detect_adapter_for_pe \
  --compression 6

# 检查执行结果
if [[ $? -eq 0 ]]; then
    echo ""
    echo "Success! Trimmed files generated:"
    echo "  R1: $OUTPUT_R1"
    echo "  R2: $OUTPUT_R2"
    
    # 显示文件大小
    echo ""
    echo "File sizes:"
    ls -lh "$OUTPUT_R1" "$OUTPUT_R2"
else
    echo "Error: fastp failed"
    exit 1
fi
