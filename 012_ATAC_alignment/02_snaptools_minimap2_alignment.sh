#!/bin/bash

# 使用说明函数
usage() {
    echo "Usage: $0 --input-dir <dir> --sample <sample_name> --output-dir <dir> [options]"
    echo ""
    echo "Required options:"
    echo "  --input-dir       Input directory containing fastq files"
    echo "  --sample          Sample name prefix"
    echo "  --output-dir      Output directory for BAM files"
    echo ""
    echo "Optional options:"
    echo "  --reference       Path to reference genome fasta (default: GRCh38-2024-A)"
    echo "  --threads         Number of threads (default: 16)"
    echo "  --aligner-path    Path to minimap2 binary (default: auto-detect)"
    echo "  --tmp-folder      Temporary folder (default: ./)"
    echo ""
    echo "Example:"
    echo "  $0 --input-dir ./processed --sample ATAC_SNARE2 --output-dir ./aligned"
    echo "  $0 --input-dir ./processed --sample ATAC_SNARE2 --output-dir ./aligned --threads 32"
    exit 1
}

# 默认参数
REFERENCE="/hwmaster/share/ZN_heart_plate/data/reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
THREADS=16
ALIGNER_PATH="/hwmaster/hujingling/miniforge3/envs/jupyter_env/bin/"
TMP_FOLDER="./"

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
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --aligner-path)
            ALIGNER_PATH="$2"
            shift 2
            ;;
        --tmp-folder)
            TMP_FOLDER="$2"
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

# 检查参考基因组是否存在
if [[ ! -f "$REFERENCE" ]]; then
    echo "Error: Reference genome not found: $REFERENCE"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 创建临时目录
mkdir -p "$TMP_FOLDER"

# 定义可能的文件扩展名模式
EXTENSIONS=("fastq.gz" "fq.gz")

# 函数：查找文件
find_file() {
    local sample=$1
    local read=$2
    local input_dir=$3
    
    for ext in "${EXTENSIONS[@]}"; do
        # 尝试多种命名模式
        for pattern in \
            "${sample}_CB_UB_${read}.${ext}" \
            "${sample}_${read}.${ext}" \
            "${sample}.${read}.${ext}" \
            "${sample}_${read:0:1}.${ext}"; do
            found=$(find "$input_dir" -maxdepth 1 -name "$pattern" | head -n 1)
            if [[ -n "$found" ]]; then
                echo "$found"
                return 0
            fi
        done
    done
    return 1
}

# 查找 R1 文件
R1_FILE=$(find_file "$SAMPLE" "R1" "$INPUT_DIR")

if [[ -z "$R1_FILE" ]]; then
    echo "Error: R1 file not found for sample: $SAMPLE"
    echo "Searched in: $INPUT_DIR"
    exit 1
fi

# 查找 R2 文件
R2_FILE=$(find_file "$SAMPLE" "R2" "$INPUT_DIR")

if [[ -z "$R2_FILE" ]]; then
    echo "Error: R2 file not found for sample: $SAMPLE"
    echo "Searched in: $INPUT_DIR"
    exit 1
fi

# 定义输出文件
OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE}_aligned.bam"

# 显示配置信息
echo "========================================="
echo "SnapTools Minimap2 Alignment"
echo "========================================="
echo "Sample:           $SAMPLE"
echo "Input directory:  $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Input files:"
echo "  R1: $R1_FILE"
echo "  R2: $R2_FILE"
echo ""
echo "Output BAM: $OUTPUT_BAM"
echo ""
echo "Parameters:"
echo "  Reference:    $REFERENCE"
echo "  Threads:      $THREADS"
echo "  Aligner path: $ALIGNER_PATH"
echo "  Tmp folder:   $TMP_FOLDER"
echo ""
echo "========================================="
echo "Starting alignment..."
echo "========================================="

# 运行 snaptools align-paired-end
snaptools align-paired-end \
  --input-reference="$REFERENCE" \
  --input-fastq1="$R1_FILE" \
  --input-fastq2="$R2_FILE" \
  --output-bam="$OUTPUT_BAM" \
  --aligner=minimap2 \
  --path-to-aligner="$ALIGNER_PATH" \
  --read-fastq-command=zcat \
  --min-cov 0 \
  --num-threads "$THREADS" \
  --if-sort True \
  --tmp-folder "$TMP_FOLDER" \
  --overwrite True

# 检查执行结果
if [[ $? -eq 0 ]]; then
    echo ""
    echo "========================================="
    echo "Alignment completed successfully!"
    echo "========================================="
    echo "Output BAM file:"
    ls -lh "$OUTPUT_BAM"
    
    # 如果 BAM index 存在，也显示
    if [[ -f "${OUTPUT_BAM}.bai" ]]; then
        echo ""
        echo "BAM index:"
        ls -lh "${OUTPUT_BAM}.bai"
    fi
    
    # 显示基本统计信息
    if command -v samtools &> /dev/null; then
        echo ""
        echo "BAM statistics:"
        samtools flagstat "$OUTPUT_BAM"
    fi
else
    echo ""
    echo "========================================="
    echo "Error: Alignment failed"
    echo "========================================="
    exit 1
fi
