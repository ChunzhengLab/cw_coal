#!/bin/bash

# ==============================================================================
# HTCondor作业包装脚本
# ==============================================================================

# ----------------------------
# 1. 环境初始化
# ----------------------------
source /storage/fdunphome/wangchunzheng/miniconda3/etc/profile.d/conda.sh
conda activate cpp_dev

# ----------------------------
# 2. 环境变量配置
# ----------------------------
# 设置编译器和清理旧ROOT环境
export CC="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-cc"
export CXX="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++"
unset ROOTSYS

# 路径清理
export PATH=$(echo ":$PATH:" | sed -E 's#:/opt/root61404/bin:#:#g' | sed 's#^:##;s#:$##')
export LD_LIBRARY_PATH=$(echo ":$LD_LIBRARY_PATH:" | sed -E 's#:/opt/root61404/lib:#:#g' | sed 's#^:##;s#:$##')

# 设置优先级路径
export PATH="$CONDA_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="/storage/fdunphome/wangchunzheng/CWCoalProject/lib64:$LD_LIBRARY_PATH"

# ----------------------------
# 3. 参数解析
# ----------------------------
INPUT_FILE=$1
JOB_ID=$2
OUTPUT_ROOT="/storage/fdunphome/wangchunzheng/CWCoalProject/bin/condor/outputs/"
SAVEDIR="$OUTPUT_ROOT/job_${JOB_ID}"  # 独立工作目录

# ----------------------------
# 4. 目录准备
# ----------------------------
mkdir -p "$SAVEDIR"

# ----------------------------
# 5. 执行主程序
# ----------------------------
echo "=== JOB STARTED [$(date)] ==="
echo "Input file: $INPUT_FILE"
echo "Job ID    : $JOB_ID"
echo "Workdir   : $SAVEDIR"

/storage/fdunphome/wangchunzheng/CWCoalProject/bin/cwcoal \
  -i "$INPUT_FILE" \
  -a KDTreeGlobal \
  -s "$SAVEDIR" \
  -r 1.5

# ----------------------------
# 6. 文件转移和清理
# ----------------------------
if [ -f "$SAVEDIR/cve_KDTreeGlobal.root" ]; then
  mv "$SAVEDIR/cve_KDTreeGlobal.root" "$OUTPUT_ROOT/cve_KDTreeGlobal_${JOB_ID}.root"
fi
if [ -f "$SAVEDIR/qa_KDTreeGlobal.root" ]; then
  mv "$SAVEDIR/qa_KDTreeGlobal.root" "$OUTPUT_ROOT/qa_KDTreeGlobal_${JOB_ID}.root"
fi
rm -rf "$SAVEDIR"

# ----------------------------
# 7. 环境清理
# ----------------------------
conda deactivate
echo "=== JOB COMPLETED [$(date)] ==="
