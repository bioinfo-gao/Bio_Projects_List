#!/bin/bash

# 1. 定义环境路径 (你的 R44 绝对路径)
R44_PATH="/home/gao/micromamba/envs/R44"

echo "========================================"
echo "   HBD High-Fidelity Analysis Pipeline"
echo "========================================"

# 2. 绕过损坏的 conda base，直接使用路径激活
source /home/gao/anaconda3/etc/profile.d/conda.sh
conda activate $R44_PATH

if [ $? -eq 0 ]; then
    echo "[SUCCESS] R44 Environment Activated."
else
    echo "[ERROR] Failed to activate R44. Please check the path."
    exit 1
fi

# 3. 运行 Python 核心分析脚本
# 假设你的主程序叫 main_hbd.py
echo "[RUNNING] Starting HBD Deduplication Analysis..."
python src/main_hbd.py

echo "[DONE] Analysis complete. Results are in results/ folder."
