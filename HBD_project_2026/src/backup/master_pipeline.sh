#!/bin/bash

# ============================================================
# HBD High-Fidelity Analysis Master Pipeline
# Created for: 2026-02-04 Project Test
# ============================================================

# 设置错误即停止
set -e

# 1. 环境激活 (针对你的 micromamba R44 路径)
R44_ENV="/home/gao/micromamba/envs/R44"
echo ">>> [STEP 0] 激活 R44 环境..."
source /home/gao/anaconda3/etc/profile.d/conda.sh
conda activate $R44_ENV

# 定义参考基因组路径 (请确保该路径下已完成 bwa index)
REF_GENOME="/home/gao/references/hg38/hg38.fa"

# 2. 数据模拟 (Python)
echo ">>> [STEP 1] 正在模拟 1000 条 HBD Fastq 数据..."
python3 -c '
import random
def generate():
    with open("sim_R1.fastq", "w") as f1, open("sim_R2.fastq", "w") as f2:
        for i in range(1000):
            # 模拟 15% 的 PCR 重复 (固定 UMI 和序列)
            if i % 100 < 15:
                umi, seq = "AAAAAA", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
            else:
                umi = "".join(random.choices("ATCG", k=6))
                seq = "".join(random.choices("ATCG", k=32))
            f1.write(f"@READ_{i}\nATGCATGC{umi}GGG{seq}\n+\n{"I"*49}\n")
            f2.write(f"@READ_{i}\n{seq[::-1]}\n+\n{"I"*32}\n")
generate()
'

# 3. 预处理与 BWA 比对 (Bash/BWA/Samtools)
echo ">>> [STEP 2] 正在进行 UMI 提取与全基因组比对..."
# 提取 UMI 并切除人工序列
python3 -c "
with open('sim_R1.fastq') as f1, open('sim_R2.fastq') as f2, open('clean_R1.fq', 'w') as o1, open('clean_R2.fq', 'w') as o2:
    for line in f1:
        id = line.strip(); seq = next(f1).strip(); next(f1); qual = next(f1).strip()
        id2 = next(f2).strip(); seq2 = next(f2).strip(); next(f2); qual2 = next(f2).strip()
        umi = seq[8:14]
        new_id = f'{id}_{umi}'
        o1.write(f'{new_id}\n{seq[17:]}\n+\n{qual[17:]}\n')
        o2.write(f'{new_id}\n{seq2}\n+\n{qual2}\n')
"

# 执行比对
bwa mem -t 8 $REF_GENOME clean_R1.fq clean_R2.fq | \
samtools view -Sb - | \
samtools sort -o sorted.bam -
samtools index sorted.bam

# 4. HBD 核心去重算法 (Python/Pysam)
echo ">>> [STEP 3] 正在执行 HBD 物理坐标 + UMI 联合去重..."
python3 -c '
import pysam
from collections import Counter
def dedup():
    inventory = Counter()
    with pysam.AlignmentFile("sorted.bam", "rb") as sam:
        for r in sam.fetch():
            if r.is_unmapped or not r.is_paired: continue
            umi = r.query_name.split("_")[-1]
            fingerprint = (r.reference_name, r.reference_start, abs(r.template_length), umi)
            inventory[fingerprint] += 1
    with open("hbd_counts.csv", "w") as f:
        f.write("chrom,pos,isize,umi,count\n")
        for k, v in inventory.items():
            f.write(f"{k[0]},{k[1]},{k[2]},{k[3]},{v}\n")
dedup()
'

# 5. 可视化结果 (R)
echo ">>> [STEP 4] 正在生成分析报告 (PDF)..."
Rscript -e '
library(ggplot2)
df <- read.csv("hbd_counts.csv")
p <- ggplot(df, aes(x=count)) + 
     geom_histogram(binwidth=1, fill="firebrick", color="white") + 
     scale_y_log10() + 
     labs(title="HBD Quantitative Report", x="Reads per Molecule (Duplicates)", y="Count") +
     theme_minimal()
ggsave("HBD_Fidelity_Analysis.pdf", p)
'

echo "============================================================"
echo "流水线运行成功！"
echo "最终定量结果: hbd_counts.csv"
echo "可视化报告: HBD_Fidelity_Analysis.pdf"
echo "============================================================"
