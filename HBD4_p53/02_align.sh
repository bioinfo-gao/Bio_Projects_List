#!/bin/bash
# 02_align.sh
REF="/home/gao/references/hg38/hg38.fa"

# 提取 UMI 并清洗 - 通过锚点GGG动态定位UMI
python3 -c "
import re

def extract_umi_and_clean_seq(seq):
    # 通过锚点GGG动态定位UMI
    # 文库结构: i5_index(8nt) + i7_index(3nt) + sample_index(8nt) + N8 + UMI(6nt) + anchor(GGG) + insert
    # 查找GGG锚点位置
    anchor_pos = seq.find('GGG')
    if anchor_pos != -1 and anchor_pos >= 14:  # 确保有足够的空间容纳前面的序列
        # UMI在GGG前6个碱基
        umi_start = anchor_pos - 6
        umi = seq[umi_start:umi_start+6]
        # 清洗后的序列从GGG之后开始
        cleaned_seq = seq[anchor_pos+3:]  # +3跳过GGG
    else:
        # 如果没找到GGG锚点，使用默认位置提取UMI（备用方案）
        umi = seq[27:33]  # 8+3+8+N8=19, +6=25, +3=28, 所以前面应该是27:33
        cleaned_seq = seq[33:]  # 默认切割位置

    return umi, cleaned_seq

with open('raw_R1.fastq') as f, open('clean_R1.fq', 'w') as o1, open('clean_R2.fq', 'w') as o2:
    with open('raw_R2.fastq') as f2:
        for line in f:
            id = line.strip()
            seq = next(f).strip()
            next(f)
            qual = next(f).strip()

            # 读取R2
            id2 = next(f2).strip()
            seq2 = next(f2).strip()
            next(f2)
            qual2 = next(f2).strip()

            # 提取UMI并清洗序列
            umi, cleaned_seq = extract_umi_and_clean_seq(seq)
            new_id = f'{id}_{umi}'

            # 写入清洗后的R1
            o1.write(f'{new_id}\n{cleaned_seq}\n+\n{qual[len(seq)-len(cleaned_seq):]}\n')
            # R2保持不变，但ID加上UMI
            o2.write(f'{id2}_{umi}\n{seq2}\n+\n{qual2}\n')
"

# BWA 比对 (2线程)
bwa mem -t 2 $REF clean_R1.fq clean_R2.fq | samtools view -Sb - | samtools sort -@ 2 -o sorted.bam -
samtools index sorted.bam