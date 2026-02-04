#!/bin/bash
# 02_align.sh
REF="/home/gao/references/hg38/hg38.fa"

# 提取索引、UMI并清洗 - 通过锚点GGG动态定位UMI
python3 -c "
import re

def extract_indices_and_umi(seq):
    # HBD文库结构: primary_index(8nt) + sample_index(8nt) + UMI(6nt) + anchor(GGG) + insert
    # 验证序列长度至少为8+8+6+3=25
    if len(seq) < 25:
        return None, None, None, seq
    
    # 提取8bp一级i5 index
    primary_index = seq[:8]
    # 提取8bp二级样本index  
    sample_index = seq[8:16]
    # 提取6bp UMI
    umi = seq[16:22]
    
    # 验证GGG锚点位置
    anchor_pos = seq.find('GGG', 22)  # GGG应该在UMI之后
    if anchor_pos == 22:  # GGG紧接在UMI后
        cleaned_seq = seq[anchor_pos+3:]  # +3跳过GGG
    else:
        # 如果GGG不在预期位置，尝试查找最近的GGG
        anchor_pos = seq.find('GGG', 16)
        if anchor_pos != -1 and anchor_pos >= 16:
            umi = seq[anchor_pos-6:anchor_pos]  # 重新提取UMI
            cleaned_seq = seq[anchor_pos+3:]
        else:
            # 备用方案：使用默认位置
            cleaned_seq = seq[22:]  # 跳过索引+UMI部分
    
    return primary_index, sample_index, umi, cleaned_seq

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

            # 提取索引、UMI并清洗序列
            primary_idx, sample_idx, umi, cleaned_seq = extract_indices_and_umi(seq)
            
            if umi is None:
                # 如果无法提取UMI，跳过此read
                continue
                
            # 新的ID格式包含所有关键信息
            new_id = f'{id}_{primary_idx}_{sample_idx}_{umi}'

            # 写入清洗后的R1
            if len(cleaned_seq) > 0:
                o1.write(f'{new_id}\n{cleaned_seq}\n+\n{qual[len(seq)-len(cleaned_seq):]}\n')
                # R2保持不变，但ID加上完整信息
                o2.write(f'{id2}_{primary_idx}_{sample_idx}_{umi}\n{seq2}\n+\n{qual2}\n')
"

# BWA 比对 (2线程)
bwa mem -t 2 $REF clean_R1.fq clean_R2.fq | samtools view -Sb - | samtools sort -@ 2 -o sorted.bam -
samtools index sorted.bam