# 01_tp53_sim.py
import random
from collections import defaultdict

def introduce_sequencing_errors(sequence, error_rate=0.01):
    """
    在序列中引入测序错误，包括替换、插入和删除
    """
    result = []
    i = 0
    while i < len(sequence):
        r = random.random()
        if r < error_rate:  # 替换错误
            # 用其他三个碱基之一替换当前碱基
            bases = ['A', 'T', 'C', 'G']
            bases.remove(sequence[i])
            result.append(random.choice(bases))
        elif r < error_rate * 1.2:  # 插入错误
            result.append(random.choice(['A', 'T', 'C', 'G']))
            # 保持i不变，重新处理当前碱基
            continue
        elif r < error_rate * 1.4:  # 删除错误
            # 不添加当前碱基，跳过它
            pass
        else:
            result.append(sequence[i])
        i += 1

    return ''.join(result)

def generate_random_umi(length=6):
    """生成随机UMI序列"""
    return ''.join(random.choices(['A', 'T', 'C', 'G'], k=length))

def generate_quality_scores(length):
    """生成Phred质量分数字符串"""
    # 模拟质量分数，范围通常在33-73对应ASCII字符'!'-'K'
    return ''.join([chr(random.randint(33, 73)) for _ in range(length)])

def generate_sample_indexes(count=8):
    """生成二级样本index，模拟不同的样本"""
    sample_indexes = []
    for i in range(count):
        # 生成8碱基的二级样本index
        sample_index = ''.join(random.choices(['A', 'T', 'C', 'G'], k=8))
        sample_indexes.append(sample_index)
    return sample_indexes

def generate_tp53_data():
    # TP53 某外显子真实序列片段
    tp53_seq = "TGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATG"

    # 根据文档，文库结构包含：
    # 8碱基的一级5端index (i5)
    # 3碱基的一级7端index (i7)
    # 8碱基的UMI锚定序列 (N8)
    # 6碱基的UMI
    # 特异性锚定序列
    # 插入序列 (目标序列)
    i5_index = "ATGCATGC"  # 8碱基的一级5端index
    i7_index = "TAC"       # 3碱基的一级7端index

    # 生成多个二级样本index，模拟不同样本的数据
    sample_indexes = generate_sample_indexes(5)  # 假设有5个不同样本

    # 创建不同大小的分子家族，模拟真实的分子分布
    families = []
    family_count = 50  # 总共创建50个不同的分子家族
    total_reads = 1000

    # 生成家族大小分布（模拟真实的分子丰度分布）
    family_sizes = []
    remaining_reads = total_reads
    for i in range(family_count - 1):
        # 随机分配read数量给每个家族，遵循长尾分布
        size = max(1, int(random.expovariate(0.5)))  # 指数分布，偏向小家族
        if size > remaining_reads:
            size = remaining_reads
        family_sizes.append(size)
        remaining_reads -= size

    family_sizes.append(remaining_reads)  # 把剩余的reads分配给最后一个家族

    # 生成每个家族的UMI和序列
    for fam_idx, fam_size in enumerate(family_sizes):
        # 从预定义的样本index中随机选择一个
        sample_index = random.choice(sample_indexes)

        # 生成6碱基UMI
        umi = generate_random_umi(6)

        # 为每个家族生成一个独特的原始序列（基于TP53序列的一个片段）
        start_pos = random.randint(0, len(tp53_seq) - 100)
        original_seq = tp53_seq[start_pos:start_pos + 100]

        # 对于家族中的每个read，模拟PCR扩增和测序错误
        for read_idx in range(fam_size):
            # 从原始序列开始，引入PCR扩增错误
            amplified_seq = introduce_sequencing_errors(original_seq, error_rate=0.005)

            # 最后引入测序错误
            final_seq = introduce_sequencing_errors(amplified_seq, error_rate=0.01)

            # 模拟随机断裂，使每个分子的read2起始位置略有不同
            # 这是关键：read2的起始位置代表随机断裂点，这是UMI家族组装的基础
            random_fragment_start = random.randint(0, len(final_seq) - 50)  # 至少50bp长度
            insert = final_seq[random_fragment_start:]

            # 构建R1序列 (包含所有文库组件)
            # 文档中的文库结构：i5_index + i7_index + sample_index + N8 + UMI + anchor + insert
            r1_seq = f"{i5_index}{i7_index}{sample_index}{'N' * 8}{umi}GGG{insert}"

            # R2是插入序列的反向互补（代表read2，用于确定随机断裂点）
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            r2_seq = ''.join([complement[base] for base in insert])[::-1]

            yield r1_seq, r2_seq

def main():
    with open("raw_R1.fastq", "w") as f1, open("raw_R2.fastq", "w") as f2:
        read_counter = 0
        for r1_seq, r2_seq in generate_tp53_data():
            # 生成质量分数
            r1_qual = generate_quality_scores(len(r1_seq))
            r2_qual = generate_quality_scores(len(r2_seq))

            # 写入FASTQ格式
            f1.write(f"@READ_{read_counter}\n{r1_seq}\n+\n{r1_qual}\n")
            f2.write(f"@READ_{read_counter}\n{r2_seq}\n+\n{r2_qual}\n")

            read_counter += 1

    print(f"Step 1: TP53 复杂数据模拟完成，共生成 {read_counter} 对 reads。")
    print("文库结构: i5_index(8nt) + i7_index(3nt) + sample_index(8nt) + UMI_anchor(N8) + UMI(6nt) + anchor(GGG) + insert_seq")
    print("分析流程: 数据拆分 -> 样本拆分 -> UMI家族统计 -> 分子组装")

if __name__ == "__main__":
    main()