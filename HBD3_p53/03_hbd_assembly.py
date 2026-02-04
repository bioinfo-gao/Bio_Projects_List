# 03_hbd_assembly.py
import pysam
from collections import Counter, defaultdict

def get_consensus(sequences):
    """简单共识组装算法：每一位取出现频率最高的碱基"""
    if not sequences: return ""
    res = ""
    for i in range(len(sequences[0])):
        column = [s[i] for s in sequences if i < len(s)]
        res += Counter(column).most_common(1)[0][0]
    return res

def deduplicate_molecules(bam_path):
    """
    HBD深度去重：利用"染色体 + 起始坐标 + 插入片段长度(ISIZE) + UMI"四位一体识别原始分子
    ISIZE代表了Read2端的随机断裂点，是区分PCR重复的关键
    """
    # 存储唯一的分子: { (chrom, pos, isize, umi) : [seq1, seq2, ...] }
    unique_molecules = defaultdict(list)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            if r.is_unmapped: continue

            # 从查询名称中提取UMI
            umi = r.query_name.split("_")[-1]

            # 获取插入片段长度(ISIZE) - 这是区分PCR重复的关键
            isize = r.template_length  # TLEN字段，即插入片段长度

            # 四位一体键：染色体 + 起始坐标 + ISIZE + UMI
            key = (r.reference_name, r.reference_start, isize, umi)

            unique_molecules[key].append(r.query_sequence)

    return unique_molecules

def assemble_families(bam_path):
    # 获取去重后的分子家族
    unique_molecules = deduplicate_molecules(bam_path)

    print(f"{'Chrom':<8} {'Position':<12} {'UMI':<8} {'ISIZE':<8} {'Reads':<6} {'Consensus_Seq_Snippet'}")
    print("-" * 80)

    with open("hbd_final_assembly_report.csv", "w") as o:
        o.write("Gene,Chrom,Pos,UMI,ISIZE,ReadCount,ConsensusSequence\n")
        for key, seqs in unique_molecules.items():
            consensus = get_consensus(seqs)
            chrom, pos, isize, umi = key
            # 演示目的：由于我们比对的是全基因组，这里我们可以手动标注 TP53 位点
            gene = "TP53" if chrom == "chr17" else "Unknown" # TP53 位于 chr17

            # 打印部分结果到屏幕
            print(f"{chrom:<8} {pos:<12} {umi:<8} {isize:<8} {len(seqs):<6} {consensus[:30]}...")

            # 写入完整报表
            o.write(f"{gene},{chrom},{pos},{umi},{isize},{len(seqs)},{consensus}\n")

    # 统计分子家族大小分布
    family_sizes = [len(seqs) for seqs in unique_molecules.values()]
    with open("molecule_family_distribution.txt", "w") as f:
        f.write("FamilySize\tCount\n")
        size_counts = {}
        for size in family_sizes:
            size_counts[size] = size_counts.get(size, 0) + 1
        for size, count in sorted(size_counts.items()):
            f.write(f"{size}\t{count}\n")

assemble_families("sorted.bam")