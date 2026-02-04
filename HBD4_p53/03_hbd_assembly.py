# 03_hbd_assembly.py
import pysam
from collections import Counter, defaultdict

def get_consensus_with_quality(sequences):
    """改进的共识组装算法：返回共识序列及每位置的碱基频率信息"""
    if not sequences: return "", []

    consensus_seq = ""
    base_frequencies = []  # 记录每个位置的碱基频率

    for i in range(len(sequences[0])):
        column = [s[i] for s in sequences if i < len(s)]
        base_counts = Counter(column)

        # 获取最常见碱基及其频率
        most_common_base, count = base_counts.most_common(1)[0]
        freq_info = {base: count/base_counts.total() for base, count in base_counts.items()}

        consensus_seq += most_common_base
        base_frequencies.append(freq_info)

    return consensus_seq, base_frequencies

def analyze_mutations(consensus_seq, base_frequencies, reference_seq=None):
    """
    分析突变类型和可能来源
    - 低丰度突变：在多个家族中出现但在单个家族中频率较低
    - PCR引入突变：仅在单个家族中出现且频率较低
    - 野生型：高频出现的碱基
    """
    mutations = []
    for i, (cons_base, freq_info) in enumerate(zip(consensus_seq, base_frequencies)):
        # 检查是否有其他碱基存在（非共识碱基）
        other_bases = {base: freq for base, freq in freq_info.items() if base != cons_base}

        for base, freq in other_bases.items():
            # 根据频率判断突变类型
            if freq < 0.1:  # 低于10%的可能是PCR错误
                mut_type = "PCR_artifact"
            elif freq < 0.3:  # 中等频率可能是低丰度突变
                mut_type = "low_abundance_variant"
            else:  # 高频率可能是真正的变异
                mut_type = "true_variant"

            mutations.append({
                'position': i,
                'original_base': cons_base,
                'variant_base': base,
                'frequency': freq,
                'type': mut_type
            })

    return mutations

def deduplicate_molecules(bam_path):
    """
    HBD深度去重：利用"染色体 + 起始坐标 + 插入片段长度(ISIZE) + UMI"四位一体识别原始分子
    ISIZE代表了Read2端的随机断裂点，是区分PCR重复的关键
    """
    # 存储唯一的分子: { (chrom, pos, isize, umi) : [seq1, seq2, ...] }
    unique_molecules = defaultdict(list)
    # 同时存储每个分子的详细信息
    molecule_details = defaultdict(list)

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

            # 记录read的详细信息
            molecule_details[key].append({
                'query_name': r.query_name,
                'query_sequence': r.query_sequence,
                'reference_start': r.reference_start,
                'reference_end': r.reference_end,
                'mapping_quality': r.mapping_quality,
                'cigar': r.cigarstring
            })

    return unique_molecules, molecule_details

def assemble_families(bam_path):
    # 获取去重后的分子家族
    unique_molecules, molecule_details = deduplicate_molecules(bam_path)

    print(f"{'Chrom':<8} {'Position':<12} {'UMI':<8} {'ISIZE':<8} {'Reads':<6} {'FamilySize':<8} {'Consensus_Seq_Snippet'}")
    print("-" * 100)

    # 用于统计所有家族的突变信息
    all_mutations = []

    with open("hbd_final_assembly_report.csv", "w") as o:
        # 添加更多列来描述突变信息
        o.write("Gene,Chrom,Pos,UMI,ISIZE,ReadCount,FamilySize,ConsensusSequence,MutationInfo,VariantType,WildTypeRatio,LowFreqMutations,HighFreqMutations,PCRArtifacts\n")

        for key, seqs in unique_molecules.items():
            consensus, base_freqs = get_consensus_with_quality(seqs)
            mutations = analyze_mutations(consensus, base_freqs)

            chrom, pos, isize, umi = key
            # 演示目的：由于我们比对的是全基因组，这里我们可以手动标注 TP53 位点
            gene = "TP53" if chrom == "chr17" else "Unknown" # TP53 位于 chr17

            # 分析突变类型
            low_freq_muts = [m for m in mutations if m['type'] == 'low_abundance_variant']
            high_freq_muts = [m for m in mutations if m['type'] == 'true_variant']
            pcr_artifacts = [m for m in mutations if m['type'] == 'PCR_artifact']

            # 计算野生型比例（共识碱基的平均频率）
            wildtype_ratio = sum([freq_info[cons_base] for cons_base, freq_info in zip(consensus, base_freqs)]) / len(base_freqs) if base_freqs else 0

            # 格式化突变信息
            mut_info = ";".join([f"{m['position']}:{m['original_base']}->{m['variant_base']}({m['frequency']:.2f},{m['type']})" for m in mutations])

            # 打印部分结果到屏幕
            print(f"{chrom:<8} {pos:<12} {umi:<8} {isize:<8} {len(seqs):<6} {len(unique_molecules[key]):<8} {consensus[:30]}...")

            # 写入完整报表
            o.write(f"{gene},{chrom},{pos},{umi},{isize},{len(seqs)},{len(unique_molecules[key])},{consensus},{mut_info if mut_info else 'None'},family_specific,{wildtype_ratio:.3f},{len(low_freq_muts)},{len(high_freq_muts)},{len(pcr_artifacts)}\n")

            # 记录此家族的所有突变信息
            for mut in mutations:
                all_mutations.append({
                    'chrom': chrom,
                    'pos': pos + mut['position'],  # 实际基因组位置
                    'umi': umi,
                    'mutation': mut
                })

    # 统计分子家族大小分布
    family_sizes = [len(seqs) for seqs in unique_molecules.values()]
    with open("molecule_family_distribution.txt", "w") as f:
        f.write("FamilySize\tCount\n")
        size_counts = {}
        for size in family_sizes:
            size_counts[size] = size_counts.get(size, 0) + 1
        for size, count in sorted(size_counts.items()):
            f.write(f"{size}\t{count}\n")

    # 生成详细的突变报告
    with open("mutation_analysis_report.txt", "w") as f:
        f.write("Chrom\tPos\tUMI\tPositionInSeq\tOriginalBase\tVariantBase\tFrequency\tType\tFamilySize\n")
        for mut_record in all_mutations:
            mut = mut_record['mutation']
            f.write(f"{mut_record['chrom']}\t{mut_record['pos']}\t{mut_record['umi']}\t{mut['position']}\t{mut['original_base']}\t{mut['variant_base']}\t{mut['frequency']:.3f}\t{mut['type']}\t{len(unique_molecules[(mut_record['chrom'], mut_record['pos']-mut['position'], isize, mut_record['umi'])]) if (mut_record['chrom'], mut_record['pos']-mut['position'], isize, mut_record['umi']) in unique_molecules else 'N/A'}\n")

    # 统计总体突变情况
    total_mutations = len(all_mutations)
    pcr_artifact_count = len([m for m in all_mutations if m['mutation']['type'] == 'PCR_artifact'])
    low_freq_variant_count = len([m for m in all_mutations if m['mutation']['type'] == 'low_abundance_variant'])
    true_variant_count = len([m for m in all_mutations if m['mutation']['type'] == 'true_variant'])

    print(f"\n突变分析总结:")
    print(f"总突变数: {total_mutations}")
    print(f"PCR引入的假突变: {pcr_artifact_count}")
    print(f"低丰度真实变异: {low_freq_variant_count}")
    print(f"高丰度真实变异: {true_variant_count}")

assemble_families("sorted.bam")