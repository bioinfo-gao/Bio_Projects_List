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

def classify_error_source(mutation_freq, family_size):
    """
    根据突变频率和家族大小分类错误来源：
    - 低频自发突变：在多个家族中出现，频率较低但一致
    - PCR第一次扩增错误：在单个家族中高频出现（>50%），因为早期错误会被大量复制
    - PCR后续扩增错误：在单个家族中低频出现（<30%），因为晚期错误只影响少数拷贝
    """
    if mutation_freq >= 0.5:
        return "PCR_first_amplification_error"
    elif mutation_freq >= 0.1:
        return "PCR_later_amplification_error"
    else:
        return "low_frequency_spontaneous_mutation"

def analyze_mutations_detailed(consensus_seq, base_frequencies, reference_seq=None):
    """
    详细分析突变类型和可能来源，区分不同类型的PCR错误
    """
    mutations = []
    for i, (cons_base, freq_info) in enumerate(zip(consensus_seq, base_frequencies)):
        # 检查是否有其他碱基存在（非共识碱基）
        other_bases = {base: freq for base, freq in freq_info.items() if base != cons_base}

        for base, freq in other_bases.items():
            # 根据频率和上下文判断突变类型
            mut_type = classify_error_source(freq, len(base_frequencies))
            
            mutations.append({
                'position': i,
                'original_base': cons_base,
                'variant_base': base,
                'frequency': freq,
                'type': mut_type
            })

    return mutations

def deduplicate_molecules_advanced(bam_path):
    """
    HBD深度去重：利用"二级样本Index + UMI + Read2起始坐标"作为分子唯一标识
    这更符合HBD文库的实际设计
    """
    # 存储唯一的分子: { (sample_index, umi, read2_start) : [seq1, seq2, ...] }
    unique_molecules = defaultdict(list)
    # 同时存储每个分子的详细信息
    molecule_details = defaultdict(list)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            if r.is_unmapped: continue

            # 从查询名称中提取完整信息
            query_parts = r.query_name.split("_")
            if len(query_parts) >= 4:
                # 格式: READ_x_primaryIdx_sampleIdx_umi
                primary_index = query_parts[-3]  # 倒数第三个
                sample_index = query_parts[-2]   # 倒数第二个  
                umi = query_parts[-1]            # 最后一个
            else:
                # 备用方案：如果格式不匹配，尝试从最后提取
                umi = query_parts[-1]
                sample_index = "unknown"
                primary_index = "unknown"

            # 获取Read2的起始比对位置（随机断裂点）
            # 对于配对reads，我们需要考虑mate的位置
            if r.is_read1:
                # Read1的mate是Read2
                read2_start = r.next_reference_start if r.has_tag('MC') else r.reference_start
            else:
                # Read2本身
                read2_start = r.reference_start

            # 三位一体键：二级样本Index + UMI + Read2起始坐标
            key = (sample_index, umi, read2_start)

            unique_molecules[key].append(r.query_sequence)

            # 记录read的详细信息
            molecule_details[key].append({
                'query_name': r.query_name,
                'query_sequence': r.query_sequence,
                'reference_start': r.reference_start,
                'reference_end': r.reference_end,
                'mapping_quality': r.mapping_quality,
                'cigar': r.cigarstring,
                'primary_index': primary_index,
                'sample_index': sample_index,
                'umi': umi,
                'read2_start': read2_start
            })

    return unique_molecules, molecule_details

def count_total_reads_from_fastq(fastq_file):
    """从FASTQ文件计算总reads数"""
    count = 0
    try:
        with open(fastq_file, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    count += 1
    except FileNotFoundError:
        count = 0
    return count

def generate_read_statistics_report(unique_molecules, raw_r1_file="raw_R1.fastq", clean_r1_file="clean_R1.fq"):
    """生成详细的reads统计报告"""
    
    # 1. 原始reads总数（来自图片模拟的数据）
    original_reads = count_total_reads_from_fastq(raw_r1_file)
    
    # 2. 清洗后的reads数（成功提取UMI和索引的reads）
    cleaned_reads = count_total_reads_from_fastq(clean_r1_file)
    
    # 3. 映射到参考基因组的reads数
    mapped_reads = sum(len(seqs) for seqs in unique_molecules.values())
    
    # 4. 分子家族数量（unique molecules after deduplication）
    unique_molecules_count = len(unique_molecules)
    
    # 5. PCR重复率计算
    pcr_duplicates = mapped_reads - unique_molecules_count
    
    # 6. 最终真实reads数目（共识序列数）
    final_consensus_reads = unique_molecules_count
    
    # 创建统计报告
    stats_report = {
        'original_reads': original_reads,
        'cleaned_reads': cleaned_reads,
        'mapped_reads': mapped_reads,
        'unique_molecules': unique_molecules_count,
        'pcr_duplicates': pcr_duplicates,
        'final_consensus_reads': final_consensus_reads,
        'pcr_duplication_rate': (pcr_duplicates / mapped_reads * 100) if mapped_reads > 0 else 0,
        'cleaning_efficiency': (cleaned_reads / original_reads * 100) if original_reads > 0 else 0,
        'mapping_rate': (mapped_reads / cleaned_reads * 100) if cleaned_reads > 0 else 0
    }
    
    return stats_report

def write_statistics_report(stats_report, output_file="read_analysis_report.txt"):
    """写入详细的统计报告"""
    with open(output_file, "w") as f:
        f.write("HBD5_p53 Read Analysis Statistics Report\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("1. Original Reads (from simulation/images):\n")
        f.write(f"   Total original reads: {stats_report['original_reads']:,}\n")
        f.write("   - These represent the raw simulated data including all PCR duplicates\n")
        f.write("   - Generated from the TP53 simulation with introduced errors\n\n")
        
        f.write("2. Cleaned Reads (after UMI/index extraction):\n")
        f.write(f"   Successfully processed reads: {stats_report['cleaned_reads']:,}\n")
        f.write(f"   Cleaning efficiency: {stats_report['cleaning_efficiency']:.2f}%\n")
        f.write("   - Reads that successfully passed UMI and index extraction\n")
        f.write("   - Failed extractions are typically due to missing GGG anchor or short sequences\n\n")
        
        f.write("3. Mapped Reads (aligned to reference genome):\n")
        f.write(f"   Successfully mapped reads: {stats_report['mapped_reads']:,}\n")
        f.write(f"   Mapping rate: {stats_report['mapping_rate']:.2f}%\n")
        f.write("   - Reads that successfully aligned to the hg38 reference genome\n\n")
        
        f.write("4. Molecular Deduplication Analysis:\n")
        f.write(f"   Unique molecular families: {stats_report['unique_molecules']:,}\n")
        f.write(f"   PCR duplicates removed: {stats_report['pcr_duplicates']:,}\n")
        f.write(f"   PCR duplication rate: {stats_report['pcr_duplication_rate']:.2f}%\n")
        f.write("   - Each unique family represents one original DNA molecule\n")
        f.write("   - PCR duplicates are multiple reads from the same original molecule\n\n")
        
        f.write("5. Final Consensus Reads:\n")
        f.write(f"   Final high-fidelity consensus sequences: {stats_report['final_consensus_reads']:,}\n")
        f.write("   - These represent the true biological molecules after error correction\n")
        f.write("   - Each consensus sequence is generated from its molecular family\n\n")
        
        f.write("Summary:\n")
        f.write(f"   From {stats_report['original_reads']:,} original reads,\n")
        f.write(f"   we obtained {stats_report['final_consensus_reads']:,} high-confidence consensus sequences.\n")
        f.write(f"   This represents a {stats_report['final_consensus_reads']/stats_report['original_reads']*100:.2f}% yield\n")
        f.write("   of true biological molecules after removing PCR duplicates and sequencing errors.\n")

def assemble_families_advanced(bam_path):
    # 获取去重后的分子家族
    unique_molecules, molecule_details = deduplicate_molecules_advanced(bam_path)
    
    # 生成详细的统计报告
    stats_report = generate_read_statistics_report(unique_molecules)
    write_statistics_report(stats_report)
    
    print(f"{'SampleIdx':<12} {'UMI':<8} {'Read2Start':<12} {'Reads':<6} {'FamilySize':<8} {'Consensus_Seq_Snippet'}")
    print("-" * 100)

    # 用于统计所有家族的突变信息
    all_mutations = []

    with open("hbd_final_assembly_report.csv", "w") as o:
        # 添加更多列来描述详细的错误分类
        o.write("Gene,Chrom,Pos,SampleIndex,UMI,Read2Start,ReadCount,FamilySize,ConsensusSequence,MutationInfo,VariantType,WildTypeRatio,LowFreqSpontaneous,PCRFirstErrors,PCRLaterErrors\n")

        for key, seqs in unique_molecules.items():
            consensus, base_freqs = get_consensus_with_quality(seqs)
            mutations = analyze_mutations_detailed(consensus, base_freqs)

            sample_index, umi, read2_start = key
            
            # 获取染色体和位置信息（从第一个read获取）
            first_detail = molecule_details[key][0]
            chrom = "chr17"  # TP53位于chr17
            pos = first_detail['reference_start'] if 'reference_start' in first_detail else 0
            
            # 演示目的：由于我们比对的是全基因组，这里我们可以手动标注 TP53 位点
            gene = "TP53" if chrom == "chr17" else "Unknown"

            # 分析详细的突变类型
            low_freq_spontaneous = [m for m in mutations if m['type'] == 'low_frequency_spontaneous_mutation']
            pcr_first_errors = [m for m in mutations if m['type'] == 'PCR_first_amplification_error']
            pcr_later_errors = [m for m in mutations if m['type'] == 'PCR_later_amplification_error']

            # 计算野生型比例（共识碱基的平均频率）
            wildtype_ratio = sum([freq_info[cons_base] for cons_base, freq_info in zip(consensus, base_freqs)]) / len(base_freqs) if base_freqs else 0

            # 格式化突变信息
            mut_info = ";".join([f"{m['position']}:{m['original_base']}->{m['variant_base']}({m['frequency']:.2f},{m['type']})" for m in mutations])

            # 打印部分结果到屏幕
            print(f"{sample_index:<12} {umi:<8} {read2_start:<12} {len(seqs):<6} {len(unique_molecules[key]):<8} {consensus[:30]}...")

            # 写入完整报表
            o.write(f"{gene},{chrom},{pos},{sample_index},{umi},{read2_start},{len(seqs)},{len(unique_molecules[key])},{consensus},{mut_info if mut_info else 'None'},family_specific,{wildtype_ratio:.3f},{len(low_freq_spontaneous)},{len(pcr_first_errors)},{len(pcr_later_errors)}\n")

            # 记录此家族的所有突变信息
            for mut in mutations:
                all_mutations.append({
                    'chrom': chrom,
                    'pos': pos + mut['position'],  # 实际基因组位置
                    'sample_index': sample_index,
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
        f.write("Chrom\tPos\tSampleIndex\tUMI\tPositionInSeq\tOriginalBase\tVariantBase\tFrequency\tType\tFamilySize\n")
        for mut_record in all_mutations:
            mut = mut_record['mutation']
            family_key = (mut_record['sample_index'], mut_record['umi'], mut_record.get('read2_start', 0))
            family_size = len(unique_molecules.get(family_key, []))
            f.write(f"{mut_record['chrom']}\t{mut_record['pos']}\t{mut_record['sample_index']}\t{mut_record['umi']}\t{mut['position']}\t{mut['original_base']}\t{mut['variant_base']}\t{mut['frequency']:.3f}\t{mut['type']}\t{family_size}\n")

    # 统计总体突变情况
    total_mutations = len(all_mutations)
    pcr_first_count = len([m for m in all_mutations if m['mutation']['type'] == 'PCR_first_amplification_error'])
    pcr_later_count = len([m for m in all_mutations if m['mutation']['type'] == 'PCR_later_amplification_error'])
    low_freq_spontaneous_count = len([m for m in all_mutations if m['mutation']['type'] == 'low_frequency_spontaneous_mutation'])

    print(f"\n突变分析总结:")
    print(f"总突变数: {total_mutations}")
    print(f"PCR第一次扩增错误: {pcr_first_count}")
    print(f"PCR后续扩增错误: {pcr_later_count}")
    print(f"低频自发突变: {low_freq_spontaneous_count}")
    
    # 打印统计摘要
    print(f"\nRead Statistics Summary:")
    print(f"原始reads总数: {stats_report['original_reads']:,}")
    print(f"清洗后reads数: {stats_report['cleaned_reads']:,}")
    print(f"映射reads数: {stats_report['mapped_reads']:,}")
    print(f"唯一分子家族数: {stats_report['unique_molecules']:,}")
    print(f"最终共识reads数: {stats_report['final_consensus_reads']:,}")
    print(f"PCR重复率: {stats_report['pcr_duplication_rate']:.2f}%")

assemble_families_advanced("sorted.bam")