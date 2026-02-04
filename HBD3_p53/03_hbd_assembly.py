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

def assemble_families(bam_path):
    # 存储家族信息: { (pos, umi) : [seq1, seq2, ...] }
    families = defaultdict(list)
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            if r.is_unmapped: continue
            umi = r.query_name.split("_")[-1]
            # 指纹：染色体, 起始位点, UMI
            key = (r.reference_name, r.reference_start, umi)
            families[key].append(r.query_sequence)
            
    print(f"{'Chrom':<8} {'Position':<12} {'UMI':<8} {'Reads':<6} {'Consensus_Seq_Snippet'}")
    print("-" * 70)
    
    with open("hbd_final_assembly_report.csv", "w") as o:
        o.write("Gene,Chrom,Pos,UMI,ReadCount,ConsensusSequence\n")
        for key, seqs in families.items():
            consensus = get_consensus(seqs)
            # 演示目的：由于我们比对的是全基因组，这里我们可以手动标注 TP53 位点
            gene = "TP53" if key[0] == "chr17" else "Unknown" # TP53 位于 chr17
            
            # 打印部分结果到屏幕
            print(f"{key[0]:<8} {key[1]:<12} {key[2]:<8} {len(seqs):<6} {consensus[:30]}...")
            
            # 写入完整报表
            o.write(f"{gene},{key[0]},{key[1]},{key[2]},{len(seqs)},{consensus}\n")

assemble_families("sorted.bam")