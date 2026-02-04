# 04_hbd_dedup.py
import pysam
from collections import Counter

def hbd_analysis():
    inventory = Counter()
    with pysam.AlignmentFile("sorted.bam", "rb") as bam:
        for r in bam.fetch():
            if r.is_unmapped or not r.is_paired: continue
            
            # 从 ID 中取回 UMI
            umi = r.query_name.split("_")[-1]
            # 物理指纹：染色体 + 起始点 + 插入片段长度 + UMI
            fingerprint = (r.reference_name, r.reference_start, abs(r.template_length), umi)
            inventory[fingerprint] += 1
            
    with open("hbd_final_results.csv", "w") as f:
        f.write("Chrom,Pos,FragSize,UMI,ReadCount\n")
        for k, v in inventory.items():
            f.write(f"{k[0]},{k[1]},{k[2]},{k[3]},{v}\n")

if __name__ == "__main__":
    hbd_analysis()
    print("Step 4 完成：HBD 去重定量表已生成。")