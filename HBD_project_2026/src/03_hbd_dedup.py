# 文件名: 03_hbd_dedup.py
import pysam
from collections import Counter

def hbd_dedup(bam_file):
    inventory = Counter()
    
    with pysam.AlignmentFile(bam_file, "rb") as sam:
        for read in sam.fetch():
            if read.is_unmapped or not read.is_paired: continue
            
            # 提取 UMI
            umi = read.query_name.split('_')[-1]
            # 获取物理位置指纹 (染色体, 起始, 片段跨度, UMI)
            fingerprint = (read.reference_name, read.reference_start, abs(read.template_length), umi)
            inventory[fingerprint] += 1
            
    # 输出去重后的矩阵
    with open("counts.csv", "w") as f:
        f.write("chrom,pos,isize,umi,count\n")
        for k, v in inventory.items():
            f.write(f"{k[0]},{k[1]},{k[2]},{k[3]},{v}\n")

if __name__ == "__main__":
    hbd_dedup("sorted.bam")
    print("Step 3: HBD 去重完成，结果存入 counts.csv。")