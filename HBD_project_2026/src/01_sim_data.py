# 文件名: 01_sim_data.py
import random
import gzip

def generate_hbd_reads(num_reads=1000):
    sample_index = "ATGCATGC" # 8bp 二级 Index
    anchor = "GGG"             # 3bp 锚点
    
    with open("sim_R1.fastq", "w") as f1, open("sim_R2.fastq", "w") as f2:
        for i in range(num_reads):
            # 模拟 PCR 重复：每 100 条里有 15 条是完全一样的
            if i % 100 < 15:
                umi = "AAAAAA"
                insert_seq = "TGCGATACGACTAGCTAGCTAGCTAGCTAGCT"
            else:
                umi = "".join(random.choices("ATCG", k=6))
                insert_seq = "".join(random.choices("ATCG", k=32))
            
            # Read1: Index + UMI + Anchor + BioSeq
            f1.write(f"@READ_{i}\n{sample_index}{umi}{anchor}{insert_seq}\n+\n{'I'*49}\n")
            # Read2: BioSeq (模拟从另一端随机断裂处读取)
            f2.write(f"@READ_{i}\n{insert_seq[::-1]}\n+\n{'I'*32}\n")

if __name__ == "__main__":
    generate_hbd_reads()
    print("Step 1: 模拟 Fastq 数据生成完成。")