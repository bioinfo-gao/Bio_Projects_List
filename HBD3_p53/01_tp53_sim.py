# 01_tp53_sim.py
import random

def generate_tp53_data():
    # TP53 某外显子真实序列片段
    tp53_seq = "TGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATG"
    S_INDEX, ANCHOR = "ATGCATGC", "GGG"
    
    with open("raw_R1.fastq", "w") as f1, open("raw_R2.fastq", "w") as f2:
        for i in range(1000):
            # 模拟 3 个超大家族，每个家族 200 条 Reads，用来展示组装效果
            if i < 600:
                family_id = i // 200
                umi = f"UMI00{family_id}"
                # 模拟随机断裂：同一家族起始位点相同，但长度略有不同（由于剪切）
                start_offset = 0 
                insert = tp53_seq[start_offset : start_offset + 100]
                # 人为在第 20 位引入一个随机测序错误 (10% 概率)
                if random.random() < 0.1:
                    list_ins = list(insert)
                    list_ins[20] = "N" 
                    insert = "".join(list_ins)
            else:
                umi = "".join(random.choices("ATCG", k=6))
                insert = "".join(random.choices("ATCG", k=100))
            
            f1.write(f"@READ_{i}\n{S_INDEX}{umi}{ANCHOR}{insert}\n+\n{'I'*len(S_INDEX+umi+ANCHOR+insert)}\n")
            f2.write(f"@READ_{i}\n{insert[::-1]}\n+\n{'I'*len(insert)}\n")

generate_tp53_data()
print("Step 1: TP53 复杂数据模拟完成。")