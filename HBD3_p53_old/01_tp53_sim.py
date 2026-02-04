# 01_advanced_sim.py
import random

def generate_complex_reads():
    P5, P7 = "AATGATACGGCGACCACCGAGATCT", "CAAGCAGAAGACGGCATACGAGAT"
    S_INDEX, ANCHOR = "ATGCATGC", "GGG"
    
    with open("raw_R1.fastq", "w") as f1, open("raw_R2.fastq", "w") as f2:
        for i in range(1000):
            # 1. 模拟片段长度不等 (150bp - 400bp)
            insert_len = random.randint(150, 400)
            
            # 2. 模拟“超大家族” (PCR 冗余)
            if i < 300: # 前 300 条属于 3 个超高表达基因
                family_id = i // 100 
                umi = f"UMI{family_id:03}"
                full_insert = "G" * insert_len # 假设高表达基因
            else:
                umi = "".join(random.choices("ATCG", k=6))
                full_insert = "".join(random.choices("ATCG", k=insert_len))
            
            # 3. 构造 Read1 (PE150): P5 + Index + UMI + Anchor + Insert...
            # 注意：测序只读 150bp
            r1_seq = (S_INDEX + umi + ANCHOR + full_insert)[:150]
            # 构造 Read2 (PE150): 从另一端反向互补读
            r2_seq = full_insert[::-1][:150]
            
            # 4. 模拟真实质量值 (不是全是 I，而是前端好后端差)
            def get_qual(length):
                return "".join([chr(random.randint(60, 74)) for _ in range(length)])

            f1.write(f"@READ_{i}\n{r1_seq}\n+\n{get_qual(len(r1_seq))}\n")
            f2.write(f"@READ_{i}\n{r2_seq}\n+\n{get_qual(len(r2_seq))}\n")

if __name__ == "__main__":
    generate_complex_reads()
    print("Step 1 完成：已生成含超大家族与真实质量值的模拟数据。")