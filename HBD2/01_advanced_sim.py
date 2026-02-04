# 01_advanced_sim.py
# Step 1: 复杂文库模拟 (Simulation)
# 核心改进：模拟 P5/P7 结构、不等长 Insert（随机断裂）、质量值衰减，以及 3 个极端的“高表达 PCR 家族”。
import random

def generate_complex_reads():
    # 模拟接头与标签
    P5, P7 = "AATGATACGGCGACCACCGAGATCT", "CAAGCAGAAGACGGCATACGAGAT"
    S_INDEX, ANCHOR = "ATGCATGC", "GGG"
    
    with open("raw_R1.fastq", "w") as f1, open("raw_R2.fastq", "w") as f2:
        for i in range(1000):
            # 1. 模拟随机剪切产生的插入片段 (150-400bp)
            insert_len = random.randint(150, 400)
            
            # 2. 模拟“超大家族”与随机分子
            if i < 300: # 模拟 3 个高表达分子，每个扩增了 100 次
                family_id = i // 100 
                umi = f"UMI{family_id:03}"
                full_insert = "G" * insert_len # 假设由于结构问题导致的低复杂度序列
            else:
                umi = "".join(random.choices("ATCG", k=6))
                full_insert = "".join(random.choices("ATCG", k=insert_len))
            
            # 3. 构造 Read1 (PE150): 包含 Index + UMI + 锚点
            r1_seq = (S_INDEX + umi + ANCHOR + full_insert)[:150]
            # 构造 Read2 (PE150): 从另一端（随机断裂处）反向读取
            r2_seq = full_insert[::-1][:150]
            
            # 4. 模拟质量值 (前端高、后端随长度增加而衰减)
            def get_qual(length):
                return "".join([chr(random.randint(60 if j > 100 else 70, 74)) for j in range(length)])

            f1.write(f"@READ_{i}\n{r1_seq}\n+\n{get_qual(len(r1_seq))}\n")
            f2.write(f"@READ_{i}\n{r2_seq}\n+\n{get_qual(len(r2_seq))}\n")

if __name__ == "__main__":
    generate_complex_reads()
    print("Step 1 完成：已生成含超大家族、不等长片段及真实质量值的 Fastq。")