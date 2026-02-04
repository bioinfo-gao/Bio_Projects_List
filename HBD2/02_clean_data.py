# 02_clean_data.py
def clean_hbd():
    with open("raw_R1.fastq") as f1, open("raw_R2.fastq") as f2, \
         open("clean_R1.fq", "w") as o1, open("clean_R2.fq", "w") as o2:
        for line in f1:
            header1 = line.strip()
            seq1 = next(f1).strip()
            plus1 = next(f1).strip()
            qual1 = next(f1).strip()
            
            header2 = next(f2).strip()
            seq2 = next(f2).strip()
            plus2 = next(f2).strip()
            qual2 = next(f2).strip()
            
            # 动态搜索 GGG 锚点
            anchor_pos = seq1.find("GGG")
            if anchor_pos >= 6:
                umi = seq1[anchor_pos-6 : anchor_pos]
                new_id = f"{header1.split()[0]}_{umi}"
                
                # 切除锚点及其之前的所有序列 (Index+UMI+GGG)
                clean_seq1 = seq1[anchor_pos+3:]
                clean_qual1 = qual1[anchor_pos+3:]
                
                o1.write(f"{new_id}\n{clean_seq1}\n+\n{clean_qual1}\n")
                o2.write(f"{new_id}\n{seq2}\n+\n{qual2}\n")

if __name__ == "__main__":
    clean_hbd()
    print("Step 2 完成：UMI 已提取至 ID，接头已切除。")