import gzip

def preprocess_hbd_fastq(r1_in, r2_in, r1_out, r2_out):
    """
    功能：提取 UMI 到 Read ID 中，并切除接头序列
    """
    with gzip.open(r1_in, 'rt') as f1, gzip.open(r2_in, 'rt') as f2, \
         open(r1_out, 'w') as o1, open(r2_out, 'w') as o2:
        
        while True:
            header1 = f1.readline().strip()
            if not header1: break
            seq1 = f1.readline().strip()
            plus1 = f1.readline().strip()
            qual1 = f1.readline().strip()
            
            header2 = f2.readline().strip()
            seq2 = f2.readline().strip()
            plus2 = f2.readline().strip()
            qual2 = f2.readline().strip()
            
            # 1. 提取 UMI (假设：8bp Index + 6bp UMI + 3bp GGG)
            # 你需要根据实际情况微调这些切片位置
            sample_index = seq1[0:8]
            umi = seq1[8:14]
            
            # 2. 修改 Read ID，将 UMI 存入其中 (BWA 会保留 ID)
            new_id = f"{header1.split()[0]}_{umi}"
            
            # 3. 切除人工标签后的序列 (从第 17 位开始是真实 RNA)
            clean_seq1 = seq1[17:]
            clean_qual1 = qual1[17:]
            
            o1.write(f"{new_id}\n{clean_seq1}\n+\n{clean_qual1}\n")
            o2.write(f"{header2.split()[0]}_{umi}\n{seq2}\n+\n{qual2}\n")

# 使用示例
# preprocess_hbd_fastq("raw_R1.fastq.gz", "raw_R2.fastq.gz", "clean_R1.fastq", "clean_R2.fastq")