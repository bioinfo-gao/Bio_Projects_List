# 文件名: 02_align.sh
# 假设你已经有了 hg38.fa 的 BWA 索引

# 1. 预处理：将 R1 中的 UMI 提取到 Read 名称中，方便后续读取
# 我们用简单脚本模拟这个过程（实际可用 sed 或 python）
python3 -c "
with open('sim_R1.fastq') as f1, open('sim_R2.fastq') as f2, open('clean_R1.fq', 'w') as o1, open('clean_R2.fq', 'w') as o2:
    for line in f1:
        id = line.strip(); seq = next(f1).strip(); next(f1); qual = next(f1).strip()
        id2 = next(f2).strip(); seq2 = next(f2).strip(); next(f2); qual2 = next(f2).strip()
        umi = seq[8:14] # 提取 6bp UMI
        new_id = f'{id}_{umi}'
        o1.write(f'{new_id}\n{seq[17:]}\n+\n{qual[17:]}\n') # 切除 Index/UMI/GGG
        o2.write(f'{new_id}\n{seq2}\n+\n{qual2}\n')
"
 
# 2. BWA # (使用 2 线程)     # (使用 8 线程) 
#  比对 lscpu  Thread(s) per core:   1     Core(s) per socket:   2 
# 注意：如果你还没索引，请先执行 bwa index hg38.fa
bwa mem -t 2 ~/references/hg38/hg38.fa clean_R1.fq clean_R2.fq | \
samtools view -Sb - | \
samtools sort -o sorted.bam -

# 3. 建立索引
samtools index sorted.bam
echo "Step 2: 比对与排序完成，生成 sorted.bam。"