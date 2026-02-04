# 02_align.sh
REF="/home/gao/references/hg38/hg38.fa"

# 提取 UMI 并清洗
python3 -c "
with open('raw_R1.fastq') as f, open('clean_R1.fq', 'w') as o1, open('clean_R2.fq', 'w') as o2:
    f2 = open('raw_R2.fastq')
    for line in f:
        id = line.strip(); seq = next(f).strip(); next(f); qual = next(f).strip()
        seq2 = next(f2); next(f2); next(f2); qual2 = next(f2)
        umi = seq[8:14]; new_id = f'{id}_{umi}'
        o1.write(f'{new_id}\n{seq[17:]}\n+\n{qual[17:]}\n')
        o2.write(f'{id}_{umi}\n{seq2}')
"

# BWA 比对 (2线程)
bwa mem -t 2 $REF clean_R1.fq clean_R2.fq | samtools view -Sb - | samtools sort -@ 2 -o sorted.bam -
samtools index sorted.bam