# 03_align.sh
REF="/home/gao/references/hg38/hg38.fa"

# -t 2 适配你的 2 线程 CPU
bwa mem -t 2 $REF clean_R1.fq clean_R2.fq | \
samtools view -Sb - | \
samtools sort -@ 2 -o sorted.bam -

samtools index sorted.bam
echo "Step 3 完成：比对结果已排序并建立索引。"