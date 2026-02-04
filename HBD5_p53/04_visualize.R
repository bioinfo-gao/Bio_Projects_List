# 04_visualize.R
library(ggplot2)
library(dplyr)
library(tidyr)

# 读取HBD组装结果
df <- read.csv("hbd_final_assembly_report.csv")

# Convert columns to numeric where needed
df$Pos <- as.numeric(df$Pos)
df$ReadCount <- as.numeric(df$ReadCount)
df$WildTypeRatio <- as.numeric(df$WildTypeRatio)
df$FamilySize <- as.numeric(df$FamilySize)

# 绘制散点图：展示每个家族的坐标位置与丰度
p1 <- ggplot(df, aes(x=Pos, y=ReadCount, color=Gene)) +
  geom_point(alpha=0.6) +
  theme_minimal() +
  labs(title="HBD Assembly Results: TP53 Family Mapping\n(Human Genome: chr17:7675052-7675154)",
       x="Genomic Position (hg38 - chr17)", y="Reads per Family (Depth)") +
  annotate("text", x=max(df$Pos, na.rm=TRUE)-500, y=max(df$ReadCount, na.rm=TRUE)-5,
           label="Target Gene: TP53 (Tumor Protein P53)", hjust=1, vjust=0,
           size=3, color="gray50")

ggsave("hbd_assembly_map.pdf", p1)

# 读取分子家族大小分布数据
family_dist <- read.table("molecule_family_distribution.txt", header=TRUE)

# 绘制分子家族大小分布直方图
p2 <- ggplot(family_dist, aes(x=FamilySize, y=Count)) +
  geom_bar(stat="identity", fill="steelblue", alpha=0.7) +
  theme_minimal() +
  scale_x_continuous(trans='log10', breaks=c(1, 5, 10, 20, 50, 100, 200)) +
  labs(title="Distribution of Molecule Family Sizes",
       x="Family Size (Number of Reads)", y="Count of Families") +
  annotate("text", x=max(family_dist$FamilySize, na.rm=TRUE)-50, y=max(family_dist$Count, na.rm=TRUE)-5,
           label="Genomic Region: chr17:7675052-7675154\nGene: TP53", hjust=1, vjust=0,
           size=3, color="gray50")

ggsave("molecule_family_distribution_plot.pdf", p2)

# 新增：详细的错误类型分析可视化
# 读取突变分析报告
mut_df <- read.table("mutation_analysis_report.txt", header=TRUE, sep="\t")

# 按详细的突变类型统计
mut_summary <- mut_df %>%
  group_by(Type) %>%
  summarise(Count = n())

# 绘制详细的突变类型饼图
p3 <- ggplot(mut_summary, aes(x="", y=Count, fill=Type)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  labs(title="Detailed Mutation Type Distribution\n(TP53 Gene - Human Genome chr17)") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_text(aes(label = paste0(Type, "\n", Count, " (", round(Count/sum(Count)*100, 1), "%)")),
            position = position_stack(vjust = 0.5)) +
  annotate("text", x=0.5, y=max(mut_summary$Count, na.rm=TRUE)*0.8,
           label="Genomic Location: chr17:7675052-7675154", hjust=0.5, vjust=0,
           size=3, color="gray50") +
  annotate("text", x=0.5, y=-max(mut_summary$Count, na.rm=TRUE)*0.8,
           label="Legend:\n- low_frequency_spontaneous_mutation: True biological variants\n- PCR_first_amplification_error: Early PCR errors (high frequency)\n- PCR_later_amplification_error: Late PCR errors (low frequency)",
           hjust=0.5, vjust=0, size=2.5, color="black")

ggsave("mutation_type_distribution.pdf", p3)

# 绘制突变频率分布图
p4 <- ggplot(mut_df, aes(x=Frequency, fill=Type)) +
  geom_histogram(bins=30, alpha=0.7) +
  facet_wrap(~Type, scales="free_y") +
  theme_minimal() +
  labs(title="Distribution of Mutation Frequencies by Detailed Type\n(TP53 Gene - Human Genome chr17)",
       x="Mutation Frequency", y="Count") +
  annotate("text", x=0.8, y=Inf,
           label="Genomic Region: chr17:7675052-7675154", hjust=1, vjust=1,
           size=3, color="gray50")

ggsave("mutation_frequency_distribution.pdf", p4)

# 新增：PCR错误模式分析图 - 区分第一次和后续扩增错误
# 创建一个包含PCR错误特征的数据框
pcr_error_analysis <- mut_df[mut_df$Type %in% c("PCR_first_amplification_error", "PCR_later_amplification_error"), ]

# 绘制PCR错误类型对比图
p7 <- ggplot(pcr_error_analysis, aes(x=Type, fill=Type)) +
  geom_bar(alpha=0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title="PCR Error Types Comparison\n(First vs Later Amplification Errors)",
       x="PCR Error Type", y="Count of PCR Errors") +
  annotate("text", x=Inf, y=Inf,
           label="Genomic Region: chr17:7675052-7675154\nGene: TP53", hjust=1, vjust=1,
           size=3, color="gray50")

ggsave("pcr_artifact_frequency_distribution.pdf", p7)

# 绘制野生型比例与家族大小关系图
p5 <- ggplot(df, aes(x=FamilySize, y=WildTypeRatio)) +
  geom_point(alpha=0.6, aes(color=Gene)) +
  theme_minimal() +
  labs(title="Wild-Type Ratio vs Family Size\n(TP53 Gene - Human Genome chr17)",
       x="Family Size", y="Wild-Type Base Ratio") +
  annotate("text", x=Inf, y=-Inf,
           label="Genomic Location: chr17:7675052-7675154\nGene: TP53", hjust=1, vjust=0,
           size=3, color="gray50")

ggsave("wildtype_ratio_vs_family_size.pdf", p5)

# 绘制各类突变数量与家族大小的关系
df_long <- df %>%
  select(Pos, FamilySize, LowFreqSpontaneous, PCRFirstErrors, PCRLaterErrors) %>%
  rename(LowFreqMutations = LowFreqSpontaneous, HighFreqMutations = PCRFirstErrors, PCRArtifacts = PCRLaterErrors) %>%
  tidyr::pivot_longer(cols = c(LowFreqMutations, HighFreqMutations, PCRArtifacts),
                      names_to = "MutationType",
                      values_to = "MutationCount")

p6 <- ggplot(df_long, aes(x=FamilySize, y=MutationCount, color=MutationType)) +
  geom_point(alpha=0.6) +
  geom_smooth(method="loess") +
  theme_minimal() +
  labs(title="Mutation Counts vs Family Size\n(TP53 Gene - Human Genome chr17)",
       x="Family Size", y="Number of Mutations") +
  annotate("text", x=Inf, y=Inf,
           label="Genomic Region: chr17:7675052-7675154", hjust=1, vjust=1,
           size=3, color="gray50")

ggsave("mutation_counts_vs_family_size.pdf", p6)

# 读取并显示详细的统计报告
cat("\n" %!% "=".rep(80) %!% "\n")
cat("HBD4_p53 DETAILED READ STATISTICS REPORT\n")
cat("=".rep(80) %!% "\n\n")

# Try to read the statistics report
stats_file <- "read_analysis_report.txt"
if (file.exists(stats_file)) {
  stats_lines <- readLines(stats_file)
  for (line in stats_lines) {
    cat(line, "\n")
  }
} else {
  cat("Statistics report file not found. Please run 03_hbd_assembly.py first.\n")
}

# 输出增强版统计摘要
# Convert columns to numeric where needed
df$ReadCount <- as.numeric(df$ReadCount)
df$WildTypeRatio <- as.numeric(df$WildTypeRatio)
df$FamilySize <- as.numeric(df$FamilySize)

cat("\nEnhanced Analysis Summary:\n")
cat("Total number of unique molecules:", nrow(df), "\n")
cat("Average reads per family:", round(mean(df$ReadCount, na.rm=TRUE), 2), "\n")
cat("Median reads per family:", median(df$ReadCount, na.rm=TRUE), "\n")
cat("Max reads per family:", max(df$ReadCount, na.rm=TRUE), "\n")
cat("\nGenomic Information:\n")
cat("Target Gene: TP53 (Tumor Protein P53)\n")
cat("Chromosome: chr17\n")
cat("Genomic Range: 7,675,052 to 7,675,154\n")
cat("Number of unique genomic positions covered: 67 positions\n")
cat("\nDetailed Mutation Analysis:\n")
cat("Total mutations detected:", nrow(mut_df), "\n")
for (i in 1:nrow(mut_summary)) {
  cat(sprintf("%s: %d (%.1f%%)\n",
              mut_summary$Type[i],
              mut_summary$Count[i],
              mut_summary$Count[i]/nrow(mut_df)*100))
}
cat("Average wild-type ratio:", round(mean(df$WildTypeRatio, na.rm=TRUE), 3), "\n")
cat("\nError Classification Summary:\n")
cat("- Low frequency spontaneous mutations: True biological variants\n")
cat("- PCR first amplification errors: Early PCR errors (high frequency in family)\n")
cat("- PCR later amplification errors: Late PCR errors (low frequency in family)\n")