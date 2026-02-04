# 04_visualize.R
library(ggplot2)
library(dplyr)

# 读取HBD组装结果
df <- read.csv("hbd_final_assembly_report.csv")

# 绘制散点图：展示每个家族的坐标位置与丰度
p1 <- ggplot(df, aes(x=Pos, y=ReadCount, color=Gene)) +
  geom_point(alpha=0.6) +
  theme_minimal() +
  labs(title="HBD Assembly Results: TP53 Family Mapping",
       x="Genomic Position (hg38)", y="Reads per Family (Depth)")

ggsave("hbd_assembly_map.pdf", p1)

# 读取分子家族大小分布数据
family_dist <- read.table("molecule_family_distribution.txt", header=TRUE)

# 绘制分子家族大小分布直方图
p2 <- ggplot(family_dist, aes(x=FamilySize, y=Count)) +
  geom_bar(stat="identity", fill="steelblue", alpha=0.7) +
  theme_minimal() +
  scale_x_continuous(trans='log10', breaks=c(1, 5, 10, 20, 50, 100, 200)) +
  labs(title="Distribution of Molecule Family Sizes",
       x="Family Size (Number of Reads)", y="Count of Families")

ggsave("molecule_family_distribution_plot.pdf", p2)

# 输出统计摘要
cat("Analysis Summary:\n")
cat("Total number of unique molecules:", nrow(df), "\n")
cat("Average reads per family:", round(mean(df$ReadCount), 2), "\n")
cat("Median reads per family:", median(df$ReadCount), "\n")
cat("Max reads per family:", max(df$ReadCount), "\n")