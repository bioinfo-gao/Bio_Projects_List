# 文件名: 04_plot_results.R
library(ggplot2)

data <- read.csv("counts.csv")

# 统计每个原始分子的副本数（PCR 重复情况）
p <- ggplot(data, aes(x=count)) +
  geom_histogram(binwidth=1, fill="steelblue", color="white") +
  scale_y_log10() +
  labs(title="HBD Molecular Redundancy",
       x="Reads per Original Molecule (PCR Copies)",
       y="Number of Molecules (Log10 Scale)") +
  theme_minimal()

ggsave("hbd_report.pdf", p)
cat("Step 4: 可视化报表 hbd_report.pdf 已生成。\n")