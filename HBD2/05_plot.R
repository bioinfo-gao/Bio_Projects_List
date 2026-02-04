# 05_plot.R
library(ggplot2)
df <- read.csv("hbd_final_results.csv")

# 绘制分子家族大小分布图
p <- ggplot(df, aes(x=ReadCount)) +
  geom_histogram(binwidth=1, fill="darkblue", color="white") +
  scale_y_log10() +
  labs(title="HBD Molecular Redundancy Analysis",
       subtitle="Identifying PCR Duplicates from 2-Thread Linux Simulation",
       x="Reads per Unique Molecule", y="Frequency (Log10)") +
  theme_minimal()

ggsave("hbd_fidelity_report.pdf", p)
cat("Step 5 完成：分析报表已生成至 hbd_fidelity_report.pdf\n")