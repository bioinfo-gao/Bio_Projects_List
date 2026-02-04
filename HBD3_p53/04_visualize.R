# 04_visualize.R
library(ggplot2)
df <- read.csv("hbd_final_assembly_report.csv")

# 绘制散点图：展示每个家族的坐标位置与丰度
p <- ggplot(df, aes(x=Pos, y=ReadCount, color=Gene)) +
  geom_point(alpha=0.6) +
  theme_minimal() +
  labs(title="HBD Assembly Results: TP53 Family Mapping",
       x="Genomic Position (hg38)", y="Reads per Family (Depth)")

ggsave("hbd_assembly_map.pdf", p)