# HBD High-Fidelity RNA-seq Analysis

本项目专门用于处理 **Bulk RNA-seq** 中 UMI 饱和问题。通过结合 **6bp UMI** 与 **Read2 随机断裂位点坐标**，实现数字 PCR 级的分子定量。

## 技术亮点
- **Dual-Level Identity**: 结合 UMI 标签与基因组物理位置，消除 PCR 扩增偏好。
- **High Sensitivity**: 支持 30M+ Reads 的快速去重。
- **Visualized Fidelity**: 自动生成分子家族分布报告。

## 快速开始
1. 克隆仓库: `git clone ...`
2. 运行脚本: `./run_hbd.sh`