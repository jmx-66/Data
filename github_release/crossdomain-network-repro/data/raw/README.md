# 原始数据说明（`data/raw/`）

本目录为全流程的**原始输入**。所有下游产物（`results/`）均可由这些文件端到端重建
（重建脚本：`code/01_build_network.py`）。逐列含义见 `docs/data_dictionary.md`。

样本命名统一为 `{组}{编号}`：组 ∈ `LL/LM/HL/HM`，编号 1–6，共 **4 组 × 6 = 24** 个样本。
分组语义：LL = 低水位·洼地，LM = 低水位·高丘，HL = 高水位·洼地，HM = 高水位·高丘。

| 文件 | 内容 | 规格 | 编码 |
|:--|:--|:--|:--|
| `16s.csv` | 细菌 16S ASV 丰度表（含 `Taxonomy`） | 26,127 ASV × 24 样本 | UTF-8 (BOM) |
| `its.csv` | 真菌 ITS ASV 丰度表（含 `Taxonomy`） | 16,009 ASV × 24 样本 | UTF-8 (BOM) |
| `18s.csv` | 原生生物 18S ASV 丰度表（含 `Taxonomy`） | 2,967 ASV × 24 样本 | UTF-8 (BOM) |
| `env_network.csv` | 环境/碳库/FTIR/酶化学计量表（样本为行） | 58 列 × 24 行 | UTF-8 (BOM) |
| `picrust2_sample_level1_relative.tsv` | 细菌功能预测（PICRUSt2，KEGG L1 相对丰度） | 6 通路 × 24 样本 | UTF-8 (BOM) |
| `funguild_sample_total_absolute.tsv` | 真菌功能预测（FUNGuild，绝对丰度） | 35 列（含 `Guild`） | UTF-8 (BOM) |

## 格式约定

1. **ASV 表**：第 1 列为 ASV 标识（`#ASV_ID`），其后 24 列为样本丰度（原始计数），最后一列为
   `Taxonomy`（形如 `k__;p__;c__;o__;f__;g__;s__`，第 2 段 `p__` 为**门**）。
2. **环境表**：第 1 列为 `Sample`，其后为各指标列；样本为行。
3. **功能表**：PICRUSt2 为 `level1` + 24 样本列；FUNGuild 含 `Guild` 列（分析取该列）。

## 数据来源与共享

- 本包随附的原始文件为**本研究实测数据**的整理结果，用于同行复现。
- ⚠️ 测序平台、引物、建库批次等**实验元信息**不在本仓库范围内；如需引用请以研究论文的
  "材料与方法"为准（本包未对这些元信息做任何推断）。
- 若涉及人类受试者或受保护物种数据，请遵循原始数据提供方的共享协议后再公开。

## 运行

```bash
python code/01_build_network.py            # 默认 0.1% 门水平过滤
python code/01_build_network.py --min-abund 0.0002   # 如需自定义过滤阈值
```