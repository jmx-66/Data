# 跨界跨类型共现网络 · WCI · SEM 可复现分析包

本仓库汇集一套完整的、可端到端复现的分析流程，用于研究**泥炭地土壤分解者系统**中
四类节点（微生物群落 / 微生物功能 / 碳库指标 / 酶化学计量）在**跨界跨类型共现网络**中的
结构性角色，并通过**加权贡献指数（WCI）**量化各节点类别的结构性关键节点富集，
最后用**分段结构方程模型（piecewise SEM）**解析环境 → 系统因子 → 网络终点的级联路径。

> 分析单元：4 组 × 6 样本 = 24 个样本（LL 低水位洼地 / LM 低水位高丘 / HL 高水位洼地 / HM 高水位高丘）。

---

## 1. 目录结构

```
crossdomain-network-repro/
├── README.md                    # 本文件
├── LICENSE                      # MIT
├── requirements.txt             # Python 依赖
├── run_all.py                   # 一键端到端复现
├── code/
│   ├── 01_build_network.py      # 跨界网络构建（门水平聚合 → Spearman → 关键种）
│   ├── 02_calculate_wci.py      # WCI 三分量计算（熵权 + Max-scaling）
│   └── 03_run_sem.py            # 级联三段 piecewise SEM
├── data/
│   ├── raw/                     # 原始输入（ASV / 环境 / 功能预测）
│   │   ├── 16s.csv  its.csv  18s.csv
│   │   ├── env_network.csv
│   │   ├── picrust2_sample_level1_relative.tsv
│   │   ├── funguild_sample_total_absolute.tsv
│   │   └── README.md            # 原始数据说明
│   ├── processed/
│   │   └── sem/                 # SEM 输入（样本级标准化表 + 上游功能表溯源）
│   └── reference/               # 归档（手稿所用）结果，供对照，不参与计算
├── results/
│   ├── tables/{network,wci,sem}/
│   └── figures/
└── docs/
    ├── methods.md               # 完整方法学（建网 / WCI / SEM）
    ├── data_dictionary.md       # 数据字典（逐文件逐列）
    ├── reproducibility.md       # 复现步骤与环境、与归档的差异
    └── decisions_log.md         # 中间决策记录（供验收）
```

---

## 2. 快速开始

```bash
# 1) 环境（Python ≥ 3.9）
pip install -r requirements.txt

# 2) 一键复现（网络 → WCI → SEM）
python run_all.py

# 或分步执行
python code/01_build_network.py     # 生成 results/tables/network/
python code/02_calculate_wci.py     # 生成 results/tables/wci/
python code/03_run_sem.py           # 生成 results/tables/sem/
```

**复现归档（手稿）口径的 WCI**：仓库保留了手稿使用的归档节点表 `data/reference/Nodes_*.csv`，
可直接用同一脚本复现手稿 WCI：

```bash
python code/02_calculate_wci.py --nodes-dir data/reference --out results/tables/wci_archived
```

---

## 3. 方法概览

| 环节 | 关键设定 |
|:--|:--|
| 节点类别 | Community（细菌/真菌/原生生物门）· Function（KEGG L1 / FUNGuild）· Carbon（碳库 24 项）· Enzyme（酶化学计量 10 项） |
| 门水平聚合 | ASV 按门求和 → 相对丰度 → 剔除平均相对丰度 < 0.1% 的罕见门 |
| 关联度量 | Spearman 秩相关 |
| 建边阈值 | p < 0.05（n=6 时等价 \|r\| ≥ 0.8286） |
| 模块划分 | greedy modularity（Newman–Girvan 模块度） |
| 关键种判定 | 介数中心性 Top 10% **或** Zi-Pi 模块角色（Connector / Module Hub / Network Hub）；**不使用 Neg 判据** |
| WCI 分量 | BWKI（全节点介数均值）+ DCI（全节点度均值）+ NKR（关键种富集率，Monte Carlo 10⁵） |
| WCI 权重 | 熵权法，对三分量**原始值**定权；归一化用 **Max-scaling** |
| SEM | 级联三段分段 OLS：环境(微地形/水位/pH) → 中间四系统 PC1 → 终点(复杂性/抵抗力/恢复力)，样本级 n=24 + 分组 n=12 |

详见 [docs/methods.md](docs/methods.md)。

---

## 4. 主要结果（本包从原始数据复现）

**网络拓扑**（0.1% 门水平过滤）：

| 组 | 节点 | 边 | 正边 | 负边 | 密度 | 平均度 | 模块度 |
|:--|--:|--:|--:|--:|--:|--:|--:|
| LL | 90 | 300 | 183 | 117 | 0.075 | 6.67 | 0.560 |
| LM | 92 | 445 | 250 | 195 | 0.106 | 9.67 | 0.391 |
| HL | 93 | 350 | 193 | 157 | 0.082 | 7.53 | 0.593 |
| HM | 88 | 290 | 143 | 147 | 0.076 | 6.59 | 0.517 |

**WCI 主导因子**：LL→微生物功能、LM→微生物功能、HL→微生物功能、HM→微生物群落。

**SEM**：见 `results/tables/sem/SEM_Cascade_Results.csv`。

> ⚠️ **与手稿归档结果的差异**：手稿 WCI 基于一套**历史归档节点表**（`data/reference/`，
> 门水平未做 0.1% 过滤、含 `Others_b`/`Unclassified` 合并节点），其主导因子为
> LL→酶、LM→碳、HL→碳、HM→功能。本包从原始数据按**当前文档化方法（0.1% 过滤）**
> 复现，网络规模与主导因子因此不同。**该差异的成因、影响与处置建议见
> [docs/reproducibility.md](docs/reproducibility.md) 与 [docs/decisions_log.md](docs/decisions_log.md)。**

---

## 5. 可复现性验证结论

- **WCI**：给定归档节点表，`02_calculate_wci.py` 与手稿结果**逐位一致**（最大绝对误差 7.8×10⁻¹⁶，主导因子完全相同）。
- **SEM**：`03_run_sem.py` 与归档级联结果**逐行一致**（104 行，β 与 R² 最大绝对误差 0，显著性判定一致率 100%）。
- **网络**：可由原始 ASV/环境/功能表端到端重建；但与归档节点表存在已知口径差异（见上）。

---

## 6. 许可与引用

代码以 MIT 许可发布（见 `LICENSE`）。数据使用请遵循原始数据提供方要求。