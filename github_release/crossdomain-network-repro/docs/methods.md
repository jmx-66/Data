# 方法学（Materials & Methods）

本文件逐环节记录跨界跨类型共现网络构建、WCI 计算与 SEM 建模的**数据来源、方法、计算过程与代码实现**。
所有设定均取自项目内已核实的生产脚本与归档结果，未做推测。

---

## 0. 数据来源总览

| 数据 | 文件（`data/raw/`） | 内容 | 样本 |
|:--|:--|:--|:--|
| 细菌群落 | `16s.csv` | 16S ASV 丰度表（含 Taxonomy） | 24 |
| 真菌群落 | `its.csv` | ITS ASV 丰度表（含 Taxonomy） | 24 |
| 原生生物群落 | `18s.csv` | 18S ASV 丰度表（含 Taxonomy） | 24 |
| 环境/碳库/酶 | `env_network.csv` | 58 列：理化、碳库、FTIR、酶化学计量 | 24 |
| 细菌功能 | `picrust2_sample_level1_relative.tsv` | PICRUSt2 KEGG L1 相对丰度 | 24 |
| 真菌功能 | `funguild_sample_total_absolute.tsv` | FUNGuild 功能预测（绝对丰度） | 24 |

样本命名：`{组}{编号}`，组为 `LL/LM/HL/HM`，编号 1–6；共 4 组 × 6 = 24 样本。
分组语义：LL = 低水位·洼地，LM = 低水位·高丘，HL = 高水位·洼地，HM = 高水位·高丘。

> **样本命名一致性**：ASV 表列名形如 `HL1`（无下划线），环境表样本列同样为 `HL1`；
> 二者必须严格一致，否则会样本错位（历史教训，见 `decisions_log.md`）。

---

## 1. 跨界跨类型共现网络构建（`code/01_build_network.py`）

### 1.1 节点定义（四类系统）

| System | 节点来源 | 聚合层级 |
|:--|:--|:--|
| **Community** | 16S / ITS / 18S | 门水平（phylum），分属 Bacteria / Fungi / Protist |
| **Function** | PICRUSt2 KEGG L1、FUNGuild | KEGG L1 通路；FUNGuild Guild（前 15 类，其余并入 `Others_ff`） |
| **Carbon** | `env_network.csv` | 24 个碳库指标（定义见下） |
| **Enzyme** | `env_network.csv` | 10 个酶化学计量指标 |

**碳库节点（24）**：`TC, MBC, DOC, DON, MBN, C/N, POM, MAOM, POM-C, POM-N, MAOM-C, MAOM-N, POM-N/MAOM-N, HI, IR_1029, IR_1630, IR_2900, IR_3400, IR_Al_Aro, IR_Aro_Poly, IR_Al_Poly, IR_Aro_OH, CUE, MBC/TC` 中实际存在且非恒定的列。

**酶节点（10）**：`C_acq(BG), BG_CBH_Ratio(CBH), N_acq(LAP), P_acq(AP), Oxidative(POX), Hydrolytic(PER), Vector_Angle, Vector_Length, Enz_CP_log, Enz_NP_log`。

### 1.2 门水平聚合与丰度过滤

1. 解析 ASV 的 Taxonomy 字符串，取第 2 段为门（`rank=1`），去除 `d__/p__` 等前缀；
   无法解析者记为 `Unclassified` 并**剔除**。
2. 按门对 ASV 丰度**求和**。
3. 转为**相对丰度**（每样本内归一）。
4. 计算每个门在样本间的**平均相对丰度**，剔除 `< 0.1%`（`min_abund=0.001`）的罕见门。

过滤效果（本包实测，24 样本）：

| 界 | 过滤前门数 | 0.1% 过滤后 |
|:--|--:|--:|
| 细菌 | 69 | 23 |
| 真菌 | 21 | 10 |
| 原生生物 | 8 | 6 |

### 1.3 建边

- 对**每组 6 个样本**单独建网。
- 节点对之间计算 **Spearman 秩相关** `r` 及 `p`。
- 保留 `p < 0.05` 的边（n=6 时对应 `|r| ≥ 0.8286`）。
- 边属性：相关系数 `Correlation`、`P_value`、正负 `Type`（Positive/Negative）、`Weight=|r|`。

### 1.4 拓扑指标与模块划分

- **度 Degree**、**介数中心性 Betweenness**、**紧密度 Closeness**（NetworkX）。
- **模块划分**：`networkx.algorithms.community.greedy_modularity_communities`（Newman–Girvan 模块度）。
- **Zi–Pi 角色**：
  - Zi = (k_in − mean(k_in, module)) / sd(k_in, module)
  - Pi = k_is / k（模块间连接比例）
  - `Zi>2.5 & Pi>0.62` → Network Hub；`Zi>2.5` → Module Hub；`Pi>0.62` → Connector；其余 Peripheral。

### 1.5 关键种判定

```
Keystone = (介数中心性 Top 10%)  ∪  (Zi-Pi 角色 ∈ {Connector, Module Hub, Network Hub})
```

- **明确不使用 Neg（负连接）判据**（本项目个性化决策，见 `decisions_log.md` §D3）。
- 输出 `Keystone`（Keystone/Others）、`Keystone_Type`（Both/Betweenness/ZiPi/Others）。

### 1.6 输出

- `results/tables/network/Nodes_{group}.csv`：节点表（22 列，含 `Group`、拓扑指标、Zi-Pi、Keystone 标签）。
- `results/tables/network/Edges_{group}.csv`：边表（含跨类别标签 `Link_Category`）。
- `results/tables/network/Network_Topology_Summary.csv`：四组拓扑汇总。

---

## 2. 加权贡献指数 WCI（`code/02_calculate_wci.py`）

### 2.1 设计目标

量化**每个节点类别（System）在整合网络中结构性关键节点富集的相对强度**，
用于判定各组的"结构性主导因子"。

> 术语约定：**主导 = 该节点类别取得最高 WCI**，即网络内结构性关键节点富集最强。
> 这是**结构性描述量**，不等价于功能重要性（见 `decisions_log.md` §D1）。

### 2.2 三个分量（16 格 = 4 组 × 4 系统）

| 分量 | 定义 | 计算 |
|:--|:--|:--|
| **BWKI** | Betweenness-Weighted Keystone Intensity | 该系统**全部节点**介数中心性之和 / 节点数（全节点均值） |
| **DCI** | Degree-Weighted Contribution Intensity | 该系统**全部节点**度之和 / 节点数（全节点均值） |
| **NKR** | Node-type Normalized Keystone Rate | 观测关键种率 / 随机期望关键种率 |

**NKR 的 Monte Carlo 置换**：对该系统，从其**他系统的节点**中随机抽取等量节点，
重复 `N_MC = 100,000` 次，得期望关键种率；`NKR = 观测率 / 期望率`。
`p = (置换中 ≥ 观测的次数 + 1)/(N_MC + 1)`；据此判 `显著富集 / 随机水平 / 显著缺失`。

> **RNG 口径**：使用 `np.random.seed(42)` + `np.random.choice`（与原脚本 pandas `.sample` 随机序列等价），
> 保证与归档结果逐位一致。

### 2.3 权重与合成

```
1) 熵权（对三分量"原始值"计算，先定权）：
     p_ij = x_ij / Σ_i x_ij
     e_j  = −(1/ln n) Σ_i p_ij ln p_ij        （n = 16）
     d_j  = 1 − e_j
     w_j  = d_j / Σ_j d_j
2) 归一化（仅统一量纲）：Max-scaling  x_norm = x / max(x)
3) 合成：WCI = Σ_j w_j · x_norm,j
```

**为什么"先定权、后归一化"**：Max-scaling 是等比缩放，不改变熵权的比例结构（`p_ij` 不变）；
而 Min-Max 含平移会破坏比例，故**不用于定权**。

本包实测熵权（`results/tables/wci/WCI_Weights.csv`）：

| 分量 | 熵 e | 差异 d | 权重 w |
|:--|--:|--:|--:|
| BWKI | 0.9948 | 0.0052 | 0.0547 |
| DCI | 0.9896 | 0.0104 | 0.1101 |
| NKR | 0.9207 | 0.0793 | 0.8352 |

> 归档（手稿）口径权重为 BWKI 0.2175 / DCI 0.2370 / NKR 0.5454（节点表不同所致，见 `reproducibility.md`）。

### 2.4 输出

- `results/tables/wci/WCI_Results.csv`：16 格 WCI、分量、组内 `Rank`、`Enrichment`。
- `results/tables/wci/WCI_Weights.csv`：熵权。

---

## 3. 级联三段结构方程模型（`code/03_run_sem.py`）

### 3.1 模型结构

```
环境层 (微地形 / 水位 / pH)
        │
        ▼
中间层 4 系统因子 PC1：碳 Env_Carbon_PC1 · 酶 Env_Enzyme_PC1
                       · 群落 Env_Community_PC1 · 功能 Env_Function_PC1
        │
        ▼
终点层 3 指标：复杂性 Complexity · 抵抗力 Robustness · 恢复力 Overall_Stability
```

### 3.2 方法

- **piecewise SEM**：分两段做**标准化（z-score）多元线性回归**。
  - 段 1：环境变量 → 每个中间因子（每个中间因子一个回归）。
  - 段 2：四个中间因子 → 每个终点指标。
- 报告**标准化系数 β**、双尾 t 检验 `p`、`R²`、`n`；环境段给出 **VIF**（多重共线性）。
- **样本级 n = 24**。

### 3.3 分组模型（5 个）

| 模型 | 样本 | 环境变量 |
|:--|:--|:--|
| 整体 | n=24 | pH、微地形二值、水位二值 |
| 微地形=高丘 | n=12 | pH、水位二值 |
| 微地形=洼地 | n=12 | pH、水位二值 |
| 水位=高 | n=12 | pH、微地形二值 |
| 水位=低 | n=12 | pH、微地形二值 |

### 3.4 输入数据

- `data/processed/sem/SEM_Cascade_Input.csv`：样本级 n=24，含 pH、Complexity、Robustness、
  Overall_Stability 与四个 `Env_*_PC1`（及各上游原始列）。
- `data/processed/sem/SEM_n24_Input.csv`：含 `WaterLevel_binary`、`Microtopo_binary` 等分组指示。

> `Env_*_PC1` 为各系统内多指标经 PCA 提取的 PC1 复合量；`Complexity/Robustness/Overall_Stability`
> 为网络拓扑与鲁棒性分析产出。其逐列来源见 `data_dictionary.md`。

### 3.5 输出

- `results/tables/sem/SEM_Cascade_Results.csv`：全部 5 模型的路径系数（104 行）。
- `results/tables/sem/SEM_Cascade_R2_Summary.csv`：各模型 R² 汇总。

---

## 4. 参考文献（方法学依据）

- Berry D, Widder S. Deciphering microbial interactions and detecting keystone species with co-occurrence networks. *Frontiers in Microbiology*, 2014, 5: 219.
- Banerjee S, Schlaeppi K, van der Heijden MGA. Keystone taxa as drivers of microbiome structure and functioning. *Nature Reviews Microbiology*, 2018, 16: 567–576.
- Faust K, Raes J. Microbial interactions: from networks to models. *Nature Reviews Microbiology*, 2012, 10: 538–550.
- Barberán A, Bates ST, Casamayor EO, Fierer N. Using network analysis to explore co-occurrence patterns in soil microbial communities. *ISME Journal*, 2012, 6: 343–351.
- de Vries FT, et al. Soil bacterial networks are less stable under drought than fungal networks. *Nature Communications*, 2018, 9: 3033.
- Kitson E, Bell NGB. The response of microbial communities to peatland drainage. *Frontiers in Microbiology*, 2020.
- Yang Q, Aghdam R, Tran PQ, Anantharaman K. Activity-informed network analysis reveals keystone microbes shaping freshwater ecosystem function. *Environmental Microbiology Reports*, 2026, 18(2): e70245.
- Cinelli M, et al. *Environment Systems and Decisions*, 2021.（熵权/多准则合成依据）
- Roszkowska E, Wachowicz T. *Entropy*, 2024, 26: 365；Roszkowska E. *Entropy*, 2026, 28: 114.（熵权定权与归一化依据）
- 跨界建网可行性实证：细菌–真菌跨界网络（*Frontiers in Microbiology*, fmicb-15-1289022）；多营养级网络（*Environmental Microbiome*, s40793-025-00788-1）。

> 注：以上条目取自项目内已核实的参考文献列表；引用格式（作者–年份）需按目标期刊调整。