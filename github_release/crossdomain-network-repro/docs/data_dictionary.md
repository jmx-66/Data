# 数据字典（Data Dictionary）

所有表格为 **UTF-8（含 BOM）** 编码的 CSV/TSV，第一行为表头。
样本列名统一为 `{组}{编号}`（`LL/LM/HL/HM` × `1..6`）。

---

## 一、原始输入（`data/raw/`）

### 1. `16s.csv` / `its.csv` / `18s.csv` —— ASV 丰度表
26 列。

| 列 | 含义 |
|:--|:--|
| `#ASV_ID` | ASV 唯一标识 |
| `LL1 … HM6` | 24 个样本的 ASV 丰度（原始 reads/计数） |
| `Taxonomy` | 分类学字符串，格式 `k__;p__;c__;o__;f__;g__;s__`，第 2 段（`p__`）为**门** |

- 行数：16S 26,127；ITS 16,009；18S 2,967。

### 2. `env_network.csv` —— 环境 / 碳库 / 酶化学计量表
58 列 × 24 行。样本为行。

| 列 | 含义 |
|:--|:--|
| `Sample` | 样本名（`LL1…HM6`） |
| `Group` | 组别（LL/LM/HL/HM） |
| `pH, SWC` | 酸碱度、土壤含水量 |
| `TC, MBC, DOC, TN, DON, MBN` | 总碳、微生物量碳、可溶有机碳、总氮、可溶有机氮、微生物量氮 |
| `C/N` | 碳氮比 |
| `POM, MAOM` | 颗粒态 / 矿物结合态有机质 |
| `POM-C, POM-N, MAOM-C, MAOM-N` | POM/MAOM 的碳、氮 |
| `POC, MAOC, POC/MAOC` | 颗粒/矿物结合态有机碳及其比 |
| `POM-N/MAOM-N` | POM-N 与 MAOM-N 之比 |
| `HI` | 疏水性指数 |
| `IR_1029, IR_1630, IR_2900, IR_3400` | FTIR 特征波数吸收（无 P 前缀） |
| `MBC/TC, POM/TC, MAOM/TC, POM/MAOM` | 各比值指标 |
| `IR_Al_Aro, IR_Aro_Poly, IR_Al_Poly, IR_Aro_Al, IR_Poly_Aro, IR_Al_OH, IR_Poly_OH, IR_Aro_OH, IR_Poly_Al` | FTIR 官能团比值指标 |
| `AP, BG, LAP, CBH, PER, POX` | 酶活（酸性磷酸酶、β-葡萄糖苷酶、亮氨酸氨肽酶、纤维二糖水解酶、过氧化物酶、酚氧化酶） |
| `C_acq, N_acq, P_acq` | 碳/氮/磷获取酶（向量法分量） |
| `Hydrolytic, Oxidative, HO_Ratio` | 水解/氧化酶及比值 |
| `BG_CBH_Ratio` | BG/CBH 比 |
| `CUE` | 微生物碳利用效率 |
| `Vector_Angle, Vector_Length` | 酶化学计量向量角度/长度 |
| `Enz_CN_log, Enz_CP_log, Enz_NP_log` | 酶化学计量对数比 |

### 3. `picrust2_sample_level1_relative.tsv` —— 细菌功能（KEGG L1）
26 列（`level1` + 24 样本）。`level1` 为 KEGG 一级通路名，值为相对丰度。

### 4. `funguild_sample_total_absolute.tsv` —— 真菌功能（FUNGuild）
35 列：`#ASV_ID`、24 个样本列、`Taxonomy, Taxon, Taxon Level, Trophic Mode, Guild,
Growth Morphology, Trait, Confidence Ranking, Notes, Citation/Source`。
分析使用 `Guild` 列（前 15 类，其余并入 `Others_ff`）。

---

## 二、网络构建输出（`results/tables/network/`，由 01 生成）

### `Nodes_{group}.csv`（22 列）

| 列 | 含义 |
|:--|:--|
| `Id, Label` | 节点标识（门名 / 功能名 / 碳库名 / 酶名） |
| `Degree, Betweenness, Closeness` | 度、介数中心性、紧密度 |
| `Kingdom` | Bacteria / Fungi / Protist / Carbon / Enzyme |
| `System` | Community / Function / Carbon / Enzyme |
| `Category` | Bac_Phylum / Fun_Phylum / Pro_Phylum / Bac_Function / Fun_Function / Carbon_Pool / Enzyme_Stoi |
| `Module` | greedy modularity 模块编号 |
| `Neg_Degree, Pos_Degree, Neg_Ratio, Neg_Importance` | 负/正边计数及比值（仅供诊断，不参与关键种判定） |
| `Zi, Pi, ZiPi_Role` | 模块内/间连接度与角色 |
| `Btw_Rank` | 介数降序排名 |
| `Keystone_Btw` | 介数 Top 10% 标记 |
| `Keystone_ZiPi` | Zi-Pi 角色关键种标记 |
| `Keystone` | 最终关键种（Keystone/Others） |
| `Keystone_Type` | Both / Betweenness / ZiPi / Others |
| `Group` | 组别 |

### `Edges_{group}.csv`

`Source, Target, Correlation, P_value, Type, Weight, Source_Kingdom, Target_Kingdom,
Source_Category, Target_Category, Link_Category`。
`Type ∈ {Positive, Negative}`；`Weight = |Correlation|`；`Link_Category` 为类别首字母组合（如 `B-B`）。

### `Network_Topology_Summary.csv`

`Group, Nodes, Edges, Positive_Edges, Negative_Edges, Density, Avg_Degree,
Clustering_Coef, Modularity`。

---

## 三、WCI 输出（`results/tables/wci/`，由 02 生成）

### `WCI_Results.csv`（16 行 = 4 组 × 4 系统）

| 列 | 含义 |
|:--|:--|
| `Group, System` | 组别、节点类别 |
| `N_Nodes, N_Keystones` | 该系统节点数、关键种数 |
| `BWKI, DCI, NKR` | 三分量原始值 |
| `P_Value, Enrichment` | NKR 置换 p 值与富集判定 |
| `BWKI_norm, DCI_norm, NKR_norm` | Max-scaling 归一值 |
| `WCI` | 加权合成值 |
| `Rank` | 组内排名（1 = 该组主导因子） |

### `WCI_Weights.csv`

`Component, Entropy, Difference, Weight`。

---

## 四、SEM 输入（`data/processed/sem/`）

### `SEM_Cascade_Input.csv`（89 列 × 24 行）

样本级宽表。关键列：

| 列 | 含义 |
|:--|:--|
| `Sample, Group` | 样本、组别 |
| `pH, SWC` | 环境 |
| `TC(mg/g) … Enz_NP_log` 等 | 各原始碳库/酶/理化指标 |
| `Nodes, Edges, Avg_Degree, Complexity` | 网络拓扑与复杂性（PCA 合成） |
| `Overall_Stability, Robustness` | 网络稳定性（恢复力）/ 抗干扰（抵抗力） |
| `NonRandom_Z` | 网络非随机性 |
| `Carbon_Stability, Microbe_Stability, Function_Stability, Enzyme_Stability` | 四系统稳定性分量 |
| `WCI_Carb, WCI_Comm, WCI_Enz, WCI_Func` | 四系统 WCI |
| `Bac_Shannon, Fun_Shannon, Prot_Shannon` | 多样性 |
| `Bac_* / Fung_*` | KEGG L1 / FUNGuild 功能相对丰度 |
| `EMF` | 生态系统多功能性 |
| `Env_Carbon_PC1, Env_Enzyme_PC1, Env_Community_PC1, Env_Function_PC1` | **四系统因子 PC1（SEM 中间层）** |

### `SEM_n24_Input.csv`（13 列 × 24 行）

| 列 | 含义 |
|:--|:--|
| `Sample, Group` | 样本、组别 |
| `WaterLevel_binary` | 水位二值（1=高，0=低） |
| `Microtopo_binary` | 微地形二值（1=高丘，0=洼地） |
| `pH, SWC` | 环境 |
| `Carbon_PC1, Enzyme_PC1, Function_PC1, Community_PC1` | 四系统因子 PC1 |
| `NPI` | 网络性能指数 |
| `WaterLevel, Microtopo` | 水位/微地形的文字标签（High/Low、Hummock/Hollow） |

### 上游功能表（溯源，**不被 `03_run_sem.py` 读取**）

以下两表为 `Env_Function_PC1` 的**上游溯源文件**（功能维度的样本级汇总），
用于说明中间层功能因子从何而来；本包 SEM 脚本不读取它们。

| 文件 | 内容 |
|:--|:--|
| `bacterial_function_KO_level2_sample.csv` | 细菌 KEGG **L2** 通路样本级相对丰度（19 列 `Bac_*` + 样本列） |
| `fungal_function_FUNGuild_sample.csv` | 真菌 FUNGuild 各功能类群样本级丰度（12 列 `Fung_*` + 样本列） |

> 注：这两表的样本列为 `LL_1` 形式（带下划线），与主流程 `LL1` 口径不同；
> 因其仅为溯源展示、不参与计算，故保留原样，未做重命名。

---

## 五、SEM 输出（`results/tables/sem/`，由 03 生成）

### `SEM_Cascade_Results.csv`（104 行）

| 列 | 含义 |
|:--|:--|
| `Stage` | 环境→中间 / 中间→终点 |
| `Target` | 因变量（中间因子或终点指标） |
| `Source` | 自变量 |
| `Beta_std` | 标准化系数 β |
| `P` | 双尾 p 值 |
| `sig` | 显著性星号（** p<0.01，* p<0.05） |
| `R2` | 模型 R² |
| `n` | 样本量 |
| `Model` | 模型标签（整体n24 / 微地形=高丘(n=12) / …） |

### `SEM_Cascade_R2_Summary.csv`

`Model, n, M:{中间因子 R²}, E:{终点 R²}`。

---

## 六、参考（归档）数据（`data/reference/`）

手稿所用的**历史归档结果**，仅供对照，不参与本包计算：

- `Nodes_{group}.csv`（23 列，含 `Neg_Rank/Keystone_Neg/Keystone_Type`；门水平未做 0.1% 过滤）
- `Edges_{group}.csv`
- `Network_Topology_Summary.csv`
- `Keystone_Statistics_v6.csv`（归档关键种统计：`Group, Total_Nodes, Total_Keystones, Keystone_Pct, Btw_Keystones, Neg_Keystones, ZiPi_Keystones, …`；含 Neg 口径，见 `decisions_log.md` §D3）
- `WCI_v3_3comp_Final_Results.csv` / `WCI_v3_3comp_Final_Weights.csv`
- `SEM_Cascade_Archived_Results.csv`
- `NOTE.md`（本目录说明）

来源与差异详见 `data/reference/NOTE.md` 与 `docs/reproducibility.md`。