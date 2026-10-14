# 中间决策记录（Decisions Log）

> **用途**：本文件记录构建本开源复现包过程中，凡**存在多种可选口径、且选择会改变数值结果**的环节。
> 每条决策给出【背景】【证据】【决定】【影响】【验收状态】五要素，供老板最终验收。
> 凡证据不足或存在未决分歧者，一律标注为 ⚠️**待验收**，不做推测性结论。

---

## 决策索引

| 编号 | 主题 | 决定 | 验收状态 |
|:--|:--|:--|:--|
| D1 | WCI 的语义定位 | 结构性描述量，不等价于功能重要性 | ✅ 已定 |
| D2 | 节点表口径（归档 vs 文档化方法） | 双轨并存：`data/reference/` 归档 + `code/01` 现管道 | ⚠️ 待验收 |
| D3 | 关键种判据是否含 Neg | 不含 Neg | ✅ 已定 |
| D4 | WCI 分量构成 | 三分量 BWKI+DCI+NKR，剔除 RI | ✅ 已定 |
| D5 | 权重与归一化方法 | 原始值熵权 + Max-scaling | ✅ 已定 |
| D6 | 关键种判定逻辑 | 介数 Top10% **∪** Zi-Pi（或，非且） | ✅ 已定 |
| D7 | NKR 的 Monte Carlo 参数 | `N_MC=100000`、`seed=42`、`np.random.choice` | ✅ 已定 |
| D8 | 模块划分与 Zi-Pi 复算 | greedy modularity；WCI 复用节点表预计算角色 | ✅ 已定 |
| D9 | 样本命名一致性 | 统一 `HL1`（无下划线） | ✅ 已定 |
| D10 | 碳库/酶节点映射 | 碳库 24 项、酶 10 项 | ✅ 已定 |
| D11 | SEM 实现方式 | 分段标准化 OLS（piecewise），非协方差 SEM | ✅ 已定 |
| D12 | 网络拓扑图是否纳入 | 不纳入（原项目无可复现绘图脚本） | ✅ 已定 |
| D13 | 目录去重与文件归置 | 删除与 `data/reference/` 逐字节重复的 `data/processed/network/` | ✅ 已定 |

---

## D1 · WCI 的语义定位：结构性描述量

**背景**：WCI 由网络拓扑中心性合成，度量的是"某节点类别在整合网络中占据结构性关键位置的程度"。
但**结构中心性 ≠ 功能重要性**——占据关键拓扑位置的类群未必对生态系统功能有等量贡献。

**证据**：
- Yang Q, et al. *Environmental Microbiology Reports*, 2026, 18(2): e70245 —— 以丰度/结构中心性识别的关键种常缺乏对应的转录活性。
- 项目内实测：WCI 主导因子排序与以 EMF 为终点的 SEM 结果**不一致**（两者为相互正交的属性：前者是网络结构描述，后者是系统因子对功能的驱动）。

**决定**：WCI 仅作为**结构性关键节点富集的描述量**；"主导因子"定义为"该节点类别取得组内最高 WCI"。
正文/文档中不将 WCI 结果表述为功能性主导。

**影响**：`docs/methods.md` §2.1、`README.md` 方法概览均据此措辞；`results/tables/wci/WCI_Results.csv`
的 `Rank` 为组内结构性排序，不含功能含义。

**验收状态**：✅ 已定（与手稿摘要/结果"由某因子主导"的定义口径一致）。

---

## D2 · 节点表口径：归档（手稿）与文档化方法双轨 ⚠️

**背景**：本包从原始数据重建网络时发现，**手稿所用归档节点表无法由任何现存脚本逐字节再生**，
两套节点表口径不同，进而 WCI 数值与主导因子不同。这是本包最需要老板拍板的一项。

**证据（已逐条核验）**：

| 项 | 归档（手稿）`data/reference/Nodes_*.csv` | 本包 `code/01` 现管道 |
|:--|:--|:--|
| 节点表列数 | 23（含 `Neg_Rank`、`Keystone_Neg`） | 22 |
| 门水平过滤 | 无 0.1% 过滤 | `min_abund=0.001`（0.1%） |
| 合并节点 | 含 `Others_b`、`Unclassified` 合并节点 | 不含 |
| 细菌门数（LL） | 58 | 23 |
| 节点数（LL） | 135 | 90 |
| 关键种判据 | 介数 Top10% ∪ ZiPi ∪ **Neg** | 介数 Top10% ∪ ZiPi |
| WCI 主导因子 | LL→酶、LM→碳、HL→碳、HM→功能 | LL→功能、LM→功能、HL→功能、HM→群落 |

- 归档门的合并规则（`Others_b`/`Unclassified`）与阈值**在任何现存脚本中均无对应实现**；
  已用多阈值（含 0.02%、0.1% 等）实测，**无任何单一阈值可精确重现归档门集合**。
- 结论：归档节点表属"来源不可再生"的历史产物，本包以 `data/reference/` 原样保留。

**决定（本包采用的权宜口径）**：**双轨并存**——
1. `code/01_build_network.py` + `code/02_calculate_wci.py` 走**当前文档化方法**（0.1% 过滤、无 Neg），保证"方法—代码—结果"自洽可复现；
2. `data/reference/` 保留归档表，并支持 `python code/02_calculate_wci.py --nodes-dir data/reference` **逐位复现手稿 WCI**（实测最大绝对误差 7.8×10⁻¹⁶，主导因子完全一致）。

**影响**：README §4 与 `docs/reproducibility.md` §4 均显式标注该差异。

**待老板决策的处置方案**：

- **方案 A（投稿一致性优先）**：以归档表为准，WCI 用 `--nodes-dir data/reference` 复现手稿数字；
  方法学注明"门水平未做 0.1% 过滤、含 Others_b/Unclassified 合并、关键种含 Neg"。
- **方案 B（方法严谨性优先）**：以现管道（0.1% 过滤、无 Neg）重跑网络与 WCI，方法—代码完全自洽，
  但主导因子与手稿现值不同，需同步修订正文。
- **方案 C**：明确区分"结构性描述量（WCI）"与"功能性结论"，接受两者主导因子不同（与 D1 一致）。

**验收状态**：⚠️ **待老板确认 A / B / C**。

---

## D3 · 关键种判据：不使用 Neg（负连接）

**背景**：归档关键种标签为"介数 Top10% ∪ Neg Top10% ∪ Zi-Pi"三者并集，
即把"负连接占比高"的节点也计为关键种。

**证据**：
- 归档表中**仅 4 个节点**靠 Neg 单项入选：LM `F_F_Plant Pathogen`、HL `MBNT15`、HM `WOR-1`、HM `MBC`。
- 剔除 Neg 后重算归档 WCI，**主导因子不变**（仍为 LL→酶、LM→碳、HL→碳、HM→功能）。

**决定**：本包关键种判定**不使用 Neg 判据**（符合项目个性化约束）。
`Neg_Degree/Pos_Degree/Neg_Ratio/Neg_Importance` 列**仅作诊断保留**，不参与判定。

**影响**：`code/01_build_network.py::identify_keystone()` 只做"介数 ∪ Zi-Pi"；
归档表中 `Keystone_Neg` 列在现管道中不再生成。

**验收状态**：✅ 已定。

---

## D4 · WCI 分量构成：剔除 RI（RobustnessImpact）

**背景**：早期 WCI 曾包含第 4 个分量 RI（RobustnessImpact，网络鲁棒性影响）。

**证据**：RI 主要由**节点占比**驱动（节点数越少 → RI 越高），而非节点类别的实际生态重要性，
会人为抬高节点数少的类别（Enzyme/Function）得分。

**决定**：最终 WCI 采用**三分量** `BWKI + DCI + NKR`，**剔除 RI**。

**影响**：`code/02_calculate_wci.py` 的 `COMPONENTS = ['BWKI','DCI','NKR']`；
`docs/methods.md` §2.2 据此描述。

**验收状态**：✅ 已定。

---

## D5 · 权重与归一化：原始值熵权 + Max-scaling

**背景**：需在"分量定权"与"量纲统一"两步之间确定算子。

**证据**：
- 熵权基于 `p_ij = x_ij / Σ_i x_ij`（比例结构）。**Max-scaling（x/max）为等比缩放，不改变 `p_ij`**，
  故"先定权、后归一化"两者兼容；
- **Min-Max 含平移**，会破坏比例结构，且当组间数值高度接近时（如某组多指标趋近 0）会产生畸变。

**决定**：
1. **熵权对三分量"原始值"计算**（先定权）；
2. 归一化用 **Max-scaling**（仅统一量纲、保比例）；
3. `WCI = Σ_j w_j · x_norm,j`。

**影响**：`code/02_calculate_wci.py::entropy_weights()` 直接吃原始值矩阵；
`WCI_Weights.csv` 中现管道权重为 BWKI 0.0547 / DCI 0.1101 / NKR 0.8352；
归档口径为 0.2175 / 0.2370 / 0.5454（节点表不同所致，见 D2）。

**验收状态**：✅ 已定。

---

## D6 · 关键种判定逻辑：并集（∪）而非交集（∩）

**背景**：曾试验"双条件"（介数 Top10% **且** degree Top10%）判定关键种。

**证据**：双条件导致关键种极度稀疏（0–2 个/系统），使 NKR 失去统计基础，
并出现 Function 主导的假性结果。

**决定**：采用**并集**：`Keystone = (介数 Top10%) ∪ (Zi-Pi ∈ {Connector, Module Hub, Network Hub})`。

**影响**：现管道各组关键种数 LL 11、LM 11、HL 10、HM 14（占比 10.8%–15.9%），量级合理。

**验收状态**：✅ 已定。

---

## D7 · NKR 的 Monte Carlo 参数

**背景**：NKR = 观测关键种率 / 随机期望关键种率，需明确置换规模与随机数口径。

**证据**：归档 WCI 在 `N_MC=100000`、`seed=42`、`np.random.choice` 口径下**逐位复现**（见 `reproducibility.md` §3.1）。

**决定**：`N_MC = 100000`、`SEED = 42`；置换从**"其他系统"节点**中随机抽取与该系统等量的节点；
`p = (置换中 ≥ 观测的次数 + 1) / (N_MC + 1)`。

**影响**：`code/02_calculate_wci.py` 默认参数；固定种子后结果完全可复现。

**验收状态**：✅ 已定。

---

## D8 · 模块划分与 Zi-Pi 复算

**背景**：Louvain/模块划分算法具随机性，不同脚本、不同图构建顺序会导致模块归属与 Zi-Pi 角色漂移
（即便同一随机种子亦可能不同）。

**决定**：
1. `code/01_build_network.py` 使用 `networkx.algorithms.community.greedy_modularity_communities`
   （Newman–Girvan 模块度，确定性）；
2. `code/02_calculate_wci.py` **直接复用节点表中预计算的 `ZiPi_Role` 与 `Betweenness`**，
   **不重算** Louvain 社区划分。

**影响**：WCI 的关键种集合与归档对齐；避免二次划分引入的漂移。

**验收状态**：✅ 已定。

---

## D9 · 样本命名一致性

**背景**：历史归档节点表生成时若用 `HL_1`（带下划线）与环境表 `HL1` 不一致，会导致**样本错位**，
使某些在组内相对丰度为 0 的节点被误判。

**决定**：本包统一使用 **`HL1`（无下划线）** 口径；
`code/01_build_network.py` 中 `get_group_samples()` 同时兼容两种命名，但内部一律归一到无下划线。

**验收状态**：✅ 已定。

---

## D10 · 碳库/酶节点映射

**背景**：需固定 Carbon/Enzyme 两类节点与 `env_network.csv` 列的对应关系。

**决定**（`code/01_build_network.py` 中 `CARBON_NODE_MAP` / `ENZYME_NODE_MAP`）：
- **碳库 24 项**：`TC, MBC, DOC, DON, MBN, C/N, POM, MAOM, POM-C, POM-N, MAOM-C, MAOM-N,
  POM-N/MAOM-N, HI, IR_1029, IR_1630, IR_2900, IR_3400, IR_Al_Aro, IR_Aro_Poly,
  IR_Al_Poly, IR_Aro_OH, CUE, MBC/TC`；
- **酶 10 项**：`C_acq, BG_CBH_Ratio, N_acq, P_acq, Oxidative, Hydrolytic,
  Vector_Angle, Vector_Length, Enz_CP_log, Enz_NP_log`；
- 各组中**零方差**的列自动剔除（故 Function 节点数随组在 15–20 间波动）。

**影响**：`docs/methods.md` §1.1、`README.md` §3 已同步为"碳库 24 项"。

**验收状态**：✅ 已定。

---

## D11 · SEM 实现方式：分段标准化 OLS（piecewise）

**背景**：级联路径"环境 → 中间四因子 → 终点三指标"可用协方差 SEM 或分段 OLS 实现。

**决定**：采用 **piecewise SEM（分段标准化 OLS）**：
- 段 1：环境变量 → 每个中间因子（各自一个回归）；
- 段 2：四个中间因子 → 每个终点指标；
- 报告标准化 β、双尾 t 检验 p、R²、n，环境段给出 VIF。

**证据**：`code/03_run_sem.py` 输出与归档 `SEM_Cascade_Archived_Results.csv`
**104 行逐行一致**（β、p、R² 最大绝对误差 0，显著性判定一致率 100%）。

**验收状态**：✅ 已定。

---

## D12 · 网络拓扑图不纳入本包

**背景**：原项目中网络可视化图（spring layout 等）**无可复现绘图脚本**；
`Nodes_*.csv` 亦不含节点坐标列。

**决定**：本包**不提供网络拓扑图绘图脚本**，`results/figures/` 留空（可由使用者自绘）。
此限制已在 `docs/reproducibility.md` §5 显式声明。

**验收状态**：✅ 已定（如老板需要出图，可另行提出，需先确认布局与配色口径）。

---

## D13 · 目录去重与文件归置

**背景**：整理仓库时发现 `data/processed/network/` 下存放的 `Nodes_*.csv`、`Edges_*.csv`
与 `data/reference/` 下同名文件**逐字节相同**（MD5 一致），实为归档网络的重复副本；
该目录还含一个未被任何脚本引用的 `Keystone_Statistics_v6.csv`。

**证据**：
- `Nodes_LL.csv` / `Nodes_HL.csv` / `Edges_LL.csv` / `Edges_HL.csv` 在 `data/processed/network/`
  与 `data/reference/` 下 **MD5 完全相同**；
- 无任何脚本或文档以 `data/processed/network/` 为输入（`02` 默认读 `results/tables/network/`）。

**决定**：
1. 将 `Keystone_Statistics_v6.csv` **移入 `data/reference/`**（作为归档证据，含 Neg 口径）；
2. **删除** `data/processed/network/` 目录（与 `data/reference/` 重复）；
3. `data/processed/` 仅保留 `sem/`（SEM 输入 + 上游功能表溯源）；
4. `results/tables/network/` 仍为 `code/01` 的默认输出位置。

**影响**：`README.md` 目录树、`docs/data_dictionary.md` §四/§六、`data/reference/NOTE.md` 已同步。

**验收状态**：✅ 已定（若老板希望保留中间产物副本，可恢复，但会与 `data/reference/` 重复）。

---

## 附：待老板验收清单

| # | 事项 | 需要老板决定/确认 |
|:--|:--|:--|
| 1 | **D2 节点表口径** | 在方案 A / B / C 中选择（影响 WCI 主导因子与正文一致性） |
| 2 | 是否将 `results/tables/wci_archived/`（由 `--nodes-dir data/reference` 生成）一并纳入仓库作为对照 | 默认不纳入，仅保留 `data/reference/` 原始归档 |
| 3 | `docs/methods.md` 参考文献的期刊格式 | 需按目标期刊调整（作者–年份制） |
| 4 | 是否需要补充网络拓扑可视化脚本（D12） | 默认不提供 |

> 本文件所有数值与差异均来自项目内已核验的脚本输出与归档文件；无证据支持的推断一律未写入。