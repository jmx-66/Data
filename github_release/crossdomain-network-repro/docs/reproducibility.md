# 复现说明（Reproducibility）

## 1. 环境

- Python ≥ 3.9（实测 3.10 / 3.11 均可）
- 依赖见 `requirements.txt`；实测版本：numpy 2.2.6、pandas 2.3.3、scipy 1.15.3、
  networkx 3.4.2、scikit-learn 1.7.2。

```bash
pip install -r requirements.txt
```

## 2. 复现步骤

```bash
python run_all.py
```

或分步：

| 步骤 | 命令 | 输入 | 输出 |
|:--|:--|:--|:--|
| 网络构建 | `python code/01_build_network.py` | `data/raw/*` | `results/tables/network/*` |
| WCI | `python code/02_calculate_wci.py` | `results/tables/network/Nodes_*.csv` | `results/tables/wci/*` |
| SEM | `python code/03_run_sem.py` | `data/processed/sem/*.csv` | `results/tables/sem/*` |

**随机性**：WCI 的 NKR 使用 `seed=42`、`N_MC=100000`；固定种子后结果可完全复现。

## 3. 可复现性验证结论（本包实测）

### 3.1 WCI —— 与归档逐位一致 ✅

用归档节点表运行同一 WCI 脚本：

```bash
python code/02_calculate_wci.py --nodes-dir data/reference --out results/tables/wci_archived
```

与 `data/reference/WCI_v3_3comp_Final_Results.csv` 比对：

| 分量 | 最大绝对误差 |
|:--|:--|
| BWKI | 0 |
| DCI | 0 |
| NKR | 0（≤4.4×10⁻¹⁶） |
| WCI | 7.8×10⁻¹⁶ |

主导因子完全一致：**LL→酶、LM→碳、HL→碳、HM→功能**。

### 3.2 SEM —— 与归档逐行一致 ✅

`03_run_sem.py` 输出与 `data/reference/SEM_Cascade_Archived_Results.csv` 比对：
104 行全部匹配，`Beta_std` 与 `R2` 最大绝对误差 **0**，显著性判定一致率 **100%**。

### 3.3 网络构建 —— 可端到端重建，但与归档存在口径差异 ⚠️

`01_build_network.py` 可由 `data/raw/` 端到端重建网络。但**归档节点表无法由本包（或任何现存脚本）逐字节再生**，
成因与影响如下。

---

## 4. 与归档（手稿）网络的口径差异

### 4.1 现象

| 项 | 归档（手稿） | 本包（文档化方法） |
|:--|:--|:--|
| 细菌门数（LL） | 58 | 23 |
| 节点数（LL） | 135 | 90 |
| 节点表列数 | 23（含 Neg 列） | 22 |
| 门水平过滤 | 无 0.1% 过滤，含 `Others_b`/`Unclassified` 合并节点 | 0.1% 平均相对丰度过滤 |
| 关键种判据 | 介数 Top10% ∪ ZiPi ∪ **Neg** | 介数 Top10% ∪ ZiPi（**无 Neg**） |
| WCI 主导因子 | LL→酶、LM→碳、HL→碳、HM→功能 | LL→功能、LM→功能、HL→功能、HM→群落 |

### 4.2 成因（已逐条核验）

1. **门水平过滤阈值不同**：本包按项目文档化方法使用 `min_abund=0.001`（0.1%）过滤，
   细菌门 69→23、真菌 21→10、原生生物 8→6。
   归档表保留了更多罕见门（细菌 58），并含 `Others_b`、`Unclassified` 两个**合并节点**，
   该合并规则与阈值在任何现存脚本中均无对应实现（已用多阈值实测，无任何阈值可精确重现归档门集合）。
2. **Neg 判据**：归档 `Keystone` 列 = 介数 Top10% ∪ Neg Top10% ∪ ZiPi 三者并集；
   其中 4 个节点仅凭 Neg 入选（LM `F_F_Plant Pathogen`、HL `MBNT15`、HM `WOR-1`、HM `MBC`）。
   本包按文档化方法**不含 Neg**。
3. **样本命名**：归档节点表生成时若用 `HL_1`（带下划线）与环境表 `HL1` 不一致，会导致样本错位；
   本包统一使用 `HL1` 无下划线口径。

### 4.3 影响

- **WCI 数值与主导因子随节点表变化**：见上表。
- **SEM 不受影响**：SEM 输入（`data/processed/sem/`）为独立的样本级标准化表，
  与本包 `03_run_sem.py` 逐行一致复现。

### 4.4 处置建议（需决策，见 `decisions_log.md` §D2）

- **方案 A（推荐用于投稿一致性）**：以归档节点表为准，即本包 `data/reference/`；
  WCI 用 `--nodes-dir data/reference` 复现手稿数字；方法学中注明"门水平未做 0.1% 过滤、
  含 Others_b/Unclassified 合并、关键种含 Neg"。
- **方案 B（推荐用于方法严谨性）**：以本包文档化方法（0.1% 过滤、无 Neg）重跑网络与 WCI，
  方法学与代码完全自洽，但主导因子与手稿现值不同，需同步修订正文。
- **方案 C**：明确区分"结构性描述量（WCI）"与"功能性结论"，接受两者主导因子不同。

---

## 5. 已知限制

1. 归档节点表来源不可再生（见 §4.2），本包以 `data/reference/` 形式保留。
2. `data/processed/sem/SEM_Cascade_Input.csv` 中的 `Env_*_PC1`、`Complexity/Robustness/Overall_Stability`
   为上游 PCA / 拓扑 / 鲁棒性分析的产出；本包将其作为**标准化样本级输入**，逐列来源见 `data_dictionary.md`。
3. 网络拓扑图（spring layout 可视化）在原项目中**无可复现脚本**，本包未包含绘图脚本；
   `Nodes_*.csv` 亦无坐标列。