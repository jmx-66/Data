# 归档（手稿）结果说明（`data/reference/`）

本目录保存**手稿所用的历史归档结果**，用途有二：

1. **对照验证**——证明本包脚本能逐位复现手稿数字（见 `docs/reproducibility.md` §3）；
2. **投稿一致性**——当需要复现手稿既有数字时，直接以本目录节点表为输入。

> ⚠️ **本目录数据不参与本包默认计算**，默认管道从 `data/raw/` 重新构建网络。

## 文件清单

| 文件 | 内容 |
|:--|:--|
| `Nodes_{group}.csv` | 归档节点表（**23 列**，含 `Neg_Rank`、`Keystone_Neg`） |
| `Edges_{group}.csv` | 归档边表 |
| `Network_Topology_Summary.csv` | 归档拓扑汇总 |
| `Keystone_Statistics_v6.csv` | 归档关键种统计（含 `Neg_Keystones` 列，即 Neg 口径证据） |
| `WCI_v3_3comp_Final_Results.csv` | 手稿 WCI 结果（16 格） |
| `WCI_v3_3comp_Final_Weights.csv` | 手稿 WCI 熵权 |
| `SEM_Cascade_Archived_Results.csv` | 手稿级联 SEM 路径系数（104 行） |

## 复现手稿数字

```bash
# WCI（手稿口径）
python code/02_calculate_wci.py --nodes-dir data/reference --out results/tables/wci_archived

# SEM（本包脚本即与归档逐行一致）
python code/03_run_sem.py
```

实测：WCI 与 `WCI_v3_3comp_Final_Results.csv` 最大绝对误差 **7.8×10⁻¹⁶**，主导因子完全一致；
SEM 与 `SEM_Cascade_Archived_Results.csv` **104 行逐行一致**（误差 0，显著性一致率 100%）。

## 与现管道（`code/01`）的口径差异

归档节点表**无法由任何现存脚本逐字节再生**，与现管道存在如下差异（详见 `docs/decisions_log.md` §D2、
`docs/reproducibility.md` §4）：

| 项 | 归档 | 现管道 |
|:--|:--|:--|
| 列数 | 23（含 Neg 相关列） | 22 |
| 门水平过滤 | 无 0.1% 过滤 | 0.1% 平均相对丰度过滤 |
| 合并节点 | 含 `Others_b`、`Unclassified` | 不含 |
| 关键种判据 | 介数 ∪ ZiPi ∪ **Neg** | 介数 ∪ ZiPi（无 Neg） |
| WCI 主导因子 | LL→酶、LM→碳、HL→碳、HM→功能 | LL→功能、LM→功能、HL→功能、HM→群落 |

归档门的合并规则与阈值在多阈值实测下**均无法精确重现**，故本目录以原样保留、不再重算。