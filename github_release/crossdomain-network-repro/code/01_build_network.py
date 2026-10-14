#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_build_network.py — 跨界跨类型共现网络构建（可复现实现）
================================================================
输入 (data/raw/)：
  16s.csv / its.csv / 18s.csv          —— 细菌/真菌/原生生物 ASV 丰度表（含 Taxonomy 列）
  env_network.csv                       —— 环境指标表（碳库 / 酶 / 理化，样本为行）
  picrust2_sample_level1_relative.tsv   —— 细菌 KEGG L1 相对丰度（PICRUSt2）
  funguild_sample_total_absolute.tsv    —— 真菌 FUNGuild 功能预测（绝对丰度）

输出 (results/tables/network/)：
  Nodes_{group}.csv     —— 节点表（拓扑指标 + Zi-Pi 角色 + 关键种标签）
  Edges_{group}.csv     —— 边表（Spearman r、p、正负、跨类别标签）
  Network_Topology_Summary.csv —— 每组网络拓扑汇总

方法要点（与手稿"材料与方法"一致）：
  1) 门水平聚合：ASV 按门求和 -> 相对丰度 -> 剔除平均相对丰度 < MIN_ABUND 的罕见门
  2) 节点四系统：Community(微生物群落) / Function(微生物功能) / Carbon(碳库) / Enzyme(酶化学计量)
  3) 关联度量：Spearman 秩相关；建边阈值 p < 0.05（n=6 时等价 |r| >= 0.8286）
  4) 模块划分：greedy modularity（Newman–Girvan 模块度）
  5) 关键种判定：介数中心性 Top 10%  或  Zi-Pi 模块角色(Connector/Module Hub/Network Hub)
     —— 明确不使用 Neg（负连接）判据

可配置参数见下方 CONFIG。运行：
    python 01_build_network.py [--min-abund 0.001]
"""
import os
import re
import sys
import argparse
import numpy as np
import pandas as pd
from scipy import stats
import networkx as nx
import warnings

warnings.filterwarnings('ignore')

# ------------------------------------------------------------------
# 路径（相对本脚本定位仓库根目录）
# ------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, 'data', 'raw')
OUT = os.path.join(ROOT, 'results', 'tables', 'network')
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------
GROUPS = ['LL', 'LM', 'HL', 'HM']
SYSTEM_TYPES = ['Community', 'Function', 'Carbon', 'Enzyme']
SPEARMAN_P = 0.05          # 建边显著性阈值
N_FUNGAL_GUILD = 15        # FUNGuild 取前 15 个 Guild，其余合并为 Others

# 碳库节点 -> env 列名映射（去冗余后 24 个碳节点）
CARBON_NODE_MAP = {
    'TC': 'TC', 'MBC': 'MBC', 'DOC': 'DOC',
    'DON': 'DON', 'MBN': 'MBN', 'CN': 'C/N',
    'POM_abs': 'POM', 'MAOM_abs': 'MAOM',
    'POM_C': 'POM-C', 'POM_N': 'POM-N',
    'MAOM_C': 'MAOM-C', 'MAOM_N': 'MAOM-N',
    'POMN_MAOMN': 'POM-N/MAOM-N', 'HI': 'HI',
    'FTIR_1029': 'IR_1029', 'FTIR_1630': 'IR_1630',
    'FTIR_2900': 'IR_2900', 'FTIR_3400': 'IR_3400',
    'IR_Al_Aro': 'IR_Al_Aro', 'IR_Aro_Poly': 'IR_Aro_Poly',
    'IR_Al_Poly': 'IR_Al_Poly', 'IR_Aro_OH': 'IR_Aro_OH',
    'CUE': 'CUE', 'MBC_TC': 'MBC/TC',
}
# 酶节点 -> env 列名映射（10 个酶化学计量节点）
ENZYME_NODE_MAP = {
    'BG': 'C_acq', 'CBH': 'BG_CBH_Ratio', 'LAP': 'N_acq', 'AP': 'P_acq',
    'POX': 'Oxidative', 'PER': 'Hydrolytic',
    'Vec_Angle': 'Vector_Angle', 'Vec_Length': 'Vector_Length',
    'Enz_CP': 'Enz_CP_log', 'Enz_NP': 'Enz_NP_log',
}


# ==================================================================
# Part 1: 读取与预处理
# ==================================================================
def read_asv(filepath, name):
    """读取 ASV 丰度表 -> (丰度矩阵[ASV x 样本], Taxonomy 数组, 样本列名)"""
    df = pd.read_csv(filepath, encoding='utf-8-sig')
    id_col = [c for c in df.columns if 'ASV_ID' in c.upper() or 'OTU_ID' in c.upper()]
    id_col = id_col[0] if id_col else df.columns[0]
    sample_cols = [c for c in df.columns if str(c).startswith(('LL', 'LM', 'HL', 'HM'))]
    tax_col = 'Taxonomy' if 'Taxonomy' in df.columns else df.columns[-1]
    abundance = df[sample_cols].copy()
    abundance.index = df[id_col].values
    abundance = abundance.apply(pd.to_numeric, errors='coerce').fillna(0)
    taxonomy = df[tax_col].values
    print(f"  {name}: {len(abundance)} ASVs, {len(sample_cols)} 样本")
    return abundance, taxonomy, sample_cols


def parse_taxonomy(tax_series, rank=1):
    """解析分类学字符串，rank=1 取门水平；无法解析记为 Unclassified"""
    out = []
    for t in tax_series:
        parts = re.split(r'[;,]', str(t))
        if len(parts) > rank:
            p = parts[rank].strip()
            p = re.sub(r'^[a-z]__', '', p)
            p = re.sub(r'^[A-Z]_\d+', '', p)
        else:
            p = 'Unclassified'
        if not p or p.lower() in ['unclassified', 'unknown', '-', '', 'none']:
            p = 'Unclassified'
        out.append(p)
    return np.array(out)


def aggregate_to_phylum(abundance, tax_series, kingdom_name, min_abund=0.001):
    """门水平聚合 + 相对丰度 + 罕见门过滤"""
    phyla = parse_taxonomy(tax_series)
    df = abundance.copy()
    df['Phylum'] = phyla
    df = df[df['Phylum'] != 'Unclassified']
    if len(df) == 0:
        print(f"  警告: {kingdom_name} 全为 Unclassified")
        return pd.DataFrame()
    agg = df.groupby('Phylum').sum()
    rel = agg.div(agg.sum(axis=0), axis=1)
    mean_abund = rel.mean(axis=1)
    n_before = len(rel)
    rel = rel[mean_abund >= min_abund]
    print(f"  {kingdom_name}: {n_before} 个门 -> 过滤({min_abund*100:.2f}%)后 {len(rel)} 个")
    return rel


def read_function(picrust_path, funguild_path):
    """读取细菌 KEGG L1 与真菌 FUNGuild 功能数据 -> (相对丰度矩阵, 相对丰度矩阵)"""
    picrust_rel = None
    if os.path.exists(picrust_path):
        p = pd.read_csv(picrust_path, sep='\t', encoding='utf-8-sig')
        cols = [c for c in p.columns if str(c).startswith(('LL', 'LM', 'HL', 'HM'))]
        idx = 'level1' if 'level1' in p.columns else p.columns[0]
        picrust_rel = p[cols].apply(pd.to_numeric, errors='coerce').fillna(0)
        picrust_rel.index = p[idx].values
        picrust_rel = picrust_rel.div(picrust_rel.sum(axis=0), axis=1)
        print(f"  细菌 KEGG L1: {len(picrust_rel)} 条")
    fun_guild_rel = None
    if os.path.exists(funguild_path):
        fg = pd.read_csv(funguild_path, sep='\t', encoding='utf-8-sig')
        scols = [c for c in fg.columns if str(c).startswith(('LL', 'LM', 'HL', 'HM'))]
        ab = fg[scols].apply(pd.to_numeric, errors='coerce').fillna(0)
        ab['Guild'] = fg['Guild'].values
        ab = ab[~ab['Guild'].isin(['-', '', 'Unassigned', None, 'nan'])]
        g = ab.groupby('Guild').sum()
        fun_guild_rel = g.div(g.sum(axis=0), axis=1)
        print(f"  真菌 FUNGuild: {len(fun_guild_rel)} 条 Guild")
    return picrust_rel, fun_guild_rel


def get_group_samples(env_index, group):
    """返回该组的 6 个样本名"""
    cand = [f'{group}_{i}' for i in range(1, 7)]
    valid = [s for s in cand if s in env_index]
    if len(valid) == 6:
        return valid
    cand2 = [f'{group}{i}' for i in range(1, 7)]
    valid2 = [s for s in cand2 if s in env_index]
    if len(valid2) == 6:
        return valid2
    matched = sorted([s for s in env_index if str(s).startswith(group)])
    return matched[:6] if len(matched) >= 6 else cand


def get_phylum_abundance(phylum_rel, group):
    if phylum_rel is None or len(phylum_rel) == 0:
        return None
    cols = [c for c in phylum_rel.columns if str(c).startswith(group)]
    cols = sorted(cols)[:6]
    if len(cols) < 6:
        print(f"    警告: {group} 门丰度样本不足 ({len(cols)})")
        return None
    return phylum_rel[cols].T.values


def get_function_abundance(func_rel, group):
    if func_rel is None:
        return None
    cols = sorted([c for c in func_rel.columns if str(c).startswith(group)])[:6]
    return func_rel[cols].T.values if len(cols) >= 6 else None


def extract_env_nodes(env_df, group, node_map):
    samples = get_group_samples(env_df.index, group)
    result = {}
    for node_name, col in node_map.items():
        if col in env_df.columns:
            vals = [float(env_df.loc[s, col]) for s in samples if s in env_df.index]
            if len(vals) == 6:
                result[node_name] = np.array(vals)
    return pd.DataFrame(result, index=samples)


# ==================================================================
# Part 2: 单组网络构建
# ==================================================================
def build_network_for_group(group, bac, fun, pro, env_df, picrust, funguild):
    bac_v = get_phylum_abundance(bac, group)
    fun_v = get_phylum_abundance(fun, group)
    pro_v = get_phylum_abundance(pro, group)
    env_carbon = extract_env_nodes(env_df, group, CARBON_NODE_MAP)
    env_enzyme = extract_env_nodes(env_df, group, ENZYME_NODE_MAP)
    bac_f = get_function_abundance(picrust, group)
    fun_f = get_function_abundance(funguild, group)
    n_samples = 6

    # 节点清单: (id, label, kingdom, system, category)
    node_info = []
    if bac_v is not None:
        for phy in bac.index:
            node_info.append((phy, phy, 'Bacteria', 'Community', 'Bac_Phylum'))
    if fun_v is not None:
        for phy in fun.index:
            node_info.append((phy, phy, 'Fungi', 'Community', 'Fun_Phylum'))
    if pro_v is not None:
        for phy in pro.index:
            node_info.append((phy, phy, 'Protist', 'Community', 'Pro_Phylum'))
    if picrust is not None:
        for fn in picrust.index:
            safe = fn.replace(' ', '_').replace('/', '_')
            node_info.append((f'B_F_{safe}', f'B_F_{safe}', 'Bacteria', 'Function', 'Bac_Function'))
    if funguild is not None:
        for gn in funguild.index[:N_FUNGAL_GUILD]:
            safe = gn.replace(' ', '_').replace('/', '_')
            node_info.append((f'F_F_{safe}', f'F_F_{safe}', 'Fungi', 'Function', 'Fun_Function'))
        if len(funguild) > N_FUNGAL_GUILD:
            node_info.append(('F_F_Others_ff', 'F_F_Others_ff', 'Fungi', 'Function', 'Fun_Function'))
    for c in env_carbon.columns:
        node_info.append((c, c, 'Carbon', 'Carbon', 'Carbon_Pool'))
    for c in env_enzyme.columns:
        node_info.append((c, c, 'Enzyme', 'Enzyme', 'Enzyme_Stoi'))

    n_nodes = len(node_info)
    if n_nodes == 0:
        return pd.DataFrame(), pd.DataFrame(), nx.Graph()

    # 丰度矩阵
    M = np.zeros((n_samples, n_nodes))
    valid = [True] * n_nodes
    for j, (nid, label, kingdom, system, cat) in enumerate(node_info):
        if system == 'Community' and kingdom == 'Bacteria' and bac_v is not None and label in list(bac.index):
            M[:, j] = bac_v[:, list(bac.index).index(label)]
        elif system == 'Community' and kingdom == 'Fungi' and fun_v is not None and label in list(fun.index):
            M[:, j] = fun_v[:, list(fun.index).index(label)]
        elif system == 'Community' and kingdom == 'Protist' and pro_v is not None and label in list(pro.index):
            M[:, j] = pro_v[:, list(pro.index).index(label)]
        elif cat == 'Bac_Function' and bac_f is not None:
            fn = label.replace('B_F_', '')
            hit = [k for k, o in enumerate(picrust.index)
                   if o.replace(' ', '_').replace('/', '_') == fn]
            if hit:
                M[:, j] = bac_f[:, hit[0]]
            else:
                valid[j] = False
        elif cat == 'Fun_Function' and fun_f is not None:
            if label == 'F_F_Others_ff':
                cols = sorted([c for c in funguild.columns if str(c).startswith(group)])[:6]
                M[:, j] = funguild.iloc[N_FUNGAL_GUILD:].sum(axis=0)[cols].values
            else:
                gn = label.replace('F_F_', '')
                hit = [o for o in funguild.index
                       if o.replace(' ', '_').replace('/', '_') == gn]
                if hit:
                    cols = sorted([c for c in funguild.columns if str(c).startswith(group)])[:6]
                    M[:, j] = funguild.loc[hit[0], cols].values
                else:
                    valid[j] = False
        elif system == 'Carbon' and label in env_carbon.columns:
            M[:, j] = env_carbon[label].values
        elif system == 'Enzyme' and label in env_enzyme.columns:
            M[:, j] = env_enzyme[label].values
        else:
            valid[j] = False

    keep = [j for j in range(n_nodes) if valid[j] and np.std(M[:, j]) > 0]
    if len(keep) < n_nodes:
        node_info = [node_info[j] for j in keep]
        M = M[:, keep]
        n_nodes = len(node_info)
    if n_nodes < 3:
        return pd.DataFrame(), pd.DataFrame(), nx.Graph()

    # Spearman 相关
    C = np.zeros((n_nodes, n_nodes))
    P = np.ones((n_nodes, n_nodes))
    for i in range(n_nodes):
        for j in range(i + 1, n_nodes):
            if np.std(M[:, i]) > 1e-10 and np.std(M[:, j]) > 1e-10:
                r, p = stats.spearmanr(M[:, i], M[:, j])
                if not np.isnan(r):
                    C[i, j] = C[j, i] = r
                    P[i, j] = P[j, i] = p

    # 建边
    edges = []
    for i in range(n_nodes):
        for j in range(i + 1, n_nodes):
            if P[i, j] < SPEARMAN_P:
                edges.append({
                    'Source': node_info[i][0], 'Target': node_info[j][0],
                    'Correlation': round(C[i, j], 4), 'P_value': round(P[i, j], 6),
                    'Type': 'Positive' if C[i, j] > 0 else 'Negative',
                    'Weight': abs(round(C[i, j], 4)),
                    'Source_Kingdom': node_info[i][2], 'Target_Kingdom': node_info[j][2],
                    'Source_Category': node_info[i][4], 'Target_Category': node_info[j][4],
                    'Link_Category': f"{node_info[i][4][0]}-{node_info[j][4][0]}",
                })
    edge_df = pd.DataFrame(edges)

    # 图
    G = nx.Graph()
    G.add_nodes_from([n[0] for n in node_info])
    if len(edge_df) > 0:
        G.add_edges_from(edge_df[['Source', 'Target']].values)

    deg = dict(G.degree())
    btw = nx.betweenness_centrality(G) if G.number_of_nodes() else {}
    clo = nx.closeness_centrality(G) if G.number_of_nodes() else {}
    neg_deg = {n: 0 for n in G.nodes()}
    pos_deg = {n: 0 for n in G.nodes()}
    if len(edge_df) > 0:
        for _, row in edge_df.iterrows():
            if row['Type'] == 'Negative':
                neg_deg[row['Source']] += 1; neg_deg[row['Target']] += 1
            else:
                pos_deg[row['Source']] += 1; pos_deg[row['Target']] += 1

    try:
        from networkx.algorithms.community import greedy_modularity_communities
        module = {}
        for i, comm in enumerate(greedy_modularity_communities(G)):
            for node in comm:
                module[node] = i + 1
    except Exception:
        module = {n: 1 for n in G.nodes()}

    rows = []
    for nid, label, kingdom, system, cat in node_info:
        r = {'Id': nid, 'Label': label,
             'Degree': deg.get(nid, 0),
             'Betweenness': round(btw.get(nid, 0), 6),
             'Closeness': round(clo.get(nid, 0), 6),
             'Kingdom': kingdom, 'System': system, 'Category': cat,
             'Module': module.get(nid, 1),
             'Neg_Degree': neg_deg.get(nid, 0), 'Pos_Degree': pos_deg.get(nid, 0)}
        r['Neg_Ratio'] = round(r['Neg_Degree'] / r['Degree'], 4) if r['Degree'] > 0 else 0
        r['Neg_Importance'] = round(r['Betweenness'] * r['Neg_Ratio'], 6)
        rows.append(r)
    node_df = pd.DataFrame(rows)
    node_df = calculate_zipi(node_df, G)
    node_df = identify_keystone(node_df)
    print(f"  {group}: {len(node_df)} 节点, {len(edge_df)} 条边")
    return node_df, edge_df, G


def calculate_zipi(node_df, G):
    """Zi(模块内连接度) / Pi(模块间连接度) -> ZiPi_Role"""
    node_df['Zi'] = 0.0
    node_df['Pi'] = 0.0
    node_df['ZiPi_Role'] = 'Peripheral'
    if G.number_of_nodes() == 0:
        return node_df
    mod_of = dict(zip(node_df['Id'], node_df['Module']))
    for idx, row in node_df.iterrows():
        nid = row['Id']
        if nid not in G:
            continue
        nbrs = list(G.neighbors(nid))
        if not nbrs:
            continue
        module = row['Module']
        k_in = sum(1 for nb in nbrs if mod_of.get(nb) == module)
        k_in_all = [G.degree(mn) for mn in node_df[node_df['Module'] == module]['Id'] if mn in G]
        mean_k = np.mean(k_in_all) if k_in_all else 0
        std_k = np.std(k_in_all) if k_in_all else 0
        node_df.at[idx, 'Zi'] = round((k_in - mean_k) / std_k if std_k > 1e-10 else 0, 4)
        k_is = sum(1 for nb in nbrs if mod_of.get(nb) != module)
        node_df.at[idx, 'Pi'] = round(k_is / len(nbrs), 4)
        zi, pi = node_df.at[idx, 'Zi'], node_df.at[idx, 'Pi']
        if zi > 2.5 and pi > 0.62:
            node_df.at[idx, 'ZiPi_Role'] = 'Network Hub'
        elif zi > 2.5:
            node_df.at[idx, 'ZiPi_Role'] = 'Module Hub'
        elif pi > 0.62:
            node_df.at[idx, 'ZiPi_Role'] = 'Connector'
    return node_df


def identify_keystone(node_df):
    """关键种 = 介数 Top10%  ∪  Zi-Pi 模块角色（不使用 Neg）"""
    if len(node_df) == 0:
        return node_df
    thr = np.percentile(node_df['Betweenness'].values, 90)
    node_df['Btw_Rank'] = node_df['Betweenness'].rank(ascending=False, method='min').astype(int)
    node_df['Keystone_Btw'] = node_df['Betweenness'].apply(
        lambda x: 'Top10pct_Btw' if x >= thr else 'Others')
    node_df['Keystone_ZiPi'] = node_df['ZiPi_Role'].apply(
        lambda x: 'Keystone' if x in ['Connector', 'Module Hub', 'Network Hub'] else 'Others')
    node_df['Keystone'] = node_df.apply(
        lambda r: 'Keystone' if (r['Keystone_Btw'] == 'Top10pct_Btw' or
                                 r['Keystone_ZiPi'] == 'Keystone') else 'Others', axis=1)
    node_df['Keystone_Type'] = node_df.apply(
        lambda r: ('Both' if (r['Keystone_Btw'] == 'Top10pct_Btw' and r['Keystone_ZiPi'] == 'Keystone')
                   else 'Betweenness' if r['Keystone_Btw'] == 'Top10pct_Btw'
                   else 'ZiPi' if r['Keystone_ZiPi'] == 'Keystone' else 'Others'), axis=1)
    return node_df


# ==================================================================
# Part 3: 主流程
# ==================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-abund', type=float, default=0.001,
                    help='门水平相对丰度过滤阈值（默认 0.001 = 0.1%%）')
    args = ap.parse_args()
    min_abund = args.min_abund

    print("=" * 66)
    print(f"01 跨界网络构建  (门水平过滤阈值 = {min_abund*100:.3f}%)")
    print("=" * 66)

    asv_16s, tax_16s, _ = read_asv(os.path.join(RAW, '16s.csv'), '16S细菌')
    asv_its, tax_its, _ = read_asv(os.path.join(RAW, 'its.csv'), 'ITS真菌')
    asv_18s, tax_18s, _ = read_asv(os.path.join(RAW, '18s.csv'), '18S原生生物')
    bac = aggregate_to_phylum(asv_16s, tax_16s, 'Bacteria', min_abund)
    fun = aggregate_to_phylum(asv_its, tax_its, 'Fungi', min_abund)
    pro = aggregate_to_phylum(asv_18s, tax_18s, 'Protist', min_abund)

    env_df = pd.read_csv(os.path.join(RAW, 'env_network.csv'), encoding='utf-8-sig')
    for key in ('Sample', 'sample'):
        if key in env_df.columns:
            env_df = env_df.set_index(key)
            break
    print(f"  环境表: {len(env_df)} 样本 × {len(env_df.columns)} 列")

    picrust, funguild = read_function(
        os.path.join(RAW, 'picrust2_sample_level1_relative.tsv'),
        os.path.join(RAW, 'funguild_sample_total_absolute.tsv'))

    topo = []
    for g in GROUPS:
        node_df, edge_df, G = build_network_for_group(g, bac, fun, pro, env_df, picrust, funguild)
        node_df['Group'] = g
        if len(node_df) > 0:
            node_df.to_csv(os.path.join(OUT, f'Nodes_{g}.csv'), index=False, encoding='utf-8-sig')
        if len(edge_df) > 0:
            edge_df.to_csv(os.path.join(OUT, f'Edges_{g}.csv'), index=False, encoding='utf-8-sig')
        n_nodes, n_edges = len(node_df), len(edge_df)
        n_pos = int((edge_df['Type'] == 'Positive').sum()) if n_edges else 0
        n_neg = int((edge_df['Type'] == 'Negative').sum()) if n_edges else 0
        density = 2 * n_edges / (n_nodes * (n_nodes - 1)) if n_nodes > 1 else 0
        avg_deg = 2 * n_edges / n_nodes if n_nodes > 0 else 0
        try:
            modularity = nx.algorithms.community.modularity(
                G, [set(node_df[node_df['Module'] == m]['Id']) for m in node_df['Module'].unique()])
        except Exception:
            modularity = 0
        clustering = nx.average_clustering(G) if G.number_of_nodes() else 0
        topo.append({'Group': g, 'Nodes': n_nodes, 'Edges': n_edges,
                     'Positive_Edges': n_pos, 'Negative_Edges': n_neg,
                     'Density': round(density, 4), 'Avg_Degree': round(avg_deg, 2),
                     'Clustering_Coef': round(clustering, 4), 'Modularity': round(modularity, 4)})

    topo_df = pd.DataFrame(topo)
    topo_df.to_csv(os.path.join(OUT, 'Network_Topology_Summary.csv'), index=False, encoding='utf-8-sig')
    print("\n拓扑汇总:")
    print(topo_df.to_string(index=False))
    print(f"\n[输出目录] {OUT}")


if __name__ == '__main__':
    main()