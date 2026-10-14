#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02_calculate_wci.py — 加权贡献指数 WCI（三分量最终版）
================================================================
输入 (默认 results/tables/network/)：
  Nodes_{group}.csv   —— 01_build_network.py 的节点表
                         必需列：Group, System, Betweenness, Degree, Keystone

输出 (results/tables/wci/)：
  WCI_Results.csv   —— 16 格（4 组 × 4 系统）WCI 值、组内排名、主导因子
  WCI_Weights.csv   —— 三分量熵权（用原始数据计算）

方法（与手稿"材料与方法"一致，见 docs/methods.md）：
  分量 1  BWKI  Betweenness-Weighted Keystone Intensity = 该系统全部节点的介数中心性均值
  分量 2  DCI   Degree-Weighted Contribution Intensity  = 该系统全部节点的度均值
  分量 3  NKR   Node-type Normalized Keystone Rate      = 观测关键种率 / 随机期望关键种率
                （Monte Carlo 置换 N_MC=100000，从其他系统随机抽取等量节点）
  权重    熵权法，直接对三分量"原始值"计算（p_ij = x_ij / Σx_ij），先定权
  归一化  Max-scaling（x / x_max），仅用于统一量纲，保比例
  合成    WCI = Σ_j w_j · x_norm,j
  主导因子 组内 WCI 最高者（= 该节点类别在整合网络中结构性关键节点富集最强）

用法：
  python 02_calculate_wci.py
  python 02_calculate_wci.py --nodes-dir ../../data/reference   # 复现归档（手稿）口径
"""
import os
import argparse
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DEFAULT_NODES = os.path.join(ROOT, 'results', 'tables', 'network')
OUT = os.path.join(ROOT, 'results', 'tables', 'wci')

GROUPS = ['LL', 'LM', 'HL', 'HM']
SYSTEMS = ['Community', 'Function', 'Carbon', 'Enzyme']
SYS_LABELS = {'Community': '微生物群落', 'Function': '微生物功能',
              'Carbon': '碳库指标', 'Enzyme': '酶化学计量'}
COMPONENTS = ['BWKI', 'DCI', 'NKR']

N_MC = 100000
SEED = 42


def load_nodes(nodes_dir):
    nodes = {}
    for g in GROUPS:
        fp = os.path.join(nodes_dir, f'Nodes_{g}.csv')
        df = pd.read_csv(fp, encoding='utf-8-sig')
        if 'Group' not in df.columns:
            df['Group'] = g
        nodes[g] = df
        nk = int((df['Keystone'] == 'Keystone').sum())
        print(f"  {g}: {len(df)} 节点, {nk} 关键种 ({nk/len(df)*100:.1f}%)")
    return nodes


def component_bwki(nodes):
    rows = []
    for g in GROUPS:
        ndf = nodes[g]
        for s in SYSTEMS:
            sub = ndf[ndf['System'] == s]
            if len(sub) == 0:
                continue
            rows.append({'Group': g, 'System': s, 'N_Nodes': len(sub),
                         'N_Keystones': int((sub['Keystone'] == 'Keystone').sum()),
                         'BWKI': sub['Betweenness'].sum() / len(sub)})
    return pd.DataFrame(rows)


def component_dci(nodes):
    rows = []
    for g in GROUPS:
        ndf = nodes[g]
        for s in SYSTEMS:
            sub = ndf[ndf['System'] == s]
            if len(sub) == 0:
                continue
            rows.append({'Group': g, 'System': s, 'DCI': sub['Degree'].sum() / len(sub)})
    return pd.DataFrame(rows)


def component_nkr(nodes, n_mc=N_MC, seed=SEED):
    """NKR = 观测关键种率 / 随机期望率；置换从"其他系统"抽等量节点。

    注：为与归档（手稿）结果逐位一致，此处沿用 np.random.seed(seed) +
    np.random.choice 的全局 RNG 口径（等价于原脚本 pandas .sample 的随机序列）。
    """
    np.random.seed(seed)
    rows = []
    for g in GROUPS:
        ndf = nodes[g].copy()
        ndf['is_keystone'] = (ndf['Keystone'] == 'Keystone').astype(int)
        for s in SYSTEMS:
            sub = ndf[ndf['System'] == s]
            n_type = len(sub)
            if n_type == 0:
                continue
            obs_rate = sub['is_keystone'].sum() / n_type
            other = ndf[ndf['System'] != s]['is_keystone'].values
            mc = np.empty(n_mc)
            rep = n_type > len(other)
            for i in range(n_mc):
                idx = np.random.choice(len(other), size=n_type, replace=rep)
                mc[i] = other[idx].mean()
            mc_mean = mc.mean()
            nkr = obs_rate / mc_mean if mc_mean > 0 else 0.0
            p = (np.sum(mc >= obs_rate) + 1) / (n_mc + 1)
            if obs_rate > np.percentile(mc, 97.5):
                enrich = '显著富集'
            elif obs_rate < np.percentile(mc, 2.5):
                enrich = '显著缺失'
            else:
                enrich = '随机水平'
            rows.append({'Group': g, 'System': s, 'NKR': nkr,
                         'P_Value': p, 'Enrichment': enrich})
            print(f"  {g} {s}: NKR={nkr:.4f}, p={p:.6f} ({enrich})")
    return pd.DataFrame(rows)


def entropy_weights(mat_raw, components):
    """熵权：p_ij = x_ij/Σx_ij（原始值比例），先定权"""
    n = mat_raw.shape[0]
    p = mat_raw / mat_raw.sum(axis=0, keepdims=True)
    p = np.clip(p, 1e-300, None)
    e = -np.sum(p * np.log(p), axis=0) / np.log(n)
    d = 1 - e
    w = d / d.sum()
    return e, d, w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nodes-dir', default=DEFAULT_NODES,
                    help='节点表目录（默认 results/tables/network）')
    ap.add_argument('--n-mc', type=int, default=N_MC)
    ap.add_argument('--seed', type=int, default=SEED)
    ap.add_argument('--out', default=OUT)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print("=" * 62)
    print("02 WCI 计算（三分量 BWKI + DCI + NKR）")
    print(f"   节点表: {args.nodes_dir}")
    print("=" * 62)
    nodes = load_nodes(args.nodes_dir)

    print("\n[1] 分量 BWKI（全节点介数均值）")
    bwki = component_bwki(nodes)
    print("\n[2] 分量 DCI（全节点度均值）")
    dci = component_dci(nodes)
    print(f"\n[3] 分量 NKR（Monte Carlo = {args.n_mc:,}）")
    nkr = component_nkr(nodes, args.n_mc, args.seed)

    wci = bwki.merge(dci, on=['Group', 'System'], how='left') \
              .merge(nkr, on=['Group', 'System'], how='left')

    print("\n[4] 熵权（原始数据）+ Max-scaling 合成")
    mat_raw = wci[COMPONENTS].values.astype(float)
    e, d, w = entropy_weights(mat_raw, COMPONENTS)
    for i, c in enumerate(COMPONENTS):
        print(f"    {c}: 熵e={e[i]:.4f}  差异d={d[i]:.4f}  权重w={w[i]:.4f}")
    print(f"    权重和 = {w.sum():.4f}")

    mat_ms = mat_raw / mat_raw.max(axis=0, keepdims=True)
    for i, c in enumerate(COMPONENTS):
        wci[f'{c}_norm'] = mat_ms[:, i]
    wci['WCI'] = sum(w[i] * mat_ms[:, i] for i in range(len(COMPONENTS)))
    wci['Rank'] = wci.groupby('Group')['WCI'].rank(ascending=False).astype(int)

    print("\n[5] 各组主导因子")
    for g in GROUPS:
        t = wci[wci['Group'] == g].sort_values('WCI', ascending=False).iloc[0]
        print(f"    {g}: {SYS_LABELS[t['System']]} (WCI={t['WCI']:.4f})")

    wci.to_csv(os.path.join(args.out, 'WCI_Results.csv'),
               index=False, encoding='utf-8-sig')
    pd.DataFrame({'Component': COMPONENTS, 'Entropy': e,
                  'Difference': d, 'Weight': w}).to_csv(
        os.path.join(args.out, 'WCI_Weights.csv'), index=False, encoding='utf-8-sig')
    print(f"\n[输出] {args.out}")
    print(wci[['Group', 'System', 'BWKI', 'DCI', 'NKR', 'WCI', 'Rank']].to_string(index=False))


if __name__ == '__main__':
    main()