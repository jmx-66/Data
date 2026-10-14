#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_run_sem.py — 级联三段分段 SEM（piecewise SEM, 标准化 OLS）
================================================================
输入 (data/processed/sem/)：
  SEM_Cascade_Input.csv  —— 样本级 n=24 数据，含
      pH, Complexity, Robustness, Overall_Stability,
      Env_Carbon_PC1, Env_Enzyme_PC1, Env_Community_PC1, Env_Function_PC1
  SEM_n24_Input.csv      —— 含 WaterLevel_binary, Microtopo_binary

输出 (results/tables/sem/)：
  SEM_Cascade_Results.csv  —— 全部模型的路径系数（β, p, R², n）
  SEM_Cascade_R2_Summary.csv

模型结构（生态级联）：
  环境(微地形/水位/pH)  →  中间四系统因子(碳/酶/群落/功能 PC1)  →  终点三指标(复杂性/抵抗力/恢复力)

方法：
  分段 OLS（piecewise SEM）：对每一段做标准化（z-score）多元线性回归，
  报告标准化系数 β、双尾 t 检验 p 值、R²、n；环境段同时给出 VIF。
  分组：整体 n=24；微地形=高丘/洼地 (各 n=12)；水位=高/低 (各 n=12)。

用法：
  python 03_run_sem.py
"""
import os
import argparse
import numpy as np
import pandas as pd
from scipy import stats as ss

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SEM_DIR = os.path.join(ROOT, 'data', 'processed', 'sem')
OUT = os.path.join(ROOT, 'results', 'tables', 'sem')

MID = {'Env_Carbon_PC1': '碳', 'Env_Enzyme_PC1': '酶',
       'Env_Community_PC1': '微生物群落', 'Env_Function_PC1': '功能'}
ENDS = {'Complexity': '复杂性', 'Robustness': '抵抗力',
        'Overall_Stability': '恢复力'}


def z(s):
    return (s - s.mean()) / s.std(ddof=0)


def vif(X):
    """预测变量矩阵的方差膨胀因子（无截距，逐列对其余列回归）"""
    out = []
    for i in range(X.shape[1]):
        y = X[:, i]
        others = np.delete(X, i, axis=1)
        Xc = np.column_stack([np.ones(len(y)), others])
        beta, *_ = np.linalg.lstsq(Xc, y, rcond=None)
        resid = y - Xc @ beta
        ss_res = float(np.sum(resid ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        out.append(float(1.0 / (1.0 - r2)) if r2 < 1 else np.inf)
    return np.array(out)


def ols(d, y, xs):
    """标准化 OLS：返回 (R², n, {x: (β, p)})"""
    zy = z(d[y])
    zxs = [z(d[x]) for x in xs]
    X = np.column_stack(zxs)
    Xc = np.column_stack([np.ones(len(d)), X])
    beta, *_ = np.linalg.lstsq(Xc, zy, rcond=None)
    resid = zy - Xc @ beta
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((zy - zy.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    n = len(d)
    dof = max(n - X.shape[1] - 1, 1)
    stderr = np.sqrt(ss_res / dof * np.linalg.inv(Xc.T @ Xc).diagonal())
    tvals = beta / stderr
    pvals = 2 * (1 - ss.t.cdf(np.abs(tvals), dof))
    return r2, n, {x: (float(b), float(pv)) for x, b, pv in zip(xs, beta[1:], pvals[1:])}


def run_model(d, env_xs):
    sig = lambda p: '**' if p < 0.01 else ('*' if p < 0.05 else '')
    zenv = np.column_stack([z(d[e]) for e in env_xs])
    vif_env = vif(zenv)
    rows = []
    for mid, mlabel in MID.items():
        r2, nn, coefs = ols(d, mid, env_xs)
        for x, (b, p) in coefs.items():
            rows.append({'Stage': '环境→中间', 'Target': mlabel, 'Source': x,
                         'Beta_std': round(b, 3), 'P': round(p, 4),
                         'sig': sig(p), 'R2': round(r2, 3), 'n': nn})
    for end, elabel in ENDS.items():
        r2, nn, coefs = ols(d, end, list(MID.keys()))
        for x, (b, p) in coefs.items():
            rows.append({'Stage': '中间→终点', 'Target': elabel, 'Source': MID[x],
                         'Beta_std': round(b, 3), 'P': round(p, 4),
                         'sig': sig(p), 'R2': round(r2, 3), 'n': nn})
    return pd.DataFrame(rows), {'n': len(d),
                                'vif_env': dict(zip(env_xs, np.round(vif_env, 3)))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sem-dir', default=SEM_DIR)
    ap.add_argument('--out', default=OUT)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    sc = pd.read_csv(os.path.join(args.sem_dir, 'SEM_Cascade_Input.csv'),
                     encoding='utf-8-sig')
    n24 = pd.read_csv(os.path.join(args.sem_dir, 'SEM_n24_Input.csv'),
                      encoding='utf-8-sig')
    for c in ['Sample']:
        sc[c] = sc[c].astype(str).str.strip()
        n24[c] = n24[c].astype(str).str.strip()

    df = sc[['Sample', 'Group', 'pH', 'Complexity', 'Robustness', 'Overall_Stability',
             'Env_Carbon_PC1', 'Env_Enzyme_PC1', 'Env_Community_PC1', 'Env_Function_PC1']].merge(
        n24[['Sample', 'WaterLevel_binary', 'Microtopo_binary', 'WaterLevel', 'Microtopo']],
        on='Sample', how='left')
    df = df.dropna(subset=['pH', 'WaterLevel_binary', 'Microtopo_binary']).reset_index(drop=True)

    print("=" * 66)
    print("03 级联三段 SEM（环境 → 中间四因子 → 终点三指标）")
    print(f"   样本级 n={len(df)} | 分段标准化 OLS")
    print("=" * 66)

    tri = {}
    def collect(d, env_xs, tag):
        rdf, meta = run_model(d, env_xs)
        tri[tag] = (rdf, meta)

    collect(df, ['pH', 'Microtopo_binary', 'WaterLevel_binary'], '整体n24')
    collect(df[df['Microtopo_binary'] == 1].reset_index(drop=True),
            ['pH', 'WaterLevel_binary'], '微地形=高丘(n=12)')
    collect(df[df['Microtopo_binary'] == 0].reset_index(drop=True),
            ['pH', 'WaterLevel_binary'], '微地形=洼地(n=12)')
    collect(df[df['WaterLevel_binary'] == 1].reset_index(drop=True),
            ['pH', 'Microtopo_binary'], '水位=高(n=12)')
    collect(df[df['WaterLevel_binary'] == 0].reset_index(drop=True),
            ['pH', 'Microtopo_binary'], '水位=低(n=12)')

    all_df = pd.concat([v[0].assign(Model=k) for k, v in tri.items()], ignore_index=True)
    all_df.to_csv(os.path.join(args.out, 'SEM_Cascade_Results.csv'),
                  index=False, encoding='utf-8-sig')

    summ = []
    for k, (rdf, meta) in tri.items():
        mids = rdf[rdf['Stage'] == '环境→中间'].groupby('Target')['R2'].first().round(3).to_dict()
        ends = rdf[rdf['Stage'] == '中间→终点'].groupby('Target')['R2'].first().round(3).to_dict()
        summ.append({'Model': k, 'n': meta['n'],
                     **{f'M:{a}': b for a, b in mids.items()},
                     **{f'E:{a}': b for a, b in ends.items()}})
    pd.DataFrame(summ).to_csv(os.path.join(args.out, 'SEM_Cascade_R2_Summary.csv'),
                              index=False, encoding='utf-8-sig')

    for k, (rdf, meta) in tri.items():
        print(f"\n########## 模型: {k}  (n={meta['n']}) ##########")
        print('环境变量 VIF:', meta['vif_env'])
        print(rdf.to_string(index=False))
    print(f"\n[输出] {args.out}")


if __name__ == '__main__':
    main()