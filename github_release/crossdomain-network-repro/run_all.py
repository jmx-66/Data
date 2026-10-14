#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_all.py — 一键端到端复现
顺序执行：01 网络构建 → 02 WCI 计算 → 03 级联 SEM
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, 'code')

STEPS = [
    ('01_build_network.py', '跨界网络构建'),
    ('02_calculate_wci.py', 'WCI 三分量计算'),
    ('03_run_sem.py', '级联三段 SEM'),
]


def main():
    for script, desc in STEPS:
        print('\n' + '=' * 70)
        print(f'>>> {desc}  ({script})')
        print('=' * 70)
        rc = subprocess.call([sys.executable, os.path.join(CODE, script)])
        if rc != 0:
            print(f'[失败] {script} 退出码 {rc}')
            sys.exit(rc)
    print('\n全部完成。结果见 results/tables/')


if __name__ == '__main__':
    main()