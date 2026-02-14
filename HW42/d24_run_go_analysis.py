#!/usr/bin/env python3
# C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d24_run_go_analysis.py
"""
GO Enrichment & Depletion Analysis for all 3 Aspects (C, P, F)
使用 Hypergeometric Test 進行富集和耗盡分析

輸入檔案:
  - d01_live_genes_only.csv (背景基因集)
  - d23_example_input.txt (輸入基因列表)
  - d22_go_features_C/P/F.csv (Feature Tables)

輸出檔案 (6個):
  - d24_go_C_enrichment.csv, d24_go_C_depletion.csv
  - d24_go_P_enrichment.csv, d24_go_P_depletion.csv
  - d24_go_F_enrichment.csv, d24_go_F_depletion.csv

統計方法: Hypergeometric Test + Benjamini-Hochberg FDR 校正
"""

import pandas as pd
import numpy as np
import os
import sys
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


def load_background_genes(background_file):
    """
    載入背景基因集 (N)
    
    參數:
        background_file: d01_live_genes_only.csv 檔案路徑
    
    返回:
        set: 所有背景基因的 ID 集合
    """
    try:
        df = pd.read_csv(background_file)
        genes = set(df['WBGene_ID'].str.strip())
        genes = {g for g in genes if g}  # 移除空字串
        print(f"✓ 背景基因集載入完成: {len(genes)} 個基因")
        return genes
    except Exception as e:
        print(f"✗ 載入背景基因集失敗: {str(e)}")
        sys.exit(1)


def load_input_genes(input_file):
    """
    載入輸入基因列表 (n)
    
    參數:
        input_file: d23_example_input.txt 檔案路徑
    
    返回:
        set: 輸入基因的 ID 集合
    """
    try:
        genes = set()
        with open(input_file, 'r', encoding='utf-8') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    genes.add(gene)
        print(f"✓ 輸入基因集載入完成: {len(genes)} 個基因")
        return genes
    except Exception as e:
        print(f"✗ 載入輸入基因集失敗: {str(e)}")
        sys.exit(1)


def calculate_statistics(N, n, M, k):
    """
    計算統計數值
    
    參數:
        N: 總基因數 (背景)
        n: 輸入基因數
        M: 特定 GO term 的基因數
        k: 輸入基因中屬於該 GO term 的數量 (Overlap)
    
    返回:
        dict: 包含各項統計數值
    """
    # 預期 k 值
    expected_k = (M * n) / N if N > 0 else 0
    
    # Fold Change
    observed_ratio = k / n if n > 0 else 0
    expected_ratio = expected_k / n if n > 0 else 0
    fold_change = observed_ratio / expected_ratio if expected_ratio > 0 else 0
    
    # Log 轉換 (处理 zero division)
    log2_fc = np.log2(fold_change) if fold_change > 0 else 0
    
    # Hypergeometric Test
    # Enrichment: 使用 sf (survival function) = 1 - cdf
    # Depletion: 使用 cdf
    pval_enrichment = hypergeom.sf(k - 1, N, M, n)  # k 及以上的概率
    pval_depletion = hypergeom.cdf(k, N, M, n)      # k 及以下的概率
    
    # -log10 轉換 (处理 zero)
    log10_pval_enrichment = -np.log10(pval_enrichment) if pval_enrichment > 0 else 0
    log10_pval_depletion = -np.log10(pval_depletion) if pval_depletion > 0 else 0
    
    return {
        'observed_k': k,
        'expected_k': expected_k,
        'observed_ratio': observed_ratio,
        'expected_ratio': expected_ratio,
        'fold_change': fold_change,
        'log2_fc': log2_fc,
        'pval_enrichment': pval_enrichment,
        'pval_depletion': pval_depletion,
        'log10_pval_enrichment': log10_pval_enrichment,
        'log10_pval_depletion': log10_pval_depletion
    }


def format_ratio(k, n, ratio):
    """格式化 Observed/Expected Ratio"""
    if n == 0:
        return "0/0 (0%)"
    pct = (k / n) * 100
    return f"{k}/{n} ({pct:.1f}%)"


def run_analysis(aspect_code, feature_file, output_dir, N, n, input_genes, background_genes):
    """
    針對特定 Aspect 進行富集和耗盡分析
    
    參數:
        aspect_code: 'C', 'P', 或 'F'
        feature_file: d22_go_features_{code}.csv 檔案路徑
        output_dir: 輸出目錄
        N: 背景基因總數
        n: 輸入基因數
        input_genes: 輸入基因集合
        background_genes: 背景基因集合
    """
    
    print(f"\n{'='*70}")
    print(f"處理 Aspect {aspect_code}...")
    print(f"{'='*70}")
    
    # 檢查特徵檔案
    if not os.path.exists(feature_file):
        print(f"⊘ 警告: 找不到特徵檔案 {feature_file}，跳過該 Aspect")
        return
    
    # 載入特徵表
    try:
        df_features = pd.read_csv(feature_file)
        print(f"✓ 特徵表載入完成: {len(df_features)} 個 GO terms")
    except Exception as e:
        print(f"✗ 載入特徵表失敗: {str(e)}")
        return
    
    # 檢查必要欄位
    if not all(col in df_features.columns for col in ['Feature Name', 'Number', 'Gene']):
        print(f"✗ 特徵表缺少必要欄位")
        return
    
    # 準備結果列表
    enrichment_results = []
    depletion_results = []
    
    # 對每個 GO term 進行分析
    for idx, row in df_features.iterrows():
        feature_name = row['Feature Name']
        M = row['Number']  # 該 GO term 中的基因數
        
        # 提取該 GO term 的基因
        feature_genes = set(g.strip() for g in str(row['Gene']).split(','))
        feature_genes = {g for g in feature_genes if g}
        
        # 計算 overlap
        k = len(input_genes & feature_genes)
        
        # 計算統計
        stats = calculate_statistics(N, n, M, k)
        
        # 富集結果
        enrichment_results.append({
            'Feature Name': feature_name,
            'Observed Ratio': format_ratio(stats['observed_k'], n, stats['observed_ratio']),
            'Expected Ratio': format_ratio(int(stats['expected_k']), n, stats['expected_ratio']),
            'Fold Change': stats['fold_change'],
            'P-value': stats['pval_enrichment'],
            'log2(Fold Change)': stats['log2_fc'],
            '-log10(P-value)': stats['log10_pval_enrichment']
        })
        
        # 耗盡結果
        depletion_results.append({
            'Feature Name': feature_name,
            'Observed Ratio': format_ratio(stats['observed_k'], n, stats['observed_ratio']),
            'Expected Ratio': format_ratio(int(stats['expected_k']), n, stats['expected_ratio']),
            'Fold Change': stats['fold_change'],
            'P-value': stats['pval_depletion'],
            'log2(Fold Change)': stats['log2_fc'],
            '-log10(P-value)': stats['log10_pval_depletion']
        })
    
    print(f"✓ 分析完成: {len(enrichment_results)} 個 GO terms")
    
    # 轉換為 DataFrame
    df_enrichment = pd.DataFrame(enrichment_results)
    df_depletion = pd.DataFrame(depletion_results)
    
    # 應用多重測試校正 (Benjamini-Hochberg FDR 和 Bonferroni)
    for df in [df_enrichment, df_depletion]:
        # FDR 校正
        reject, fdr_pvals, _, _ = multipletests(
            df['P-value'], 
            method='fdr_bh',
            is_sorted=False,
            returnsorted=False
        )
        df['FDR P-value'] = fdr_pvals
        
        # Bonferroni 校正
        df['Bonferroni P-value'] = np.minimum(df['P-value'] * len(df), 1.0)
        
        # -log10 FDR
        df['-log10(FDR P-value)'] = df['FDR P-value'].apply(
            lambda x: -np.log10(x) if x > 0 else 0
        )
    
    # 重新排列欄位順序
    col_order = [
        'Feature Name',
        'Observed Ratio',
        'Expected Ratio',
        'Fold Change',
        'P-value',
        'FDR P-value',
        'Bonferroni P-value',
        'log2(Fold Change)',
        '-log10(P-value)',
        '-log10(FDR P-value)'
    ]
    
    df_enrichment = df_enrichment[col_order].sort_values('P-value')
    df_depletion = df_depletion[col_order].sort_values('P-value')
    
    # 統計顯著結果 (P < 0.05)
    sig_enrichment = (df_enrichment['P-value'] < 0.05).sum()
    sig_depletion = (df_depletion['P-value'] < 0.05).sum()
    sig_enrichment_fdr = (df_enrichment['FDR P-value'] < 0.05).sum()
    sig_depletion_fdr = (df_depletion['FDR P-value'] < 0.05).sum()
    
    # 輸出結果
    enrichment_file = os.path.join(output_dir, f'd24_go_{aspect_code}_enrichment.csv')
    depletion_file = os.path.join(output_dir, f'd24_go_{aspect_code}_depletion.csv')
    
    df_enrichment.to_csv(enrichment_file, index=False)
    df_depletion.to_csv(depletion_file, index=False)
    
    print(f"\n📊 Aspect {aspect_code} 分析結果摘要:")
    print(f"  ├─ 富集 (Enrichment):")
    print(f"  │  ├─ P < 0.05: {sig_enrichment}")
    print(f"  │  └─ FDR < 0.05: {sig_enrichment_fdr}")
    print(f"  ├─ 耗盡 (Depletion):")
    print(f"  │  ├─ P < 0.05: {sig_depletion}")
    print(f"  │  └─ FDR < 0.05: {sig_depletion_fdr}")
    print(f"  ├─ 輸出檔案:")
    print(f"  │  ├─ {os.path.basename(enrichment_file)}")
    print(f"  │  └─ {os.path.basename(depletion_file)}")
    
    return df_enrichment, df_depletion


def main():
    # 設定檔案路徑
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    
    background_file = os.path.join(script_dir, 'd01_live_genes_only.csv')
    input_file = os.path.join(script_dir, 'd23_example_input.txt')
    
    output_dir = script_dir
    
    print(f"開始 GO 富集和耗盡分析...")
    print(f"")
    
    # 載入背景基因和輸入基因
    background_genes = load_background_genes(background_file)
    input_genes = load_input_genes(input_file)
    
    N = len(background_genes)  # 總基因數
    n = len(input_genes)       # 輸入基因數
    
    print(f"\n參數:")
    print(f"  N (背景基因總數): {N}")
    print(f"  n (輸入基因數): {n}")
    print()
    
    # 對每個 Aspect 進行分析
    for aspect_code in ['C', 'P', 'F']:
        feature_file = os.path.join(script_dir, f'd22_go_features_{aspect_code}.csv')
        run_analysis(aspect_code, feature_file, output_dir, N, n, input_genes, background_genes)
    
    print(f"\n{'='*70}")
    print(f"✓ GO 分析完成！")
    print(f"{'='*70}")
    print(f"\n生成的檔案:")
    print(f"  ├─ d24_go_C_enrichment.csv / d24_go_C_depletion.csv")
    print(f"  ├─ d24_go_P_enrichment.csv / d24_go_P_depletion.csv")
    print(f"  └─ d24_go_F_enrichment.csv / d24_go_F_depletion.csv")


if __name__ == '__main__':
    main()
