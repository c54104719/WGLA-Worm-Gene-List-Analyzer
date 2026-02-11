"""
完整的生物信息學富集分析報告生成腳本
執行超幾何分佈檢定 (Hypergeometric Test) - Depletion Only
生成包含多重檢驗校正、比值和對數轉換的完整報告
"""

import pandas as pd
import numpy as np
import math
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# ========== 檔案路徑設定 ==========
population_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d01_live_genes_only.csv"
feature_table_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d12_fgg_feature_table.csv"
input_genes_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d13_example_input.txt"
output_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d14_fgg_depletion_full.csv"

print("=" * 80)
print("完整富集分析報告生成器")
print("超幾何分佈檢定 (Hypergeometric Test) - Depletion Only")
print("=" * 80)

# ========== 步驟1: 讀取背景母體 ==========
print("\n[步驟1] 讀取背景母體 (Live 基因)...")
population_df = pd.read_csv(population_file)
N = len(population_df)  # 總基因數
print(f"  ✓ 背景母體總數 N = {N}")

# ========== 步驟2: 讀取輸入基因清單 ==========
print(f"\n[步驟2] 讀取輸入基因清單...")
try:
    with open(input_genes_file, 'r', encoding='utf-8') as f:
        input_genes = set(line.strip() for line in f if line.strip())
except FileNotFoundError:
    print(f"  ✗ 找不到輸入檔案: {input_genes_file}")
    exit(1)

n = len(input_genes)
print(f"  ✓ 輸入基因數 n = {n}")

# ========== 步驟3: 讀取特徵表 ==========
print(f"\n[步驟3] 讀取特徵表...")
try:
    feature_df = pd.read_csv(feature_table_file)
except FileNotFoundError:
    print(f"  ✗ 找不到特徵表: {feature_table_file}")
    exit(1)

print(f"  ✓ 特徵數量: {len(feature_df)}")

# ========== 步驟4: 執行去耗盡分析並計算所有統計量 ==========
print(f"\n[步驟4] 執行超幾何分佈檢定並計算統計量 (Depletion)...")

results = []
pvalues_list = []

for idx, row in feature_df.iterrows():
    feature_name = row['Gene list']
    M = int(row['Number'])  # 特徵群成員數
    
    # 解析基因 ID (逗號分隔)
    gene_list_str = row['Gene']
    if pd.isna(gene_list_str) or gene_list_str == '':
        continue
    
    feature_genes = set(g.strip() for g in str(gene_list_str).split(','))
    
    # 計算交集
    overlap_genes = input_genes & feature_genes
    k = len(overlap_genes)
    
    # ===== 計算比值 =====
    # Expected Ratio = (n * M / N) / n = M / N
    expected_count = (n * M) / N
    expected_ratio = M / N if N > 0 else 0
    
    # Observed Ratio = k / n
    observed_ratio = k / n if n > 0 else 0
    
    # ===== Fold Change =====
    if expected_ratio > 0:
        fold_change = observed_ratio / expected_ratio
    else:
        fold_change = 0
    
    # ===== P-value (Depletion) =====
    # cdf(k, N, M, n) = P(X <= k)
    depletion_pvalue = hypergeom.cdf(k, N, M, n)
    pvalues_list.append(depletion_pvalue)
    
    # ===== 格式化比值為字串 =====
    exp_count_rounded = round(expected_count, 2)
    exp_ratio_pct = round((expected_ratio * 100), 2)
    expected_ratio_str = f"{exp_count_rounded}/{n} ({exp_ratio_pct}%)"
    
    obs_ratio_pct = round((observed_ratio * 100), 2)
    observed_ratio_str = f"{k}/{n} ({obs_ratio_pct}%)"
    
    # 儲存結果
    results.append({
        'Feature Name': feature_name,
        'Expected Ratio': expected_ratio_str,
        'Observed Ratio': observed_ratio_str,
        'Depletion Fold Change': fold_change,
        'P-value': depletion_pvalue,
    })
    
    if (idx + 1) % 10 == 0:
        print(f"  • 已處理 {idx + 1}/{len(feature_df)} 個特徵")

print(f"  ✓ 檢定完成: 共 {len(results)} 個特徵")

# ========== 步驟5: 多重檢驗校正 ==========
print(f"\n[步驟5] 應用多重檢驗校正...")

pvalues_array = np.array([r['P-value'] for r in results])

# Benjamini-Hochberg (FDR) 校正
fdr_rejected, fdr_pvalues, _, _ = multipletests(pvalues_array, method='fdr_bh')

# Bonferroni 校正
bonf_rejected, bonf_pvalues, _, _ = multipletests(pvalues_array, method='bonferroni')

for i, result in enumerate(results):
    result['FDR depleted P-value'] = fdr_pvalues[i]
    result['Bonferroni depleted P-value'] = bonf_pvalues[i]

print(f"  ✓ FDR 校正完成")
print(f"  ✓ Bonferroni 校正完成")

# ========== 步驟6: 計算對數轉換 ==========
print(f"\n[步驟6] 計算對數轉換值...")

for result in results:
    # log2(Fold Change)
    fc = result['Depletion Fold Change']
    if fc > 0:
        result['log2(Fold Change)'] = math.log2(fc)
    else:
        result['log2(Fold Change)'] = -1000  # 處理 zero division
    
    # -log10(P-value)
    pval = result['P-value']
    if pval > 0:
        result['-log10(P-value)'] = -math.log10(pval)
    else:
        result['-log10(P-value)'] = 1000  # P-value 為 0 時
    
    # -log10(FDR P-value)
    fdr_pval = result['FDR depleted P-value']
    if fdr_pval > 0:
        result['-log10(FDR P-value)'] = -math.log10(fdr_pval)
    else:
        result['-log10(FDR P-value)'] = 1000

print(f"  ✓ 對數轉換完成")

# ========== 步驟7: 建立結果 DataFrame 並排序 ==========
print(f"\n[步驟7] 整理結果並排序...")

results_df = pd.DataFrame(results)

# 按 P-value 升序排序
results_df = results_df.sort_values('P-value', ascending=True).reset_index(drop=True)

# 重新排列欄位順序
column_order = [
    'Feature Name',
    'Expected Ratio',
    'Observed Ratio',
    'Depletion Fold Change',
    'P-value',
    'FDR depleted P-value',
    'Bonferroni depleted P-value',
    'log2(Fold Change)',
    '-log10(P-value)',
    '-log10(FDR P-value)'
]

results_df = results_df[column_order]

print(f"  ✓ 結果已排序 (按 P-value)")

# ========== 步驟8: 輸出結果 ==========
print(f"\n[步驟8] 儲存結果...")

results_df.to_csv(output_file, index=False, encoding='utf-8')
print(f"  ✓ 完整報告已儲存: {output_file}")

# ========== 摘要資訊 ==========
print("\n" + "=" * 80)
print("分析完成摘要")
print("=" * 80)
print(f"背景母體 N:              {N} 個基因")
print(f"輸入基因 n:              {n} 個基因")
print(f"特徵群總數:              {len(results_df)} 個")
print(f"\nDepletion 統計:")
print(f"  • 最小 P-value:        {results_df['P-value'].min():.4e}")
print(f"  • 最大 P-value:        {results_df['P-value'].max():.4f}")
print(f"  • p < 0.05:            {(results_df['P-value'] < 0.05).sum()} 個")
print(f"  • p < 0.01:            {(results_df['P-value'] < 0.01).sum()} 個")
print(f"  • p < 0.001:           {(results_df['P-value'] < 0.001).sum()} 個")
print(f"\nFDR 校正 (Benjamini-Hochberg):")
print(f"  • FDR < 0.05:          {(results_df['FDR depleted P-value'] < 0.05).sum()} 個")
print(f"  • FDR < 0.01:          {(results_df['FDR depleted P-value'] < 0.01).sum()} 個")
print(f"\nBonferroni 校正:")
print(f"  • Bonf < 0.05:         {(results_df['Bonferroni depleted P-value'] < 0.05).sum()} 個")
print(f"  • Bonf < 0.01:         {(results_df['Bonferroni depleted P-value'] < 0.01).sum()} 個")
print(f"\nFold Change 統計:")
print(f"  • 最小:                {results_df['Depletion Fold Change'].min():.4f}")
print(f"  • 最大:                {results_df['Depletion Fold Change'].max():.4f}")
print(f"  • 平均:                {results_df['Depletion Fold Change'].mean():.4f}")
print("=" * 80)

# 顯示前幾個最顯著的特徵
print("\n前 10 個最顯著的特徵 (按 P-value):")
print(results_df.head(10)[['Feature Name', 'Observed Ratio', 'Depletion Fold Change', 'P-value', 'FDR depleted P-value']].to_string(index=False))

print("\n✓ 完整去耗盡分析報告生成完成！")
