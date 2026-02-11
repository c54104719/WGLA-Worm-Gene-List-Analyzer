"""
最終版本：完整富集與耗盡分析腳本
執行超幾何分佈檢定 (Hypergeometric Test) - Enrichment & Depletion
在單次運行中生成兩個完整報告檔案
"""

import pandas as pd
import numpy as np
import math
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# ========== 檔案路徑設定 ==========
population_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d01_live_genes_only.csv"
feature_table_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d12_fgg_55feature_table.csv"
input_genes_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d13_example_input.txt"
enrichment_output = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d14_fgg_exampleOutput_enrichment.csv"
depletion_output = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d14_fgg_exampleOutput_depletion.csv"

print("=" * 90)
print("最終版本：完整富集與耗盡分析腳本")
print("超幾何分佈檢定 (Hypergeometric Test)")
print("=" * 90)

# ========== 步驟1: 讀取背景母體 ==========
print("\n[步驟1] 讀取背景母體 (Live 基因)...")
try:
    population_df = pd.read_csv(population_file)
    N = len(population_df)  # 總基因數
    print(f"  ✓ 背景母體總數 N = {N}")
except FileNotFoundError:
    print(f"  ✗ 找不到背景母體檔案: {population_file}")
    exit(1)

# ========== 步驟2: 讀取輸入基因清單 ==========
print(f"\n[步驟2] 讀取輸入基因清單...")
try:
    with open(input_genes_file, 'r', encoding='utf-8') as f:
        input_genes = set(line.strip() for line in f if line.strip())
except FileNotFoundError:
    print(f"  ✗ 找不到輸入基因清單: {input_genes_file}")
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

# ========== 步驟4: 執行分析並計算統計量 ==========
print(f"\n[步驟4] 執行超幾何分佈檢定並計算統計量...")

enrichment_results = []
depletion_results = []
enrichment_pvalues = []
depletion_pvalues = []

for idx, row in feature_df.iterrows():
    feature_name = row['Gene list']
    M = int(row['Number'])  # 特徵群成員數
    
    # 解析基因 ID (逗號分隔，去除空白)
    gene_list_str = row['Gene']
    if pd.isna(gene_list_str) or gene_list_str == '':
        continue
    
    feature_genes = set(g.strip() for g in str(gene_list_str).split(',') if g.strip())
    
    # 計算交集
    overlap_genes = input_genes & feature_genes
    k = len(overlap_genes)
    
    # ===== 計算比值和 Fold Change =====
    expected_count = (n * M) / N
    expected_ratio = M / N if N > 0 else 0
    observed_ratio = k / n if n > 0 else 0
    
    if expected_ratio > 0:
        fold_change = observed_ratio / expected_ratio
    else:
        fold_change = 0
    
    # ===== 格式化比值為字串 =====
    exp_count_rounded = round(expected_count, 2)
    exp_ratio_pct = round((expected_ratio * 100), 2)
    expected_ratio_str = f"{exp_count_rounded}/{n} ({exp_ratio_pct}%)"
    
    obs_ratio_pct = round((observed_ratio * 100), 2)
    observed_ratio_str = f"{k}/{n} ({obs_ratio_pct}%)"
    
    # ===== P-value (Enrichment) =====
    enrichment_pvalue = hypergeom.sf(k - 1, N, M, n)
    enrichment_pvalues.append(enrichment_pvalue)
    
    # ===== P-value (Depletion) =====
    depletion_pvalue = hypergeom.cdf(k, N, M, n)
    depletion_pvalues.append(depletion_pvalue)
    
    # 儲存 Enrichment 結果
    enrichment_results.append({
        'Feature Name': feature_name,
        'Observed Ratio': observed_ratio_str,
        'Expected Ratio': expected_ratio_str,
        'Fold Change': fold_change,
        'P-value': enrichment_pvalue,
    })
    
    # 儲存 Depletion 結果
    depletion_results.append({
        'Feature Name': feature_name,
        'Observed Ratio': observed_ratio_str,
        'Expected Ratio': expected_ratio_str,
        'Fold Change': fold_change,
        'P-value': depletion_pvalue,
    })
    
    if (idx + 1) % 10 == 0:
        print(f"  • 已處理 {idx + 1}/{len(feature_df)} 個特徵")

print(f"  ✓ 檢定完成: 共 {len(enrichment_results)} 個特徵")

# ========== 步驟5: 多重檢驗校正 ==========
print(f"\n[步驟5] 應用多重檢驗校正...")

# Enrichment 校正
enrichment_pvalues_array = np.array(enrichment_pvalues)
erich_fdr_rejected, erich_fdr_pvalues, _, _ = multipletests(enrichment_pvalues_array, method='fdr_bh')
erich_bonf_rejected, erich_bonf_pvalues, _, _ = multipletests(enrichment_pvalues_array, method='bonferroni')

# Depletion 校正
depletion_pvalues_array = np.array(depletion_pvalues)
depl_fdr_rejected, depl_fdr_pvalues, _, _ = multipletests(depletion_pvalues_array, method='fdr_bh')
depl_bonf_rejected, depl_bonf_pvalues, _, _ = multipletests(depletion_pvalues_array, method='bonferroni')

for i in range(len(enrichment_results)):
    enrichment_results[i]['FDR P-value'] = erich_fdr_pvalues[i]
    enrichment_results[i]['Bonferroni P-value'] = erich_bonf_pvalues[i]
    
    depletion_results[i]['FDR P-value'] = depl_fdr_pvalues[i]
    depletion_results[i]['Bonferroni P-value'] = depl_bonf_pvalues[i]

print(f"  ✓ FDR 與 Bonferroni 校正完成")

# ========== 步驟6: 計算對數轉換 ==========
print(f"\n[步驟6] 計算對數轉換值...")

for result in enrichment_results:
    # log2(Fold Change)
    fc = result['Fold Change']
    if fc > 0:
        result['log2(Fold Change)'] = math.log2(fc)
    else:
        result['log2(Fold Change)'] = -1000
    
    # -log10(P-value)
    pval = result['P-value']
    if pval > 0:
        result['-log10(P-value)'] = -math.log10(pval)
    else:
        result['-log10(P-value)'] = 1000
    
    # -log10(FDR P-value)
    fdr_pval = result['FDR P-value']
    if fdr_pval > 0:
        result['-log10(FDR P-value)'] = -math.log10(fdr_pval)
    else:
        result['-log10(FDR P-value)'] = 1000

for result in depletion_results:
    # log2(Fold Change)
    fc = result['Fold Change']
    if fc > 0:
        result['log2(Fold Change)'] = math.log2(fc)
    else:
        result['log2(Fold Change)'] = -1000
    
    # -log10(P-value)
    pval = result['P-value']
    if pval > 0:
        result['-log10(P-value)'] = -math.log10(pval)
    else:
        result['-log10(P-value)'] = 1000
    
    # -log10(FDR P-value)
    fdr_pval = result['FDR P-value']
    if fdr_pval > 0:
        result['-log10(FDR P-value)'] = -math.log10(fdr_pval)
    else:
        result['-log10(FDR P-value)'] = 1000

print(f"  ✓ 對數轉換完成")

# ========== 步驟7: 建立 DataFrame 並排序 ==========
print(f"\n[步驟7] 整理結果並排序...")

enrichment_df = pd.DataFrame(enrichment_results)
enrichment_df = enrichment_df.sort_values('P-value', ascending=True).reset_index(drop=True)

depletion_df = pd.DataFrame(depletion_results)
depletion_df = depletion_df.sort_values('P-value', ascending=True).reset_index(drop=True)

# 重新排列欄位順序
column_order = [
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

enrichment_df = enrichment_df[column_order]
depletion_df = depletion_df[column_order]

print(f"  ✓ 結果已排序 (按 P-value)")

# ========== 步驟8: 輸出結果 ==========
print(f"\n[步驟8] 儲存結果檔案...")

enrichment_df.to_csv(enrichment_output, index=False, encoding='utf-8')
print(f"  ✓ Enrichment 報告已儲存: {enrichment_output}")

depletion_df.to_csv(depletion_output, index=False, encoding='utf-8')
print(f"  ✓ Depletion 報告已儲存: {depletion_output}")

# ========== 摘要資訊 ==========
print("\n" + "=" * 90)
print("分析完成摘要")
print("=" * 90)
print(f"背景母體 N:              {N} 個基因")
print(f"輸入基因 n:              {n} 個基因")
print(f"特徵群總數:              {len(enrichment_df)} 個")

print(f"\n【Enrichment 統計】")
print(f"  • 最小 P-value:        {enrichment_df['P-value'].min():.4e}")
print(f"  • 最大 P-value:        {enrichment_df['P-value'].max():.4f}")
print(f"  • p < 0.05:            {(enrichment_df['P-value'] < 0.05).sum()} 個")
print(f"  • p < 0.01:            {(enrichment_df['P-value'] < 0.01).sum()} 個")
print(f"  • p < 0.001:           {(enrichment_df['P-value'] < 0.001).sum()} 個")
print(f"  • FDR < 0.05:          {(enrichment_df['FDR P-value'] < 0.05).sum()} 個")
print(f"  • Bonf < 0.05:         {(enrichment_df['Bonferroni P-value'] < 0.05).sum()} 個")
print(f"  • Fold Change 平均:    {enrichment_df['Fold Change'].mean():.4f}")

print(f"\n【Depletion 統計】")
print(f"  • 最小 P-value:        {depletion_df['P-value'].min():.4e}")
print(f"  • 最大 P-value:        {depletion_df['P-value'].max():.4f}")
print(f"  • p < 0.05:            {(depletion_df['P-value'] < 0.05).sum()} 個")
print(f"  • p < 0.01:            {(depletion_df['P-value'] < 0.01).sum()} 個")
print(f"  • p < 0.001:           {(depletion_df['P-value'] < 0.001).sum()} 個")
print(f"  • FDR < 0.05:          {(depletion_df['FDR P-value'] < 0.05).sum()} 個")
print(f"  • Bonf < 0.05:         {(depletion_df['Bonferroni P-value'] < 0.05).sum()} 個")
print(f"  • Fold Change 平均:    {depletion_df['Fold Change'].mean():.4f}")

print("=" * 90)

# 顯示前幾個最顯著的特徵
print("\n【Enrichment - 前 10 個最顯著特徵 (按 P-value)】")
print(enrichment_df.head(10)[['Feature Name', 'Observed Ratio', 'Fold Change', 'P-value']].to_string(index=False))

print("\n【Depletion - 前 10 個最顯著特徵 (按 P-value)】")
print(depletion_df.head(10)[['Feature Name', 'Observed Ratio', 'Fold Change', 'P-value']].to_string(index=False))

print("\n✓ 分析完成！已生成 2 個報告檔案")
print(f"  • {enrichment_output}")
print(f"  • {depletion_output}")
