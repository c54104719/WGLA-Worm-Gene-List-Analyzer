"""
超幾何分佈檢定 (Hypergeometric Test)
針對功能基因群(FGG)進行 Enrichment/Depletion 分析
"""

import pandas as pd
from scipy.stats import hypergeom
import numpy as np

# ========== 檔案路徑設定 ==========
# 背景母體 (Population): 顯示49164個Live基因的清單
population_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d01_live_genes_only.csv"

# 特徵表
feature_table_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d12_fgg_feature_table.csv"

# 輸入基因清單
input_genes_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d13_example_input.txt"

# 輸出檔案
enrichment_output = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d14_fgg_enrichment.csv"
depletion_output = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d14_fgg_depletion.csv"

print("=" * 70)
print("超幾何分佈檢定 (Hypergeometric Test)")
print("功能基因群 (FGG) 富集分析")
print("=" * 70)

# ========== 步驟1: 讀取背景母體 ==========
print("\n[步驟1] 讀取背景母體 (Live 基因)...")
population_df = pd.read_csv(population_file)
N = len(population_df)  # 總基因數
print(f"  ✓ 背景母體總數 N = {N}")

# 建立基因 ID 集合
population_genes = set(population_df.iloc[:, 0])  # 第一欄是基因 ID

# ========== 步驟2: 讀取輸入基因清單 ==========
print(f"\n[步驟2] 讀取輸入基因清單 ({input_genes_file})...")
try:
    with open(input_genes_file, 'r', encoding='utf-8') as f:
        input_genes = set(line.strip() for line in f if line.strip())
except FileNotFoundError:
    print(f"  ✗ 找不到輸入檔案: {input_genes_file}")
    print("  請確保 d13_example.txt 存在於 HW42 資料夾")
    exit(1)

n = len(input_genes)
print(f"  ✓ 輸入基因數 n = {n}")

# 計算真實輸入數 (只保留在母體中的基因)
valid_input_genes = input_genes & population_genes
valid_n = len(valid_input_genes)
print(f"  ✓ 有效輸入基因 (在母體中) = {valid_n}")

if valid_n == 0:
    print("  ✗ 錯誤: 沒有輸入基因在背景母體中")
    exit(1)

# ========== 步驟3: 讀取特徵表 ==========
print(f"\n[步驟3] 讀取特徵表 ({feature_table_file})...")
try:
    feature_df = pd.read_csv(feature_table_file)
except FileNotFoundError:
    print(f"  ✗ 找不到特徵表: {feature_table_file}")
    print("  請確保 d12_fgg_feature_table.csv 存在於 HW42 資料夾")
    exit(1)

print(f"  ✓ 特徵數量: {len(feature_df)}")
print(f"  ✓ 欄位: {list(feature_df.columns)}")

# ========== 步驟4: 執行超幾何分佈檢定 ==========
print(f"\n[步驟4] 執行超幾何分佈檢定...")

enrichment_results = []
depletion_results = []

for idx, row in feature_df.iterrows():
    feature_name = row['Gene list']
    M = int(row['Number'])  # 特徵群成員數
    
    # 解析基因 ID (逗號分隔)
    gene_list_str = row['Gene']
    if pd.isna(gene_list_str) or gene_list_str == '':
        continue
    
    feature_genes = set(g.strip() for g in str(gene_list_str).split(','))
    
    # 計算交集
    overlap_genes = valid_input_genes & feature_genes
    k = len(overlap_genes)
    
    # 計算期望值 (Expected Count)
    expected_count = (n * M) / N
    
    # 超幾何分佈參數
    # N: 總體大小
    # M: 總體中成功狀態的數量
    # n: 抽樣大小
    # k: 抽樣中成功的數量
    
    # Enrichment p-value: P(X >= k) = sf(k-1)
    # sf(k-1) 等同於 1 - cdf(k-1)
    enrichment_pvalue = hypergeom.sf(k - 1, N, M, n)
    
    # Depletion p-value: P(X <= k) = cdf(k)
    depletion_pvalue = hypergeom.cdf(k, N, M, n)
    
    # 儲存 Enrichment 結果
    enrichment_results.append({
        'Feature Name': feature_name,
        'Overlap Count': k,
        'Expected Count': f"{expected_count:.4f}",
        'p-value': enrichment_pvalue
    })
    
    # 儲存 Depletion 結果
    depletion_results.append({
        'Feature Name': feature_name,
        'Overlap Count': k,
        'Expected Count': f"{expected_count:.4f}",
        'p-value': depletion_pvalue
    })
    
    if (idx + 1) % 10 == 0:
        print(f"  • 已處理 {idx + 1}/{len(feature_df)} 個特徵")

print(f"  ✓ 檢定完成")

# ========== 步驟5: 建立結果 DataFrame 並排序 ==========
print(f"\n[步驟5] 整理結果...")

enrichment_results_df = pd.DataFrame(enrichment_results)
enrichment_results_df = enrichment_results_df.sort_values('p-value')

depletion_results_df = pd.DataFrame(depletion_results)
depletion_results_df = depletion_results_df.sort_values('p-value')

print(f"  ✓ Enrichment 結果: {len(enrichment_results_df)} 個特徵")
print(f"  ✓ Depletion 結果: {len(depletion_results_df)} 個特徵")

# ========== 步驟6: 輸出結果 ==========
print(f"\n[步驟6] 儲存結果...")

enrichment_results_df.to_csv(enrichment_output, index=False, encoding='utf-8')
print(f"  ✓ Enrichment 結果: {enrichment_output}")

depletion_results_df.to_csv(depletion_output, index=False, encoding='utf-8')
print(f"  ✓ Depletion 結果: {depletion_output}")

# ========== 摘要資訊 ==========
print("\n" + "=" * 70)
print("分析完成摘要")
print("=" * 70)
print(f"背景母體 N:           {N} 個基因")
print(f"輸入基因 n:           {valid_n} 個基因")
print(f"特徵群總數:           {len(feature_df)} 個")
print(f"\nEnrichment 結果:")
print(f"  • 最小 p-value:     {enrichment_results_df['p-value'].min():.4e}")
print(f"  • 最大 p-value:     {enrichment_results_df['p-value'].max():.4f}")
print(f"  • p < 0.05:         {(enrichment_results_df['p-value'] < 0.05).sum()} 個")
print(f"  • p < 0.01:         {(enrichment_results_df['p-value'] < 0.01).sum()} 個")
print(f"\nDepletion 結果:")
print(f"  • 最小 p-value:     {depletion_results_df['p-value'].min():.4e}")
print(f"  • 最大 p-value:     {depletion_results_df['p-value'].max():.4f}")
print(f"  • p < 0.05:         {(depletion_results_df['p-value'] < 0.05).sum()} 個")
print(f"  • p < 0.01:         {(depletion_results_df['p-value'] < 0.01).sum()} 個")
print("=" * 70)

# 顯示前幾個顯著的特徵
print("\nEnrichment 前 5 個顯著特徵:")
print(enrichment_results_df.head().to_string(index=False))

print("\nDepletion 前 5 個顯著特徵:")
print(depletion_results_df.head().to_string(index=False))

print("\n✓ 分析完成!")
