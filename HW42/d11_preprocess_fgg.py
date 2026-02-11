"""
清洗 GLH4.master.21_23G.csv 資料
步驟1: 篩選 Live 基因 (Row Filtering)
步驟2: 移除數值型欄位 (Column Filtering)
"""

import pandas as pd

# ========== 檔案路徑設定 ==========
whitelist_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d01_live_genes_only.csv"
raw_data_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\GLH4.master.21_23G.csv"
output_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d11_fgg_cleaned.csv"

print("=" * 60)
print("資料清洗程式 - GLH4.master.21_23G.csv")
print("=" * 60)

# ========== 步驟1: 讀取白名單 (Live 基因) ==========
print("\n[步驟1] 讀取 Live 基因白名單...")
whitelist_df = pd.read_csv(whitelist_file)
live_gene_ids = set(whitelist_df['WBGene_ID'])
print(f"  ✓ 白名單中共有 {len(live_gene_ids)} 個 Live 基因")

# ========== 步驟2: 讀取原始資料 ==========
print("\n[步驟2] 讀取原始資料...")
raw_df = pd.read_csv(raw_data_file)
print(f"  ✓ 原始資料: {raw_df.shape[0]} 行 × {raw_df.shape[1]} 列")

# ========== 步驟3: 篩選 Live 基因 (Row Filtering) ==========
print("\n[步驟3] 篩選 Live 基因 (保留在白名單中的基因)...")
# 使用 row_names 欄位作為基因 ID
filtered_df = raw_df[raw_df['row_names'].isin(live_gene_ids)]
print(f"  ✓ 篩選後: {filtered_df.shape[0]} 行 × {filtered_df.shape[1]} 列")
print(f"  ✓ 移除了 {raw_df.shape[0] - filtered_df.shape[0]} 行非 Live 基因")

# ========== 步驟4: 移除數值型欄位 (Column Filtering) ==========
print("\n[步驟4] 檢查並移除數值型欄位...")

# 找出需要保留的欄位 (非全數值型)
columns_to_keep = []
numeric_columns_removed = []

for col in filtered_df.columns:
    # 檢查該欄位是否全部都是浮點數
    # 先排除 NaN 值，然後檢查剩餘的值
    non_null_values = filtered_df[col].dropna()
    
    if len(non_null_values) == 0:
        # 如果整欄都是空值，保留它
        columns_to_keep.append(col)
        continue
    
    # 檢查是否全部都是數值型 (float 或 int)
    is_all_numeric = all(isinstance(val, (int, float)) and not isinstance(val, bool) 
                         for val in non_null_values)
    
    if is_all_numeric:
        # 如果全部是數值，移除該欄位
        numeric_columns_removed.append(col)
    else:
        # 保留非全數值的欄位 (包含 TRUE/FALSE, 字串等)
        columns_to_keep.append(col)

# 只保留非數值型欄位
cleaned_df = filtered_df[columns_to_keep]

print(f"  ✓ 保留欄位: {len(columns_to_keep)} 個")
print(f"  ✓ 移除數值型欄位: {len(numeric_columns_removed)} 個")
print(f"  ✓ 最終資料: {cleaned_df.shape[0]} 行 × {cleaned_df.shape[1]} 列")

# ========== 步驟5: 輸出結果 ==========
print("\n[步驟5] 儲存清洗後的資料...")
cleaned_df.to_csv(output_file, index=False, encoding='utf-8')
print(f"  ✓ 已儲存至: {output_file}")

# ========== 摘要資訊 ==========
print("\n" + "=" * 60)
print("清洗完成摘要")
print("=" * 60)
print(f"原始資料:     {raw_df.shape[0]} 行 × {raw_df.shape[1]} 列")
print(f"清洗後資料:   {cleaned_df.shape[0]} 行 × {cleaned_df.shape[1]} 列")
print(f"移除基因數:   {raw_df.shape[0] - cleaned_df.shape[0]} 行")
print(f"移除欄位數:   {len(numeric_columns_removed)} 列")
print("=" * 60)

# 顯示保留的欄位名稱
print("\n保留的欄位清單:")
for i, col in enumerate(columns_to_keep, 1):
    print(f"  {i}. {col}")

print("\n✓ 資料清洗完成!")
