"""
反整理 (Unpivot/Transpose) 布尔值特徵，建立特徵表
將每個 True/False 欄位轉換為包含基因列表的特徵
"""

import pandas as pd

# ========== 檔案路徑設定 ==========
input_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d11_fgg_cleaned.csv"
output_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d12_fgg_55feature_table.csv"

print("=" * 70)
print("反整理特徵表生成程式")
print("=" * 70)

# ========== 步驟1: 讀取清洗後的資料 ==========
print("\n[步驟1] 讀取清洗後的資料...")
df = pd.read_csv(input_file)
print(f"  ✓ 資料維度: {df.shape[0]} 行 × {df.shape[1]} 列")

# 確認基因 ID 欄位
gene_id_column = 'row_names'
if gene_id_column not in df.columns:
    print(f"  ✗ 找不到 '{gene_id_column}' 欄位")
    print(f"  可用欄位: {list(df.columns[:5])}")
    exit(1)

print(f"  ✓ 基因 ID 欄位: {gene_id_column}")

# ========== 步驟2: 識別布爾值欄位 ==========
print("\n[步驟2] 識別布爾值 (True/False) 欄位...")

boolean_columns = []

for col in df.columns:
    # 跳過基因 ID 欄位和其他識別欄位
    if col in ['row_names', 'row_names2', 'row_names3']:
        continue
    
    # 獲取非空值
    non_null_values = df[col].dropna().unique()
    
    # 檢查是否只包含 True/False (可能是布爾型或字串型)
    # 轉換為字串後檢查
    str_values = set(str(v).upper() for v in non_null_values)
    
    # 判斷是否為布爾值欄位
    if str_values.issubset({'TRUE', 'FALSE', 'True', 'False', 'true', 'false'}):
        boolean_columns.append(col)

print(f"  ✓ 找到 {len(boolean_columns)} 個布爾值欄位")

if len(boolean_columns) == 0:
    print("  ✗ 沒有找到布爾值欄位，程式結束")
    exit(1)

# 顯示前幾個布爾值欄位
print(f"  範例欄位: {boolean_columns[:5]}")

# ========== 步驟3: 反整理 (建立特徵表) ==========
print("\n[步驟3] 反整理布爾值欄位，建立特徵表...")

feature_list = []

for col in boolean_columns:
    # 找出該欄位為 True 的所有基因 ID
    # 處理可能的布爾型或字串型
    true_mask = df[col].astype(str).str.upper() == 'TRUE'
    true_genes = df.loc[true_mask, gene_id_column].tolist()
    
    # 計算數量
    count = len(true_genes)
    
    # 將基因 ID 用逗號連接
    gene_string = ','.join(true_genes)
    
    # 添加到特徵列表
    feature_list.append({
        'Gene list': col,
        'Number': count,
        'Gene': gene_string
    })

# 建立 DataFrame
feature_df = pd.DataFrame(feature_list)

print(f"  ✓ 生成 {len(feature_df)} 個特徵")

# ========== 步驟4: 輸出結果 ==========
print("\n[步驟4] 儲存特徵表...")
feature_df.to_csv(output_file, index=False, encoding='utf-8')
print(f"  ✓ 已儲存至: {output_file}")

# ========== 摘要資訊 ==========
print("\n" + "=" * 70)
print("特徵表生成完成摘要")
print("=" * 70)
print(f"輸入資料:     {df.shape[0]} 個基因 × {df.shape[1]} 個欄位")
print(f"布爾值欄位:   {len(boolean_columns)} 個")
print(f"生成特徵:     {len(feature_df)} 個")
print("=" * 70)

# 顯示統計資訊
print("\n特徵統計:")
print(f"  最大成員數: {feature_df['Number'].max()}")
print(f"  最小成員數: {feature_df['Number'].min()}")
print(f"  平均成員數: {feature_df['Number'].mean():.2f}")

# 顯示前幾個特徵
print("\n前 5 個特徵預覽:")
print(feature_df.head().to_string(index=False))

print("\n✓ 特徵表生成完成!")
