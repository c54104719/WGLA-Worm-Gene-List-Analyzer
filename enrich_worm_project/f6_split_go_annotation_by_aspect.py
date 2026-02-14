# C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\enrich_worm_project\f6_split_go_annotation_by_aspect.py
#!/usr/bin/env python3
"""
Split go annotation file by Aspect (C, P, F)
輸入: f5_go_annotation_with_terms.csv
輸出: 三個CSV檔案
  - f6_go_annotation_aspect_C.csv (Cellular Component)
  - f6_go_annotation_aspect_P.csv (Biological Process)
  - f6_go_annotation_aspect_F.csv (Molecular Function)
"""

import pandas as pd
import os

# 設定檔案路徑
input_file = 'f5_go_annotation_with_terms.csv'
output_dir = '.'

# 檢查輸入檔案是否存在
if not os.path.exists(input_file):
    print(f"錯誤: 找不到檔案 {input_file}")
    exit(1)

# 讀取CSV檔案
print(f"正在讀取 {input_file}...")
df = pd.read_csv(input_file)

print(f"總共讀取 {len(df)} 筆資料")
print(f"欄位: {list(df.columns)}")

# 檢查Aspect欄位
if 'Aspect' not in df.columns:
    print("錯誤: 檔案中找不到 'Aspect' 欄位")
    exit(1)

print(f"\nAspect值的分佈:")
aspect_counts = df['Aspect'].value_counts().sort_index()
print(aspect_counts)

# 根據Aspect分割資料
df_C = df[df['Aspect'] == 'C']
df_P = df[df['Aspect'] == 'P']
df_F = df[df['Aspect'] == 'F']

# 輸出三個CSV檔案
output_files = {
    'f6_go_annotation_aspect_C.csv': (df_C, 'C (Cellular Component)'),
    'f6_go_annotation_aspect_P.csv': (df_P, 'P (Biological Process)'),
    'f6_go_annotation_aspect_F.csv': (df_F, 'F (Molecular Function)')
}

print("\n正在生成輸出檔案:")
for output_file, (data, aspect_name) in output_files.items():
    output_path = os.path.join(output_dir, output_file)
    data.to_csv(output_path, index=False)
    print(f"  {output_file}: {len(data)} 筆資料 ({aspect_name})")

print("\n完成！")
print(f"\n檔案已儲存於: {os.path.abspath(output_dir)}")
