#!/usr/bin/env python3
# C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d22_split_and_feature_go.py
"""
Split GO annotation by Aspect and Generate Feature Tables
從 d21_go_cleaned.csv 分割資料並生成Feature Table

輸入: d21_go_cleaned.csv
輸出: 三個CSV檔案
  - d22_go_features_C.csv (Cellular Component)
  - d22_go_features_P.csv (Biological Process)
  - d22_go_features_F.csv (Molecular Function)

Feature Table 結構:
  - Feature Name: GO ID (例如 GO:0005737)  <-- [修改] 現在改成用 ID
  - Number: Count of unique genes for this GO term
  - Gene: Comma-separated list of WBGene IDs
"""

import pandas as pd
import os
import sys

def load_go_names(obo_file):
    """
    Load GO term names from OBO file
    參數: obo_file: gene_ontology.*.obo 檔案路徑
    返回: dict: {GO_ID: GO_Name}
    """
    go_names = {}
    
    if not os.path.exists(obo_file):
        print(f"警告: 找不到 OBO 檔案 {obo_file}，將使用 GO ID 做為名稱")
        return go_names
    
    try:
        with open(obo_file, 'r', encoding='utf-8') as f:
            current_id = None
            for line in f:
                line = line.strip()
                if line.startswith('id: '):
                    current_id = line.split(': ')[1]
                elif line.startswith('name: ') and current_id:
                    go_name = line.split(': ', 1)[1]
                    go_names[current_id] = go_name
        
        print(f"✓ 載入 {len(go_names)} 個 GO term 名稱")
        return go_names
        
    except Exception as e:
        print(f"警告: 讀取 OBO 檔案時發生錯誤: {str(e)}")
        return go_names


def generate_feature_table(df, aspect_name, aspect_code, go_names, min_count=5):
    """
    Generate feature table from GO annotation data
    """
    print(f"\n處理 {aspect_name} ({aspect_code})...")
    print(f"  原始資料筆數: {len(df)}")
    
    # 確保有必要的欄位
    required_cols = ['WormBase ID', 'GO ID']
    if not all(col in df.columns for col in required_cols):
        print(f"  ✗ 缺少必要欄位: {required_cols}")
        return None
    
    # 分組並聚合
    # 使用 groupby('GO ID') 確保每個 ID 只有一列
    grouped = df.groupby('GO ID').agg({
        'WormBase ID': lambda ids: ','.join(ids.unique())  # 去重並用逗號連接
    }).reset_index()
    
    # 計算每個 GO term 的基因數量
    grouped['Number'] = grouped['WormBase ID'].apply(lambda x: len(x.split(',')))
    
    # ---------------- [修改區域開始] ----------------
    # 生成 Feature Name
    
    # [舊邏輯] 使用 GO Name (已註解，保留備用)
    # grouped['Feature Name'] = grouped['GO ID'].apply(
    #     lambda goid: go_names.get(goid, goid)
    # )

    # [新邏輯] 直接使用 GO ID 當作 Feature Name
    grouped['Feature Name'] = grouped['GO ID']
    # ---------------- [修改區域結束] ----------------
    
    # 重新排列和篩選欄位
    feature_table = grouped[['Feature Name', 'Number', 'WormBase ID']].copy()
    feature_table.columns = ['Feature Name', 'Number', 'Gene']
    
    # 應用最小數量過濾
    before_filter = len(feature_table)
    feature_table = feature_table[feature_table['Number'] >= min_count]
    after_filter = len(feature_table)
    
    print(f"  過濾前 GO terms: {before_filter}")
    print(f"  最小數量過濾 (>= {min_count}): {before_filter - after_filter} 筆移除")
    print(f"  過濾後 GO terms: {after_filter}")
    
    # 排序 (按基因數量降序)
    feature_table = feature_table.sort_values('Number', ascending=False).reset_index(drop=True)
    
    return feature_table


def main():
    # 設定檔案路徑
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    
    # 輸入檔 (請確認這個路徑是否正確，如果不對請自行修改)
    input_file = os.path.join(script_dir, 'd21_go_cleaned.csv')
    
    # OBO 檔 (如果不需要轉名字其實可以不用讀，但留著也無妨)
    obo_file = os.path.join(project_root, 'HW4PreprocessRawData', 'gene_ontology.WS298.obo')
    
    output_files = {
        'C': os.path.join(script_dir, 'd22_go_features_C.csv'),
        'P': os.path.join(script_dir, 'd22_go_features_P.csv'),
        'F': os.path.join(script_dir, 'd22_go_features_F.csv')
    }
    
    print(f"開始分割 GO annotation 並生成 Feature Tables...")
    print(f"輸入檔案: {input_file}")
    print()
    
    # 檢查輸入檔案
    if not os.path.exists(input_file):
        print(f"✗ 錯誤: 找不到輸入檔案 {input_file}")
        sys.exit(1)
    
    # 讀取清理後的 GO annotation 資料
    try:
        df = pd.read_csv(input_file)
        print(f"✓ 讀取成功")
        print(f"  總筆數: {len(df)}")
        print(f"  欄位: {list(df.columns)}")
        print()
    except Exception as e:
        print(f"✗ 讀取檔案時發生錯誤: {str(e)}")
        sys.exit(1)
    
    # 檢查必要欄位
    if 'Aspect' not in df.columns:
        print(f"✗ 錯誤: 檔案中找不到 'Aspect' 欄位")
        sys.exit(1)
    
    print(f"Aspect 值分佈:")
    aspect_dist = df['Aspect'].value_counts().sort_index()
    for aspect, count in aspect_dist.items():
        print(f"  {aspect}: {count}")
    print()
    
    # 載入 GO term 名稱 (雖然現在 Feature Name 不用它，但也許之後要查表)
    go_names = load_go_names(obo_file)
    print()
    
    # 分割資料並生成 Feature Tables
    aspect_map = {
        'C': ('Cellular Component', df[df['Aspect'] == 'C']),
        'P': ('Biological Process', df[df['Aspect'] == 'P']),
        'F': ('Molecular Function', df[df['Aspect'] == 'F'])
    }
    
    results = {}
    for aspect_code, (aspect_name, aspect_df) in aspect_map.items():
        if len(aspect_df) == 0:
            print(f"警告: {aspect_name} ({aspect_code}) 沒有資料")
            continue
        
        feature_table = generate_feature_table(aspect_df, aspect_name, aspect_code, go_names)
        if feature_table is not None:
            results[aspect_code] = feature_table
    
    # 寫入輸出檔案
    print(f"\n生成輸出檔案:")
    for aspect_code, output_file in output_files.items():
        if aspect_code in results:
            feature_table = results[aspect_code]
            feature_table.to_csv(output_file, index=False)
            print(f"  ✓ {os.path.basename(output_file)}: {len(feature_table)} 個 Feature")
        else:
            print(f"  ⊘ {os.path.basename(output_file)}: 跳過 (無資料)")
    
    print(f"\n✓ 完成！")


if __name__ == '__main__':
    main()