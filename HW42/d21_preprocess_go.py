#!/usr/bin/env python3
# C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d21_preprocess_go.py
"""
前處理 WormBase GO annotation 檔案（含白名單過濾）
從 gene_association.WS298.wb.c_elegans 提取所需欄位，輸出為 CSV

輸入: 
  - gene_association.WS298.wb.c_elegans (GAF format)
  - d01_live_genes_only.csv (whitelist)
輸出: d21_go_cleaned.csv

保留欄位:
- WormBase ID
- Symbol
- GO ID
- Reference
- Evidence Code
- Aspect

新增功能:
- 白名單過濾: 只保留在 d01_live_genes_only.csv 中的基因
"""

import csv
import sys
import os

def load_whitelist(whitelist_file):
    """
    載入白名單基因 ID
    
    參數:
        whitelist_file: 白名單檔案路徑
    
    返回:
        set: 包含所有有效 WBGene ID 的集合
    """
    whitelist = set()
    
    try:
        with open(whitelist_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                wbgene_id = row['WBGene_ID'].strip()
                if wbgene_id:  # 確保不是空字串
                    whitelist.add(wbgene_id)
        
        print(f"✓ 白名單載入完成")
        print(f"  白名單基因數: {len(whitelist)}")
        print()
        
        return whitelist
        
    except FileNotFoundError:
        print(f"✗ 錯誤: 找不到白名單檔案 {whitelist_file}")
        sys.exit(1)
    except Exception as e:
        print(f"✗ 載入白名單時發生錯誤: {str(e)}")
        sys.exit(1)

def preprocess_gaf_file(input_file, output_file, whitelist):
    """
    讀取 GAF 檔案並提取所需欄位，去除重複行，並依白名單過濾
    
    GAF 格式欄位索引:
    0: DB
    1: DB Object ID (WormBase ID)
    2: DB Object Symbol
    3: Qualifier
    4: GO ID
    5: DB Reference
    6: Evidence Code
    7: With (or) From
    8: Aspect
    9+: 其他欄位
    
    參數:
        input_file: 輸入 GAF 檔案路徑
        output_file: 輸出 CSV 檔案路徑
        whitelist: 白名單基因 ID 集合
    """
    
    rows_read = 0
    rows_skipped_whitelist = 0
    unique_rows = set()  # 用來儲存唯一行
    
    try:
        with open(input_file, 'r', encoding='utf-8') as infile:
            # 讀取 GAF 檔案（tab 分隔）
            reader = csv.reader(infile, delimiter='\t')
            
            for row in reader:
                rows_read += 1
                
                # 跳過註釋行（以 ! 開頭）
                if row and row[0].startswith('!'):
                    continue
                
                # 確保行有足夠的欄位
                if len(row) < 9:
                    print(f"警告: 第 {rows_read} 行欄位不足，已跳過")
                    continue
                
                # 提取 WormBase ID
                wormbase_id = row[1]
                
                # 白名單過濾：如果 ID 不在白名單中，跳過此行
                if wormbase_id not in whitelist:
                    rows_skipped_whitelist += 1
                    continue
                
                # 提取所需的 6 個欄位
                symbol = row[2]            # Symbol
                go_id = row[4]             # GO ID
                reference = row[5]         # Reference
                evidence_code = row[6]     # Evidence Code
                aspect = row[8]            # Aspect
                
                # 製成 tuple 加入 set（自動去重）
                row_tuple = (wormbase_id, symbol, go_id, reference, evidence_code, aspect)
                unique_rows.add(row_tuple)
        
        # 寫入輸出檔案
        rows_written = 0
        with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
            writer = csv.writer(outfile)
            
            # 寫入表頭
            writer.writerow(['WormBase ID', 'Symbol', 'GO ID', 'Reference', 'Evidence Code', 'Aspect'])
            
            # 將去重後的行寫入檔案（排序以確保一致性）
            for row_tuple in sorted(unique_rows):
                writer.writerow(row_tuple)
                rows_written += 1
        
        print(f"✓ 前處理完成（已去重複 + 白名單過濾）")
        print(f"  讀取行數: {rows_read}")
        print(f"  白名單過濾移除: {rows_skipped_whitelist}")
        print(f"  去重後保留: {rows_written}")
        print(f"  輸出檔案: {output_file}")
        
    except FileNotFoundError:
        print(f"✗ 錯誤: 找不到輸入檔案 {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"✗ 錯誤: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    # 取得腳本所在目錄
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)  # enrich_project 目錄
    
    # 設定檔案路徑（使用絕對路徑）
    whitelist_file = os.path.join(script_dir, 'd01_live_genes_only.csv')
    input_file = os.path.join(script_dir, 'gene_association.WS298.wb.c_elegans')
    output_file = os.path.join(script_dir, 'd21_go_cleaned.csv')
    
    print(f"開始前處理（含白名單過濾）...")
    print(f"白名單檔案: {whitelist_file}")
    print(f"輸入檔案: {input_file}")
    print()
    
    # 載入白名單
    whitelist = load_whitelist(whitelist_file)
    
    # 前處理 GAF 檔案
    preprocess_gaf_file(input_file, output_file, whitelist)
