#!/usr/bin/env python3
# C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW4PreprocessRawData\f4_preprocess_go_annotation.py
"""
前處理 WormBase GO annotation 檔案
從 gene_association.WS298.wb.c_elegans 提取所需欄位，輸出為 CSV

輸入: gene_association.WS298.wb.c_elegans (GAF format)
輸出: f4_go_annotation_processed.csv

保留欄位:
- WormBase ID
- Symbol
- GO ID
- Reference
- Evidence Code
- Aspect
"""

import csv
import sys

def preprocess_gaf_file(input_file, output_file):
    """
    讀取 GAF 檔案並提取所需欄位，去除重複行
    
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
    """
    
    rows_read = 0
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
                
                # 提取所需的 6 個欄位
                wormbase_id = row[1]      # WormBase ID
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
        
        duplicates_removed = rows_read - 7 - rows_written  # 7 是註釋行
        print(f"✓ 前處理完成（已去重複）")
        print(f"  讀取行數: {rows_read}")
        print(f"  去重前: {rows_read - 7}")
        print(f"  去重後: {rows_written}")
        print(f"  移除重複數: {duplicates_removed}")
        print(f"  輸出檔案: {output_file}")
        
    except FileNotFoundError:
        print(f"✗ 錯誤: 找不到輸入檔案 {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"✗ 錯誤: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    # 設定輸入輸出檔案
    input_file = 'gene_association.WS298.wb.c_elegans'
    output_file = 'f4_go_annotation_processed.csv'
    
    print(f"開始前處理...")
    print(f"輸入: {input_file}")
    print()
    
    preprocess_gaf_file(input_file, output_file)
