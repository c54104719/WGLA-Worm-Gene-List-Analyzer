#!/usr/bin/env python3
# C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d23_generate_example_input.py
"""
Generate example input file for GO enrichment analysis
從 d22_go_features_C.csv 選擇最大的 GO term，提取基因，進行隨機採樣和噪音注入

輸入: d22_go_features_C.csv
輸出: d23_example_input.txt (一行一個基因ID)

邏輯:
1. 讀取 Feature Table
2. 找到基因數最多的 GO term
3. 隨機選擇該 GO term 中 50% 的基因（上限 100 個）
4. 添加 20 個隨機噪音基因 (mock IDs)
5. 輸出為文本檔案，並顯示選中的 GO term
"""

import pandas as pd
import os
import sys
import random


def generate_mock_genes(count=20, start_id=999000):
    """
    生成隨機的 mock 基因 ID (噪音)
    
    參數:
        count: 噪音基因數量
        start_id: 起始 ID
    
    返回:
        list: Mock 基因 ID 列表
    """
    mock_genes = []
    for i in range(count):
        mock_id = f"WBGene{start_id + i}"
        mock_genes.append(mock_id)
    
    return mock_genes


def main():
    # 設定檔案路徑
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    input_file = os.path.join(script_dir, 'd22_go_features_C.csv')
    output_file = os.path.join(script_dir, 'd23_example_input.txt')
    
    print(f"開始生成範例輸入檔...")
    print(f"輸入檔案: {input_file}")
    print()
    
    # 檢查輸入檔案
    if not os.path.exists(input_file):
        print(f"✗ 錯誤: 找不到輸入檔案 {input_file}")
        sys.exit(1)
    
    # 讀取 Feature Table
    try:
        df = pd.read_csv(input_file)
        print(f"✓ 讀取成功")
        print(f"  總 GO terms: {len(df)}")
        print()
    except Exception as e:
        print(f"✗ 讀取檔案時發生錯誤: {str(e)}")
        sys.exit(1)
    
    # 檢查必要欄位
    required_cols = ['Feature Name', 'Number', 'Gene']
    if not all(col in df.columns for col in required_cols):
        print(f"✗ 錯誤: 缺少必要欄位 {required_cols}")
        print(f"  現有欄位: {list(df.columns)}")
        sys.exit(1)
    
    # 找到基因數最多的 GO term
    max_idx = df['Number'].idxmax()
    max_row = df.loc[max_idx]
    
    feature_name = max_row['Feature Name']
    go_count = max_row['Number']
    gene_string = max_row['Gene']
    
    print(f"✓ 找到最大的 GO term:")
    print(f"  Feature Name: {feature_name}")
    print(f"  GO 基因數: {go_count}")
    print()
    
    # 提取基因列表
    genes = [g.strip() for g in gene_string.split(',')]
    genes = [g for g in genes if g]  # 移除空字串
    
    print(f"  提取來自該 GO term 的基因數: {len(genes)}")
    
    # 隨機選擇 50% 的基因（上限 100 個）
    sample_size = min(len(genes) // 2, 100)
    selected_genes = random.sample(genes, sample_size)
    
    print(f"  隨機選擇 50% 的基因 (上限 100): {len(selected_genes)}")
    
    # 添加噪音基因
    noise_genes = generate_mock_genes(count=20)
    all_genes = selected_genes + noise_genes
    
    print(f"  添加噪音基因: {len(noise_genes)}")
    print(f"  總基因數 (含噪音): {len(all_genes)}")
    print()
    
    # 隨機打亂順序
    random.shuffle(all_genes)
    
    # 寫入輸出檔案
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            for gene in all_genes:
                f.write(f"{gene}\n")
        
        print(f"✓ 輸出檔案生成完成")
        print(f"  輸出檔案: {output_file}")
        print(f"  基因總數: {len(all_genes)}")
        print()
        print(f"{'='*60}")
        print(f"💡 重要提示:")
        print(f"{'='*60}")
        print(f"選中的 GO term (預期應顯示 Enriched):")
        print(f"  Feature Name: {feature_name}")
        print(f"  GO term 原始基因數: {go_count}")
        print(f"  輸入檔案中該 GO term 的基因數: {len(selected_genes)}")
        print(f"             (該 GO term 原始基因的 50%)")
        print(f"  估計 P-value: 應該很小 (顯著富集)")
        print(f"{'='*60}")
        
    except Exception as e:
        print(f"✗ 寫入檔案時發生錯誤: {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    main()
