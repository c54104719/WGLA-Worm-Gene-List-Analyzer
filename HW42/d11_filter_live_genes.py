"""
筛选只保留 Live 状态的基因数据
输入: c_elegans.PRJNA13758.WS298.geneOtherIDs.txt
输出: d11_live_genes_only.csv
"""

import csv

# 输入文件路径
input_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW2GO_term\c_elegans.PRJNA13758.WS298.geneOtherIDs.txt"

# 输出文件路径
output_file = r"C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW42\d11_live_genes_only.csv"

# 读取输入文件并筛选 Live 基因
live_genes = []

print(f"正在读取文件: {input_file}")

with open(input_file, 'r', encoding='utf-8') as f:
    for line in f:
        # 以 Tab 分隔
        columns = line.strip().split('\t')
        
        # 检查第二列是否为 "Live"
        if len(columns) >= 2 and columns[1] == "Live":
            live_genes.append(columns)

print(f"总共找到 {len(live_genes)} 个 Live 基因")

# 写入输出文件
print(f"正在写入文件: {output_file}")

with open(output_file, 'w', encoding='utf-8', newline='') as f:
    writer = csv.writer(f)
    
    # 写入表头
    writer.writerow(['WBGene_ID', 'Status', 'Sequence_Name', 'Gene_Name', 'CELE_ID'])
    
    # 写入所有 Live 基因数据
    writer.writerows(live_genes)

print(f"完成！已保存 {len(live_genes)} 个 Live 基因到 {output_file}")
