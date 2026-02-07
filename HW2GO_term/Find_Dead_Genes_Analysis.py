import csv
import os

def analyze_non_live_genes(gaf_path, background_path, output_report):
    # 1. 讀取母體清單，但這次我們關注 "非 Live" 的基因
    dead_genes = set()      # 存 ID
    dead_gene_info = {}     # 存 ID -> (Status, Symbol) 的對照，方便查看
    
    live_count = 0
    
    print(f"[1/3] 正在讀取母體清單: {background_path} ...")
    
    with open(background_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            columns = line.split('\t')
            if len(columns) >= 2:
                gene_id = columns[0]
                status = columns[1]
                symbol = columns[2] if len(columns) > 2 else "N/A" # 有些檔案有第三欄 Symbol
                
                if status == 'Live':
                    live_count += 1
                else:
                    # 抓出所有非 Live (Dead, Suppressed, etc.)
                    dead_genes.add(gene_id)
                    dead_gene_info[gene_id] = f"{status} | {symbol}"
    
    print(f"      - Live 基因數: {live_count}")
    print(f"      - 非 Live 基因數 (目標): {len(dead_genes)}")

    # 2. 掃描 GAF 檔，看看這些 Dead Genes 有沒有 GO 特徵
    print(f"[2/3] 正在掃描 GAF 檔案: {gaf_path} ...")
    
    # 記錄：哪些 Dead Gene 竟然有 GO 特徵？
    # 格式: {Dead_Gene_ID: [GO_ID_1, GO_ID_2...]}
    dead_genes_with_go = {}
    
    with open(gaf_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('!'):
                continue
            
            columns = line.strip().split('\t')
            if len(columns) < 5:
                continue

            gene_id = columns[1] # WBGeneID
            go_id = columns[4]   # GO:xxxxxxx
            
            # 檢查：這個 GAF 裡的基因，是否在我們的 "非 Live 名單" 中？
            if gene_id in dead_genes:
                if gene_id not in dead_genes_with_go:
                    dead_genes_with_go[gene_id] = set()
                dead_genes_with_go[gene_id].add(go_id)

    # 3. 輸出報告
    print(f"[3/3] 分析完成！正在寫入報告: {output_report}")
    
    with open(output_report, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        
        # 寫入檔頭
        writer.writerow(['Gene_ID', 'Status_and_Symbol', 'Has_GO_Term?', 'GO_Term_Count', 'GO_IDs'])
        
        # 為了報告完整，我們列出所有發現的 Dead Genes
        # 優先列出「有 GO 特徵」的，比較有研究價值
        sorted_dead_ids = sorted(list(dead_genes), key=lambda x: (x not in dead_genes_with_go, x))
        
        count_dead_with_go = 0
        
        for gene_id in sorted_dead_ids:
            info = dead_gene_info.get(gene_id, "Unknown")
            
            if gene_id in dead_genes_with_go:
                has_go = "Yes"
                go_list = dead_genes_with_go[gene_id]
                go_count = len(go_list)
                go_ids_str = ",".join(sorted(list(go_list)))
                count_dead_with_go += 1
            else:
                has_go = "No"
                go_count = 0
                go_ids_str = ""
            
            writer.writerow([gene_id, info, has_go, go_count, go_ids_str])

    print("-" * 30)
    print(f"總結報告：")
    print(f"1. 非 Live 基因總數: {len(dead_genes)}")
    print(f"2. 其中擁有 GO 特徵的數量: {count_dead_with_go} (這些就是你之前程式沒抓到的)")
    print(f"3. 詳細清單請查看: {output_report}")
    print("-" * 30)

# --- 執行設定 ---
# 請確認檔名與你電腦的一致
input_gaf_file = 'gene_association.WS298.wb.c_elegans'       
background_list_file = 'c_elegans.PRJNA13758.WS298.geneOtherIDs.txt' 
output_report_file = 'Non_Live_Genes_Report.csv' 

try:
    analyze_non_live_genes(input_gaf_file, background_list_file, output_report_file)
except FileNotFoundError as e:
    print(f"錯誤：找不到檔案 - {e}")