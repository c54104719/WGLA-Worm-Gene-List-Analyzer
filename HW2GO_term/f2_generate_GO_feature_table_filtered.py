import csv
import os

def load_background_genes(filepath):
    """
    讀取母體清單 (c_elegans.PRJNA13758.WS298.geneOtherIDs.txt)。
    邏輯修改：
    1. 使用 Tab (\t) 切割欄位。
    2. 檢查第二欄 (Index 1) 是否為 'Live'。
    3. 只有 'Live' 的基因 ID (第一欄) 才會被加入。
    """
    valid_genes = set()
    dead_count = 0
    
    if not os.path.exists(filepath):
        print(f"[警告] 找不到母體清單檔案 '{filepath}'，將跳過過濾步驟 (所有基因都會被納入)。")
        return None

    print(f"[進度] 正在讀取母體清單: {filepath} ...")

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # 根據你提供的範例，檔案是用 Tab 分隔的
            columns = line.split('\t')
            
            # 確保至少有兩欄 (ID 和 Status)
            if len(columns) >= 2:
                gene_id = columns[0] # 例如 WBGene00306080
                status = columns[1]  # 例如 Live 或 Dead
                
                # 關鍵過濾條件：只抓取 Live 的基因
                if status == 'Live':
                    valid_genes.add(gene_id)
                else:
                    dead_count += 1
            
    print(f"[資訊] 成功讀取母體清單。")
    print(f"       - 有效母體 (Live): {len(valid_genes)} 個 (這就是你的 D 值)")
    print(f"       - 排除基因 (Dead/Other): {dead_count} 個")
    
    return valid_genes

def process_gaf_to_feature_table(gaf_path, background_path, output_path):
    # 1. 讀取並篩選母體清單 (只留 Live)
    valid_genes = load_background_genes(background_path) # 這裡的background是txt檔，裡面有Live/Dead的資訊
    
    # 用來儲存 GO 特徵：Key = GO_ID, Value = Set of WormBase_IDs
    go_features = {}
    
    # 用來統計被過濾掉的例外 (不在 Live 名單中的)
    excluded_count = 0
    excluded_examples = []

    print(f"[進度] 正在讀取 GAF 原始檔案: {gaf_path} ...") # 這裡的gaf_path是從WormBase下載的gene_association.WS298.wb.c_elegans檔案，裡面包含了GO ID和對應的WormBase ID
    
    with open(gaf_path, 'r', encoding='utf-8') as f:
        for line in f:
            # 跳過註解行
            if line.startswith('!'):
                continue
            
            columns = line.strip().split('\t')
            # 確保欄位足夠
            if len(columns) < 5:
                continue

            # --- 關鍵欄位抓取 (GAF 2.1) ---
            # Col 2 (Index 1) = DB_Object_ID (WBGeneID)
            # Col 5 (Index 4) = GO ID
            
            gene_id = columns[1]      
            go_id = columns[4]
            
            # --- 母體過濾機制 ---
            # 確認這個基因是否在 "Live" 名單中
            if valid_genes is not None:
                if gene_id not in valid_genes:
                    excluded_count += 1
                    if len(excluded_examples) < 5: # 只記錄前 5 個範例
                        excluded_examples.append(gene_id)
                    continue

            # 初始化並加入
            if go_id not in go_features:
                go_features[go_id] = set()
            
            go_features[go_id].add(gene_id)

    # 顯示過濾結果
    if valid_genes is not None:
        print(f"[統計] 共排除 {excluded_count} 筆資料 (因為該基因不是 Live 或不在清單中)。")
        if excluded_examples:
            print(f"       範例例外 ID (可能是 Dead 或 Pseudogene): {', '.join(excluded_examples)} ...")

    print(f"[進度] 資料處理完畢，共找到 {len(go_features)} 個 GO 特徵，正在寫入 CSV...")

    # 2. 寫入結果 CSV
    # 格式：Feature_Name (GO ID), Count, Members (WormBase IDs)
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        
        # 寫入標題
        writer.writerow(['Feature_Name', 'Count', 'Members'])
        
        for go_id, genes in go_features.items():
            # 轉回 list 並排序
            sorted_genes = sorted(list(genes))
            count = len(sorted_genes)
            
            # 用逗號串接 ID
            members_str = ",".join(sorted_genes)
            
            # 寫入一行
            writer.writerow([go_id, count, members_str])

    print(f"[完成] 轉換成功！結果已儲存至：{output_path}")

# --- 執行設定 (請確認檔名正確) ---

# 1. 原始 GAF 檔案 (從 WormBase 下載的)
input_gaf_file = 'gene_association.WS298.wb.c_elegans'       

# 2. 母體清單檔案 (包含 Live/Dead 資訊的那個 txt)
background_list_file = 'c_elegans.PRJNA13758.WS298.geneOtherIDs.txt' 

# 3. 輸出結果
output_csv_file = 'GO_feature_table_WBID.csv' 

# 執行主程式
try:
    process_gaf_to_feature_table(input_gaf_file, background_list_file, output_csv_file)
except FileNotFoundError as e:
    print(f"錯誤：找不到檔案 - {e}")
except Exception as e:
    print(f"發生未預期的錯誤：{e}")