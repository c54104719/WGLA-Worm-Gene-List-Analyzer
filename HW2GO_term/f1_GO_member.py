import csv

def process_gaf_to_feature_table(input_path, output_path):
    # 1. 建立字典來儲存資料：Key 是 GO ID，Value 是基因集合 (用 set 自動去重)
    # 影片提到：要將每個 GO 對應的基因統計出來，並做 unique 處理 [[1]](https://app.heptabase.com/84c1ad40-c9a3-45ea-b0af-5e83d6c372b2/card/be56416f-12c9-420a-8786-0937cc17305b?playbackStart=1705.385&playbackEnd=1768.8429999999998)
    go_features = {}

    print("正在讀取原始檔案...")
    
    with open(input_path, 'r', encoding='utf-8') as f:
        for line in f:
            # 跳過註解行 (以 ! 開頭)
            if line.startswith('!'):
                continue
            
            # GAF 檔案是用 Tab 分隔的
            columns = line.strip().split('\t')
            
            # 確保欄位足夠多，避免讀到空行報錯
            if len(columns) < 5:
                continue

            # 根據 GAF 2.1 格式標準 & 影片內容：
            # 第 3 欄 (Index 2) 是 Gene Symbol (如 aap-1) [[3]](https://app.heptabase.com/84c1ad40-c9a3-45ea-b0af-5e83d6c372b2/card/be56416f-12c9-420a-8786-0937cc17305b?playbackStart=1632.068&playbackEnd=1705.375)
            # 第 5 欄 (Index 4) 是 GO ID (如 GO:0005942)
            gene_symbol = columns[2]
            go_id = columns[4]
            
            # 初始化該 GO ID (如果尚未存在)
            if go_id not in go_features:
                go_features[go_id] = set()
            
            # 加入基因 (set 會自動處理重複，不用寫 if 判斷)
            go_features[go_id].add(gene_symbol)

    print(f"資料讀取完畢，共找到 {len(go_features)} 個 GO 特徵，正在寫入 CSV...")

    # 2. 寫入結果 CSV
    # 影片要求的格式：GO ID, 成員數, 成員名單(逗號分隔) [[4]](https://app.heptabase.com/84c1ad40-c9a3-45ea-b0af-5e83d6c372b2/card/be56416f-12c9-420a-8786-0937cc17305b?playbackStart=1.625&playbackEnd=66.481)[[2]](https://app.heptabase.com/84c1ad40-c9a3-45ea-b0af-5e83d6c372b2/card/be56416f-12c9-420a-8786-0937cc17305b?playbackStart=1768.843&playbackEnd=1841.002)
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        
        # 寫入標題 (Header)
        writer.writerow(['Feature_ID', 'Count', 'Members'])
        
        for go_id, genes in go_features.items():
            # 將 set 轉回 list 並排序，確保輸出順序固定
            sorted_genes = sorted(list(genes))
            
            # 計算數量
            count = len(sorted_genes)
            
            # 將基因名單串接成字串，例如 "aap-1, ag-1"
            members_str = ", ".join(sorted_genes)
            
            # 寫入一行
            writer.writerow([go_id, count, members_str])

    print(f"轉換完成！結果已儲存至：{output_path}")

# --- 執行設定 ---
# 請將這裡的檔名改成你實際下載的檔名
input_filename = 'gene_association.WS298.wb.c_elegans'  # WormBase 下載的原始檔
output_filename = 'GO_feature_table.csv' # 輸出的結果檔

# 執行函式
try:
    process_gaf_to_feature_table(input_filename, output_filename)
except FileNotFoundError:
    print(f"錯誤：找不到檔案 '{input_filename}'，請確認檔案是否在同一資料夾內。")