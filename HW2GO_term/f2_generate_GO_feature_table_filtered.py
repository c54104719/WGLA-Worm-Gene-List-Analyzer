import csv
import os


def load_background_genes(filepath):
    """
    讀取母體清單 (c_elegans.PRJNA13758.WS298.geneOtherIDs.txt)。
    邏輯修改：
    1. 使用 Tab (\t) 切割欄位。
    2. 檢查第二欄 (Index 1) 是否為 'Live'。
    3. 只保留 'Live' 的基因，其他狀態排除
    """
    live_genes = set()
    live_count = 0
    non_live_types = {}

    if not os.path.exists(filepath):
        print(f"[警告] 找不到母體清單檔案 '{filepath}'，無法進行過濾。")
        return None

    print(f"[進度] 正在讀取母體清單: {filepath} ...")

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            columns = line.split('\t')
            if len(columns) < 2:
                continue

            gene_id = columns[0]
            status = columns[1]

            if status == 'Live':
                live_genes.add(gene_id)
                live_count += 1
            else:
                non_live_types[status] = non_live_types.get(status, 0) + 1

    print("[資訊] 成功讀取母體清單。")
    print(f"       - Live 基因: {live_count} 個 (保留)")
    if non_live_types:
        for status_type, count in non_live_types.items():
            print(f"       - 排除狀態 {status_type}: {count} 個")

    return live_genes


def process_gaf_by_aspect(gaf_path, background_path, output_prefix):
    """
    讀取 GAF 並按 Aspect 分流成三個表格 (P/F/C)

    Args:
        gaf_path: GAF 檔案路徑
        background_path: Live 清單檔案路徑
        output_prefix: 輸出檔案的前綴 (如 'GO_Term_non_live')
    """
    live_genes = load_background_genes(background_path)
    if live_genes is None:
        print("❌ 無法讀取背景清單，程式中止。")
        return

    aspect_features = {
        'P': {},  # Biological Process
        'F': {},  # Molecular Function
        'C': {}   # Cellular Component
    }

    total_lines = 0
    processed_lines = 0
    excluded_count = 0
    aspect_counts = {'P': 0, 'F': 0, 'C': 0, 'other': 0}

    print(f"\n[進度] 正在讀取 GAF 原始檔案: {gaf_path} ...")

    with open(gaf_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('!'):
                continue

            total_lines += 1
            columns = line.strip().split('\t')
            if len(columns) < 9:
                continue

            gene_id = columns[1]
            go_id = columns[4]
            aspect = columns[8]  # P, F, 或 C

            if gene_id not in live_genes:
                excluded_count += 1
                continue

            processed_lines += 1

            if aspect not in aspect_features:
                aspect_counts['other'] += 1
                continue

            aspect_counts[aspect] += 1

            if go_id not in aspect_features[aspect]:
                aspect_features[aspect][go_id] = set()

            aspect_features[aspect][go_id].add(gene_id)

    print("\n[統計] GAF 檔案處理結果：")
    print(f"       - 總行數: {total_lines}")
    print(f"       - 成功處理: {processed_lines} 筆")
    print(f"       - 排除 (Live 基因): {excluded_count} 筆")
    print("       - 按 Aspect 分佈：")
    print(f"         • P (Biological Process): {aspect_counts['P']} 筆")
    print(f"         • F (Molecular Function): {aspect_counts['F']} 筆")
    print(f"         • C (Cellular Component): {aspect_counts['C']} 筆")
    if aspect_counts['other'] > 0:
        print(f"         • 其他: {aspect_counts['other']} 筆")

    aspect_names = {
        'P': 'BiologicalProcess',
        'F': 'MolecularFunction',
        'C': 'CellularComponent'
    }

    for aspect_char, aspect_name in aspect_names.items():
        output_path = f"{output_prefix}_aspect_{aspect_char}_{aspect_name}_live_only.csv"
        write_feature_table(aspect_features[aspect_char], output_path)

    print("\n[完成] 所有分流表格已生成完畢！")


def write_feature_table(go_features, output_path):
    """
    寫入 CSV 格式的 GO 特徵表

    格式：Domains, number, ID
    """
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Domains', 'number', 'ID'])

        for go_id in sorted(go_features.keys()):
            genes = go_features[go_id]
            sorted_genes = sorted(list(genes))
            count = len(sorted_genes)
            members_str = ",".join(sorted_genes)
            writer.writerow([go_id, count, members_str])

    print(f"   ✓ 輸出檔案: {output_path} ({len(go_features)} 個 GO terms)")


if __name__ == "__main__":
    input_gaf_file = 'gene_association.WS298.wb.c_elegans'
    background_list_file = 'c_elegans.PRJNA13758.WS298.geneOtherIDs.txt'
    output_prefix = 'GO_Term_live_only'

    try:
        process_gaf_by_aspect(input_gaf_file, background_list_file, output_prefix)
    except FileNotFoundError as e:
        print(f"❌ 錯誤：找不到檔案 - {e}")
    except Exception as e:
        print(f"❌ 發生未預期的錯誤：{e}")
