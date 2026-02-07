"""
GO 項詞富集分析工具
輸入：基因列表 + GO 特徵表 (分別為 C/F/P 三種類型)
輸出：6 個結果 CSV 檔案 (Enriched 和 Depleted 各 3 個)

計算方法：Fisher's Exact Test + 多重檢定校正 (FDR Benjamini-Hochberg, Bonferroni)
"""

import csv
import math
import scipy.stats
from statsmodels.stats.multitest import multipletests


def load_input_genes(filepath):
    """讀取輸入的基因列表"""
    genes = set()
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            gene_id = line.strip()
            if gene_id:
                genes.add(gene_id)
    return genes


def load_feature_table(filepath):
    """
    讀取 GO 特徵表
    格式：Domains, number, ID
    回傳：dict {GO_ID: set(基因IDs)}
    """
    features = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            go_id = row['Domains']
            id_str = row['ID']
            
            # 解析基因列表 (可能是逗號分隔或引號包含)
            if id_str.startswith('"') and id_str.endswith('"'):
                id_str = id_str[1:-1]
            
            gene_list = set(id.strip() for id in id_str.split(',') if id.strip())
            features[go_id] = gene_list
    
    return features


def calculate_background_size(feature_table):
    """計算背景基因總數 (所有不重複的基因)"""
    all_genes = set()
    for genes in feature_table.values():
        all_genes.update(genes)
    return len(all_genes)


def perform_fisher_exact_test(input_genes, feature_table):
    """
    對每個 GO term 進行 Fisher's Exact Test
    
    回傳：list of dict，包含所有計算結果
    """
    results = []
    
    B = len(input_genes)  # 輸入基因數
    D = calculate_background_size(feature_table)  # 背景基因總數
    
    for go_id, go_genes in feature_table.items():
        C = len(go_genes)  # 該 GO term 的基因數
        A = len(input_genes.intersection(go_genes))  # 交集
        
        # 計算 ratios
        observed_ratio = A / B if B > 0 else 0
        expected_ratio = C / D if D > 0 else 0
        fold_change = observed_ratio / expected_ratio if expected_ratio > 0 else 0
        
        # 計算 log2(FC)
        log2_fc = math.log2(fold_change) if fold_change > 0 else -1000.0
        
        # Fisher's Exact Test
        # 2x2 表格: [[A, C-A], [B-A, D-C-B+A]]
        cont_table = [[A, C-A], [B-A, D-C-B+A]]
        
        # Enrichment (基因過度表示)
        _, p_enrich = scipy.stats.fisher_exact(cont_table, alternative='greater')
        
        # Depletion (基因低表示)
        _, p_deplete = scipy.stats.fisher_exact(cont_table, alternative='less')
        
        # 格式化 ratio 字串
        exp_str = f"{C}/{D} ({expected_ratio:.3%})"
        obs_str = f"{A}/{B} ({observed_ratio:.1%})"
        
        results.append({
            'go_id': go_id,
            'A': A,
            'B': B,
            'C': C,
            'D': D,
            'observed_ratio': observed_ratio,
            'expected_ratio': expected_ratio,
            'fold_change': fold_change,
            'log2_fc': log2_fc,
            'exp_str': exp_str,
            'obs_str': obs_str,
            'p_enrich': p_enrich,
            'p_deplete': p_deplete
        })
    
    return results, B, D


def apply_multiple_test_correction(results):
    """
    應用 FDR 和 Bonferroni 校正
    """
    # 提取 p-values
    p_enrich_list = [r['p_enrich'] for r in results]
    p_deplete_list = [r['p_deplete'] for r in results]
    
    # FDR 校正 (Benjamini-Hochberg)
    fdr_enrich = multipletests(p_enrich_list, method='fdr_bh')[1]
    fdr_deplete = multipletests(p_deplete_list, method='fdr_bh')[1]
    
    # Bonferroni 校正
    bonf_enrich = multipletests(p_enrich_list, method='bonferroni')[1]
    bonf_deplete = multipletests(p_deplete_list, method='bonferroni')[1]
    
    # 加入校正後的 p-values
    for i, result in enumerate(results):
        result['p_enrich_fdr'] = fdr_enrich[i]
        result['p_enrich_bonf'] = bonf_enrich[i]
        result['p_deplete_fdr'] = fdr_deplete[i]
        result['p_deplete_bonf'] = bonf_deplete[i]
    
    return results


def calculate_log10_p(p_value):
    """計算 -log10(p-value)"""
    if p_value <= 0:
        return -0.0
    return -math.log10(p_value)


def write_enrichment_result(results, output_path):
    """輸出富集分析結果 (過度表示)"""
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        
        # Header
        header = [
            'Domain Name',
            'Expected Ratio',
            'Observed Ratio',
            'Enrichment Fold Change',
            'P-value',
            'FDR enriched P-value',
            'Bonferroni enriched P-value',
            'log2(Enrichment Fold Change)',
            '-log10(P-value)',
            '-log10(FDR P-value)',
            '-log10(Bonferroni P-value)'
        ]
        writer.writerow(header)
        
        # 按 GO ID 排序
        for result in sorted(results, key=lambda x: x['go_id']):
            log10_p = calculate_log10_p(result['p_enrich'])
            log10_p_fdr = calculate_log10_p(result['p_enrich_fdr'])
            log10_p_bonf = calculate_log10_p(result['p_enrich_bonf'])
            
            row = [
                result['go_id'],
                result['exp_str'],
                result['obs_str'],
                result['fold_change'],
                result['p_enrich'],
                result['p_enrich_fdr'],
                result['p_enrich_bonf'],
                result['log2_fc'],
                log10_p,
                log10_p_fdr,
                log10_p_bonf
            ]
            writer.writerow(row)


def write_depletion_result(results, output_path):
    """輸出富集分析結果 (低表示)"""
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        
        # Header (注意：FDR 和 Bonferroni 改為 depleted)
        header = [
            'Domain Name',
            'Expected Ratio',
            'Observed Ratio',
            'Enrichment Fold Change',
            'P-value',
            'FDR depleted P-value',
            'Bonferroni depleted P-value',
            'log2(Enrichment Fold Change)',
            '-log10(P-value)',
            '-log10(FDR P-value)',
            '-log10(Bonferroni P-value)'
        ]
        writer.writerow(header)
        
        # 按 GO ID 排序
        for result in sorted(results, key=lambda x: x['go_id']):
            log10_p = calculate_log10_p(result['p_deplete'])
            log10_p_fdr = calculate_log10_p(result['p_deplete_fdr'])
            log10_p_bonf = calculate_log10_p(result['p_deplete_bonf'])
            
            row = [
                result['go_id'],
                result['exp_str'],
                result['obs_str'],
                result['fold_change'],
                result['p_deplete'],
                result['p_deplete_fdr'],
                result['p_deplete_bonf'],
                result['log2_fc'],
                log10_p,
                log10_p_fdr,
                log10_p_bonf
            ]
            writer.writerow(row)


def run_enrichment_analysis(input_file, feature_table_file, output_prefix):
    """
    執行完整的富集分析流程
    
    Args:
        input_file: 輸入基因列表
        feature_table_file: GO 特徵表
        output_prefix: 輸出檔案的前綴 (如 'GO_Term_C')
    """
    print(f"[開始] 正在執行 GO 富集分析...")
    print(f"  輸入檔案: {input_file}")
    print(f"  特徵表: {feature_table_file}")
    
    # 1. 載入數據
    input_genes = load_input_genes(input_file)
    feature_table = load_feature_table(feature_table_file)
    
    print(f"  輸入基因數: {len(input_genes)}")
    print(f"  GO terms 數: {len(feature_table)}")
    
    # 2. 執行 Fisher's Exact Test
    results, B, D = perform_fisher_exact_test(input_genes, feature_table)
    print(f"  背景基因總數: {D}")
    
    # 3. 多重檢定校正
    results = apply_multiple_test_correction(results)
    
    # 4. 輸出結果
    enriched_output = f"{output_prefix}_enriched_result.csv"
    depleted_output = f"{output_prefix}_depleted_result.csv"
    
    write_enrichment_result(results, enriched_output)
    write_depletion_result(results, depleted_output)
    
    print(f"✓ 富集結果: {enriched_output}")
    print(f"✓ 低表示結果: {depleted_output}")
    print(f"[完成] GO 富集分析完畢！")


# --- 執行設定 ---
if __name__ == "__main__":
    # 輸入基因列表
    input_gene_file = "go_term_test_input.txt"
    
    # 三種 GO 類型分析
    go_types = {
        'C': 'GO_Term_domain_type_C_live_only.csv',
        'F': 'GO_Term_domain_type_F_live_only.csv',
        'P': 'GO_Term_domain_type_P_live_only.csv'
    }
    
    # 逐個執行三種類型の分析
    for go_type, feature_file in go_types.items():
        try:
            run_enrichment_analysis(
                input_gene_file,
                feature_file,
                f"GO_Term_{go_type}"
            )
        except FileNotFoundError as e:
            print(f"❌ 錯誤：找不到檔案 - {e}")
        except Exception as e:
            print(f"❌ 發生錯誤：{e}")
        
        print()  # 空行分隔
