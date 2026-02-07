"""
GO 項詞富集分析工具函式集
可供 tool2_view 調用
"""

import math
import scipy.stats
from statsmodels.stats.multitest import multipletests


def perform_go_fisher_exact_test(input_genes, feature_table):
    """
    對 GO 特徵表進行 Fisher's Exact Test
    
    Args:
        input_genes: set of gene IDs (輸入基因)
        feature_table: dict {GO_ID: set(基因IDs)} (GO 特徵表)
    
    Returns:
        results: list of dict，包含每個 GO term 的計算結果
        B: 輸入基因數
        D: 背景基因總數
    """
    results = []
    
    B = len(input_genes)  # 輸入基因數
    
    # 計算背景基因總數
    all_genes = set()
    for genes in feature_table.values():
        all_genes.update(genes)
    D = len(all_genes)  # 背景基因總數
    
    # 逐個 GO term 進行計算
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
    應用 FDR (Benjamini-Hochberg) 和 Bonferroni 校正
    
    Args:
        results: 原始計算結果
    
    Returns:
        results: 加入校正後 p-values 的結果
    """
    # 提取 p-values
    p_enrich_list = [r['p_enrich'] for r in results]
    p_deplete_list = [r['p_deplete'] for r in results]
    
    # FDR 校正
    fdr_enrich = multipletests(p_enrich_list, method='fdr_bh')[1]
    fdr_deplete = multipletests(p_deplete_list, method='fdr_bh')[1]
    
    # Bonferroni 校正
    bonf_enrich = multipletests(p_enrich_list, method='bonferroni')[1]
    bonf_deplete = multipletests(p_deplete_list, method='bonferroni')[1]
    
    # 加入校正結果
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
