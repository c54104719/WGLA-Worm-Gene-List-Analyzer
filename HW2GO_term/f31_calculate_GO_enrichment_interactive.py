"""
交互式 GO 富集分析程序
用户输入基因列表和选择 C/F/P aspect，计算 GO 富集结果
"""

import csv
import math
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


# 常量定义
TOTAL_GENES = 49164  # 背景基因总数（Live genes）
ASPECTS = {
    'C': ('GO_Term_live_only_aspect_C_CellularComponent_live_only.csv', 'Cellular Component'),
    'F': ('GO_Term_live_only_aspect_F_MolecularFunction_live_only.csv', 'Molecular Function'),
    'P': ('GO_Term_live_only_aspect_P_BiologicalProcess_live_only.csv', 'Biological Process'),
}

# 结果表头
RESULT_HEADER = [
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


def load_input_genes(input_file):
    """
    从文件或直接输入加载基因列表
    """
    genes = set()
    try:
        with open(input_file, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    genes.add(gene)
    except FileNotFoundError:
        print(f"文件 {input_file} 不存在")
        return None
    
    print(f"载入了 {len(genes)} 个基因")
    return genes


def load_go_feature_table(feature_file):
    """
    加载 GO 特征表（Domains, number, ID）
    返回字典：{go_id: set(genes_in_background)}
    """
    go_genes = {}
    try:
        with open(feature_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                go_id = row['Domains'].strip()
                genes_str = row['ID'].strip().strip('"')  # 移除可能的引号
                genes = set(gene.strip() for gene in genes_str.split(','))
                go_genes[go_id] = genes
    except FileNotFoundError:
        print(f"文件 {feature_file} 不存在")
        return None
    
    print(f"从 {feature_file} 载入了 {len(go_genes)} 个 GO 项目")
    return go_genes


def calculate_enrichment(input_genes, go_genes):
    """
    计算 GO 富集分析
    
    参数:
    - input_genes: 用户输入的基因集合
    - go_genes: 背景中的 GO 项目到基因映射 {go_id: set(genes)}
    
    返回: 富集结果列表
    """
    results = []
    input_count = len(input_genes)
    
    for go_id, bg_genes in go_genes.items():
        # 计算 2x2 列联表
        # 背景中有该 GO 的基因数
        bg_with_go = len(bg_genes)
        # 背景中没有该 GO 的基因数
        bg_without_go = TOTAL_GENES - bg_with_go
        
        # 输入中有该 GO 的基因数
        input_with_go = len(input_genes & bg_genes)
        # 输入中没有该 GO 的基因数
        input_without_go = input_count - input_with_go
        
        # Fisher's Exact Test
        oddsratio, pvalue = fisher_exact(
            [[input_with_go, input_without_go],
             [bg_with_go, bg_without_go]]
        )
        
        # 计算比例和富集倍数
        expected_ratio = bg_with_go / TOTAL_GENES
        observed_ratio = input_with_go / input_count if input_count > 0 else 0
        fold_change = observed_ratio / expected_ratio if expected_ratio > 0 else 0
        
        # log2 变换
        log2_fc = math.log2(fold_change) if fold_change > 0 else -1000.0
        
        # -log10 变换
        log10_p = -math.log10(pvalue) if pvalue > 0 else 0
        
        results.append({
            'go_id': go_id,
            'expected_ratio': expected_ratio,
            'expected_ratio_str': f"{bg_with_go}/{TOTAL_GENES} ({expected_ratio*100:.4f}%)",
            'observed_ratio': observed_ratio,
            'observed_ratio_str': f"{input_with_go}/{input_count} ({observed_ratio*100:.4f}%)" if input_count > 0 else f"0/{input_count} (0.0%)",
            'fold_change': fold_change,
            'pvalue': pvalue,
            'log2_fc': log2_fc,
            'log10_p': log10_p,
        })
    
    return results


def apply_correction(results):
    """
    应用 FDR (Benjamini-Hochberg) 和 Bonferroni 多重检验校正
    """
    pvalues = [r['pvalue'] for r in results]
    
    # FDR 校正
    reject_fdr, pvals_fdr, _, _ = multipletests(pvalues, method='fdr_bh')
    
    # Bonferroni 校正
    reject_bonf, pvals_bonf, _, _ = multipletests(pvalues, method='bonferroni')
    
    # 将校正后的 p-value 加入结果
    for i, result in enumerate(results):
        result['pvalue_fdr'] = pvals_fdr[i]
        result['pvalue_bonf'] = pvals_bonf[i]
        result['log10_p_fdr'] = -math.log10(pvals_fdr[i]) if pvals_fdr[i] > 0 else 0
        result['log10_p_bonf'] = -math.log10(pvals_bonf[i]) if pvals_bonf[i] > 0 else 0
    
    return results


def save_results(results, output_file, aspect_name):
    """
    将结果保存到 CSV 文件
    """
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(RESULT_HEADER)
        
        for r in results:
            writer.writerow([
                r['go_id'],
                r['expected_ratio_str'],
                r['observed_ratio_str'],
                f"{r['fold_change']:.4f}",
                f"{r['pvalue']:.6f}",
                f"{r['pvalue_fdr']:.6f}",
                f"{r['pvalue_bonf']:.6f}",
                f"{r['log2_fc']:.4f}",
                f"{r['log10_p']:.4f}",
                f"{r['log10_p_fdr']:.4f}",
                f"{r['log10_p_bonf']:.4f}",
            ])
    
    print(f"结果已保存到 {output_file}")


def main():
    """
    主程序流程
    """
    print("=" * 60)
    print("GO 富集分析程序")
    print("=" * 60)
    
    # 1. 获取输入文件
    input_file = input("\n请输入基因列表文件路径 (文本文件，每行一个 WBGene ID): ").strip()
    input_genes = load_input_genes(input_file)
    
    if input_genes is None or len(input_genes) == 0:
        print("无法加载基因列表，程序退出")
        return
    
    # 2. 获取 aspect 选择
    print("\n可选择的 aspect:")
    print("  C - Cellular Component (细胞成分)")
    print("  F - Molecular Function (分子功能)")
    print("  P - Biological Process (生物过程)")
    
    aspect_input = input("\n请选择 aspect (多个请用逗号分隔，如: C,F,P): ").strip().upper()
    selected_aspects = [a.strip() for a in aspect_input.split(',') if a.strip() in ASPECTS]
    
    if not selected_aspects:
        print("未选择有效的 aspect，程序退出")
        return
    
    print(f"\n已选择: {', '.join([ASPECTS[a][1] for a in selected_aspects])}")
    
    # 3. 对每个 aspect 进行分析
    for aspect in selected_aspects:
        feature_file, aspect_name = ASPECTS[aspect]
        print(f"\n{'='*60}")
        print(f"分析 {aspect_name}...")
        print(f"{'='*60}")
        
        # 加载特征表
        go_genes = load_go_feature_table(feature_file)
        if go_genes is None:
            print(f"跳过 {aspect_name}")
            continue
        
        # 计算富集
        results = calculate_enrichment(input_genes, go_genes)
        print(f"计算了 {len(results)} 个 GO 项目的富集统计")
        
        # 应用多重检验校正
        results = apply_correction(results)
        
        # 保存结果
        output_file = f"GO_Term_enrichment_result_{aspect}.csv"
        save_results(results, output_file, aspect_name)
    
    print(f"\n{'='*60}")
    print("分析完成！")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
