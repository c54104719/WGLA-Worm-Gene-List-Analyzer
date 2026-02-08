"""
GO 富集分析程序（命令行版本）
使用方法: python calculate_GO_enrichment.py <input_gene_file> <aspect> [aspect2] [aspect3]
示例: python calculate_GO_enrichment.py genes.txt C F P
"""

import csv
import math
import sys
import os
from pathlib import Path
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# Fix encoding for Windows
if sys.stdout.encoding != 'utf-8':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# Constants
TOTAL_GENES = 49164
ASPECTS = {
    'C': (r'C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\enrich_worm_project\GO_Term_domain_type_C_live_only.csv', 'Cellular Component'),
    'F': (r'C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\enrich_worm_project\GO_Term_domain_type_F_live_only.csv', 'Molecular Function'),
    'P': (r'C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\enrich_worm_project\GO_Term_domain_type_P_live_only.csv', 'Biological Process'),
}

# Result header
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
    """Load gene list from file."""
    genes = set()
    input_path = Path(input_file)
    
    if not input_path.exists():
        print(f"Error: File {input_file} not found", file=sys.stderr)
        return None
    
    try:
        with open(input_path, 'r', encoding='utf-8') as f:
            for line in f:
                gene = line.strip()
                if gene and not gene.startswith('#'):
                    genes.add(gene)
    except Exception as e:
        print(f"Error: Cannot read file {input_file}: {e}", file=sys.stderr)
        return None
    
    print(f"[OK] Loaded {len(genes)} genes")
    return genes


def load_go_feature_table(feature_file):
    """Load GO feature table (Domains, number, ID)."""
    go_genes = {}
    feature_path = Path(feature_file)
    
    if not feature_path.exists():
        print(f"Error: File {feature_file} not found", file=sys.stderr)
        return None
    
    try:
        with open(feature_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                go_id = row['Domains'].strip()
                genes_str = row['ID'].strip().strip('"')
                genes = set(gene.strip() for gene in genes_str.split(',') if gene.strip())
                go_genes[go_id] = genes
    except Exception as e:
        print(f"Error: Cannot read file {feature_file}: {e}", file=sys.stderr)
        return None
    
    print(f"[OK] Loaded {len(go_genes)} GO terms from {Path(feature_file).name}")
    return go_genes


def calculate_enrichment(input_genes, go_genes):
    """Calculate GO enrichment and depletion analysis."""
    results = []
    input_count = len(input_genes)
    
    for go_id, bg_genes in go_genes.items():
        # 2x2 contingency table construction
        input_count = len(input_genes)
        
        # A: genes in input AND have GO term
        input_with_go = len(input_genes & bg_genes)
        # B - A: genes in input but WITHOUT GO term
        input_without_go = input_count - input_with_go
        
        # C - A: genes in background WITHOUT input that have GO term
        bg_with_go = len(bg_genes) - input_with_go
        # D - B - (C - A): genes in background WITHOUT input and WITHOUT GO term
        bg_without_go = TOTAL_GENES - input_count - bg_with_go
        
        # Fisher's Exact Test (both enrichment and depletion)
        # Enrichment (one-tailed, greater)
        _, pvalue_enrich = fisher_exact(
            [[input_with_go, input_without_go],
             [bg_with_go, bg_without_go]],
            alternative='greater'
        )
        
        # Depletion (one-tailed, less)
        _, pvalue_deplete = fisher_exact(
            [[input_with_go, input_without_go],
             [bg_with_go, bg_without_go]],
            alternative='less'
        )
        
        # Calculate ratios and fold change
        # Expected ratio based on total background with GO term
        feature_count = len(bg_genes)
        expected_ratio = feature_count / TOTAL_GENES
        observed_ratio = input_with_go / input_count if input_count > 0 else 0
        fold_change = observed_ratio / expected_ratio if expected_ratio > 0 else 0
        
        # log2 transformation
        log2_fc = math.log2(fold_change) if fold_change > 0 else -1000.0
        
        # -log10 transformation
        log10_p_enrich = -math.log10(pvalue_enrich) if pvalue_enrich > 0 else 0
        log10_p_deplete = -math.log10(pvalue_deplete) if pvalue_deplete > 0 else 0
        
        results.append({
            'go_id': go_id,
            'expected_ratio': expected_ratio,
            'expected_ratio_str': f"{feature_count}/{TOTAL_GENES} ({expected_ratio*100:.4f}%)",
            'observed_ratio': observed_ratio,
            'observed_ratio_str': f"{input_with_go}/{input_count} ({observed_ratio*100:.4f}%)" if input_count > 0 else f"0/{input_count} (0.0%)",
            'fold_change': fold_change,
            'pvalue_enrich': pvalue_enrich,
            'pvalue_deplete': pvalue_deplete,
            'log2_fc': log2_fc,
            'log10_p_enrich': log10_p_enrich,
            'log10_p_deplete': log10_p_deplete,
        })
    
    return results
    
    return results


def apply_correction(results):
    """Apply FDR and Bonferroni multiple test correction for both enrichment and depletion."""
    pvalues_enrich = [r['pvalue_enrich'] for r in results]
    pvalues_deplete = [r['pvalue_deplete'] for r in results]
    
    # FDR correction
    reject_fdr_enrich, pvals_fdr_enrich, _, _ = multipletests(pvalues_enrich, method='fdr_bh')
    reject_fdr_deplete, pvals_fdr_deplete, _, _ = multipletests(pvalues_deplete, method='fdr_bh')
    
    # Bonferroni correction
    reject_bonf_enrich, pvals_bonf_enrich, _, _ = multipletests(pvalues_enrich, method='bonferroni')
    reject_bonf_deplete, pvals_bonf_deplete, _, _ = multipletests(pvalues_deplete, method='bonferroni')
    
    # Add corrected p-values to results
    for i, result in enumerate(results):
        # Enrichment corrections
        result['pvalue_enrich_fdr'] = pvals_fdr_enrich[i]
        result['pvalue_enrich_bonf'] = pvals_bonf_enrich[i]
        result['log10_p_enrich_fdr'] = -math.log10(pvals_fdr_enrich[i]) if pvals_fdr_enrich[i] > 0 else 0
        result['log10_p_enrich_bonf'] = -math.log10(pvals_bonf_enrich[i]) if pvals_bonf_enrich[i] > 0 else 0
        
        # Depletion corrections
        result['pvalue_deplete_fdr'] = pvals_fdr_deplete[i]
        result['pvalue_deplete_bonf'] = pvals_bonf_deplete[i]
        result['log10_p_deplete_fdr'] = -math.log10(pvals_fdr_deplete[i]) if pvals_fdr_deplete[i] > 0 else 0
        result['log10_p_deplete_bonf'] = -math.log10(pvals_bonf_deplete[i]) if pvals_bonf_deplete[i] > 0 else 0
    
    return results


def save_results(results, output_file, aspect_name, result_type='enrichment'):
    """Save results to CSV file.
    
    Args:
        results: Analysis results
        output_file: Output file path
        aspect_name: GO aspect name (for logging)
        result_type: 'enrichment' or 'depletion'
    """
    try:
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(RESULT_HEADER)
            
            for r in results:
                if result_type == 'enrichment':
                    pvalue = r['pvalue_enrich']
                    pvalue_fdr = r['pvalue_enrich_fdr']
                    pvalue_bonf = r['pvalue_enrich_bonf']
                    log10_p = r['log10_p_enrich']
                    log10_p_fdr = r['log10_p_enrich_fdr']
                    log10_p_bonf = r['log10_p_enrich_bonf']
                else:  # depletion
                    pvalue = r['pvalue_deplete']
                    pvalue_fdr = r['pvalue_deplete_fdr']
                    pvalue_bonf = r['pvalue_deplete_bonf']
                    log10_p = r['log10_p_deplete']
                    log10_p_fdr = r['log10_p_deplete_fdr']
                    log10_p_bonf = r['log10_p_deplete_bonf']
                
                writer.writerow([
                    r['go_id'],
                    r['expected_ratio_str'],
                    r['observed_ratio_str'],
                    f"{r['fold_change']:.14f}",
                    f"{pvalue:.16e}",
                    f"{pvalue_fdr:.16e}",
                    f"{pvalue_bonf:.16e}",
                    f"{r['log2_fc']:.16e}",
                    f"{log10_p:.16e}",
                    f"{log10_p_fdr:.16e}",
                    f"{log10_p_bonf:.16e}",
                ])
        
        print(f"[OK] Results saved: {output_file} ({len(results)} rows)")
        return True
    except Exception as e:
        print(f"Error: Cannot save results file: {e}", file=sys.stderr)
        return False


def main():
    """Main program flow."""
    if len(sys.argv) < 3:
        print("Usage:")
        print("  python calculate_GO_enrichment.py <input_gene_file> <aspect> [aspect2] [aspect3]")
        print("\nParameters:")
        print("  input_gene_file: Gene list file (text file, one WBGene ID per line)")
        print("  aspect:          Aspect type (C, F, P)")
        print("\nExample:")
        print("  python calculate_GO_enrichment.py genes.txt C")
        print("  python calculate_GO_enrichment.py genes.txt C F P")
        sys.exit(1)
    
    input_file = sys.argv[1]
    selected_aspects = [a.upper() for a in sys.argv[2:] if a.upper() in ASPECTS]
    
    if not selected_aspects:
        print("Error: No valid aspects selected (C, F, P)", file=sys.stderr)
        sys.exit(1)
    
    # Load gene list
    print("=" * 60)
    print("GO Enrichment Analysis Program")
    print("=" * 60)
    print()
    
    input_genes = load_input_genes(input_file)
    if input_genes is None or len(input_genes) == 0:
        print("Error: Cannot load gene list", file=sys.stderr)
        sys.exit(1)
    
    print(f"[OK] Selected aspects: {', '.join([ASPECTS[a][1] for a in selected_aspects])}")
    print()
    
    # Analyze each aspect
    success_count = 0
    for aspect in selected_aspects:
        feature_file, aspect_name = ASPECTS[aspect]
        print(f"Processing {aspect_name}...")
        
        # Load feature table
        go_genes = load_go_feature_table(feature_file)
        if go_genes is None:
            print(f"[WARN] Cannot load {feature_file}, skipping")
            continue
        
        # Calculate enrichment
        results = calculate_enrichment(input_genes, go_genes)
        print(f"[OK] Calculated enrichment for {len(results)} GO terms")
        
        # Apply multiple test correction
        results = apply_correction(results)
        print(f"[OK] Applied multiple test correction (FDR, Bonferroni)")
        
        # Save enrichment results
        output_file_enrich = f"GO_Term_enrichment_result_{aspect}.csv"
        if save_results(results, output_file_enrich, aspect_name, result_type='enrichment'):
            success_count += 1
        
        # Save depletion results
        output_file_deplete = f"GO_Term_depletion_result_{aspect}.csv"
        if save_results(results, output_file_deplete, aspect_name, result_type='depletion'):
            success_count += 1
        
        print()
    
    print("=" * 60)
    expected_files = len(selected_aspects) * 2  # 2 files per aspect (enrichment + depletion)
    if success_count == expected_files:
        print(f"[OK] Analysis completed! Generated {success_count} result files")
        print("=" * 60)
        sys.exit(0)
    else:
        print(f"[WARN] Analysis partially completed ({success_count}/{expected_files} succeeded)")
        print("=" * 60)
        sys.exit(1)


if __name__ == '__main__':
    main()
