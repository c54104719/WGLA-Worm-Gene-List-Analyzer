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
    'C': ('GO_Term_live_only_aspect_C_CellularComponent_live_only.csv', 'Cellular Component'),
    'F': ('GO_Term_live_only_aspect_F_MolecularFunction_live_only.csv', 'Molecular Function'),
    'P': ('GO_Term_live_only_aspect_P_BiologicalProcess_live_only.csv', 'Biological Process'),
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
    """Calculate GO enrichment analysis."""
    results = []
    input_count = len(input_genes)
    
    for go_id, bg_genes in go_genes.items():
        # 2x2 contingency table
        bg_with_go = len(bg_genes)
        bg_without_go = TOTAL_GENES - bg_with_go
        
        input_with_go = len(input_genes & bg_genes)
        input_without_go = input_count - input_with_go
        
        # Fisher's Exact Test
        oddsratio, pvalue = fisher_exact(
            [[input_with_go, input_without_go],
             [bg_with_go, bg_without_go]]
        )
        
        # Calculate ratios and fold change
        expected_ratio = bg_with_go / TOTAL_GENES
        observed_ratio = input_with_go / input_count if input_count > 0 else 0
        fold_change = observed_ratio / expected_ratio if expected_ratio > 0 else 0
        
        # log2 transformation
        log2_fc = math.log2(fold_change) if fold_change > 0 else -1000.0
        
        # -log10 transformation
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
    """Apply FDR and Bonferroni multiple test correction."""
    pvalues = [r['pvalue'] for r in results]
    
    # FDR correction
    reject_fdr, pvals_fdr, _, _ = multipletests(pvalues, method='fdr_bh')
    
    # Bonferroni correction
    reject_bonf, pvals_bonf, _, _ = multipletests(pvalues, method='bonferroni')
    
    for i, result in enumerate(results):
        result['pvalue_fdr'] = pvals_fdr[i]
        result['pvalue_bonf'] = pvals_bonf[i]
        result['log10_p_fdr'] = -math.log10(pvals_fdr[i]) if pvals_fdr[i] > 0 else 0
        result['log10_p_bonf'] = -math.log10(pvals_bonf[i]) if pvals_bonf[i] > 0 else 0
    
    return results


def save_results(results, output_file, aspect_name):
    """Save results to CSV file."""
    try:
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
        
        # Save results
        output_file = f"GO_Term_enrichment_result_{aspect}.csv"
        if save_results(results, output_file, aspect_name):
            success_count += 1
        
        print()
    
    print("=" * 60)
    if success_count == len(selected_aspects):
        print(f"[OK] Analysis completed! Generated {success_count} result files")
        print("=" * 60)
        sys.exit(0)
    else:
        print(f"[WARN] Analysis partially completed ({success_count}/{len(selected_aspects)} succeeded)")
        print("=" * 60)
        sys.exit(1)


if __name__ == '__main__':
    main()
