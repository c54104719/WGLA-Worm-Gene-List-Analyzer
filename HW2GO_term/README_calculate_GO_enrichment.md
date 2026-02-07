# GO Enrichment Analysis Program Usage Guide

## Overview
`calculate_GO_enrichment.py` is a command-line tool that performs GO (Gene Ontology) enrichment analysis based on user-provided gene lists. It supports analysis for three GO aspects: Cellular Component (C), Molecular Function (F), and Biological Process (P).

## Features
- **Multi-aspect analysis**: Analyze any combination of C, F, and P aspects
- **Fisher's Exact Test**: Performs Fisher's exact test for statistical significance
- **Multiple test correction**: Applies FDR (Benjamini-Hochberg) and Bonferroni corrections
- **Comprehensive output**: Generates detailed results with fold changes and p-values

## Installation Requirements
```bash
pip install scipy pandas statsmodels
```

## Input Files Required
The program requires three feature table files in the same directory:
- `GO_Term_live_only_aspect_C_CellularComponent_live_only.csv` (Cellular Component)
- `GO_Term_live_only_aspect_F_MolecularFunction_live_only.csv` (Molecular Function)  
- `GO_Term_live_only_aspect_P_BiologicalProcess_live_only.csv` (Biological Process)

Each feature table should have columns:
- `Domains`: GO term ID (e.g., GO:0000015)
- `number`: Number of genes with this term
- `ID`: Comma-separated list of WormBase gene IDs

## Input Gene List Format
Create a text file with one WormBase gene ID per line:
```
WBGene00002324
WBGene00003164
WBGene00004343
...
```

## Usage

### Basic Syntax
```bash
python calculate_GO_enrichment.py <input_gene_file> <ASPECT> [ASPECT2] [ASPECT3]
```

### Parameters
- `<input_gene_file>`: Path to text file containing WormBase gene IDs
- `<ASPECT>`: One or more aspect codes: C (Cellular Component), F (Molecular Function), P (Biological Process)

### Examples

**Analyze single aspect:**
```bash
python calculate_GO_enrichment.py genes.txt C
```

**Analyze multiple aspects:**
```bash
python calculate_GO_enrichment.py genes.txt C F P
```

**Analyze F and P only:**
```bash
python calculate_GO_enrichment.py genes.txt F P
```

## Output Files
The program generates CSV result files for each selected aspect:
- `GO_Term_enrichment_result_C.csv` - Cellular Component results
- `GO_Term_enrichment_result_F.csv` - Molecular Function results
- `GO_Term_enrichment_result_P.csv` - Biological Process results

## Output Format
Each result file contains the following columns:

| Column | Description |
|--------|-------------|
| Domain Name | GO term ID |
| Expected Ratio | Background frequency of the GO term |
| Observed Ratio | Frequency in the input gene list |
| Enrichment Fold Change | Observed/Expected ratio |
| P-value | Fisher's exact test p-value |
| FDR enriched P-value | Benjamini-Hochberg corrected p-value |
| Bonferroni enriched P-value | Bonferroni corrected p-value |
| log2(Enrichment Fold Change) | Log2 transformation of fold change |
| -log10(P-value) | -log10 transformation of p-value |
| -log10(FDR P-value) | -log10 transformation of FDR p-value |
| -log10(Bonferroni P-value) | -log10 transformation of Bonferroni p-value |

## Example Output (First 3 rows of C aspect results)
```
Domain Name,Expected Ratio,Observed Ratio,Enrichment Fold Change,P-value,FDR enriched P-value,Bonferroni enriched P-value,log2(Enrichment Fold Change),-log10(P-value),-log10(FDR P-value),-log10(Bonferroni P-value)
GO:0000015,1/49164 (0.0020%),0/32 (0.0000%),0.0000,1.000000,1.000000,1.000000,-1000.0000,-0.0000,-0.0000,-0.0000
GO:0000109,1/49164 (0.0020%),0/32 (0.0000%),0.0000,1.000000,1.000000,1.000000,-1000.0000,-0.0000,-0.0000,-0.0000
GO:0000110,3/49164 (0.0061%),0/32 (0.0000%),0.0000,1.000000,1.000000,1.000000,-1000.0000,-0.0000,-0.0000,-0.0000
```

## Statistical Methods

### Fisher's Exact Test
Uses a 2x2 contingency table to test for significant association between GO terms and the input gene list:
- Cell A: Input genes with the GO term
- Cell B: Input genes without the GO term
- Cell C: Background genes with the GO term
- Cell D: Background genes without the GO term

### FDR Correction (Benjamini-Hochberg)
Controls the false discovery rate across all GO terms tested. Use this for exploratory analysis.

### Bonferroni Correction
More conservative correction controlling family-wise error rate. Use for confirmatory analysis.

## Background Population
- **Background genes**: 49,164 live C. elegans genes
- **Background source**: WormBase Gene Ontology Annotation File (GAF)

## Notes
- The program only considers "Live" genes from the WormBase annotation
- All three feature table files must be present in the working directory
- Input gene IDs must be WormBase gene IDs (WBGene format)
- The program handles large gene lists efficiently
