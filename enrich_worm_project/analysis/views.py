from django.shortcuts import render
import json
import pandas as pd
import scipy.stats
from statsmodels.stats.multitest import multipletests
import math
import csv

# 導入 GO 富集分析工具函式
from .go_enrichment_utils import perform_go_fisher_exact_test, apply_multiple_test_correction, calculate_log10_p

# 設定全域變數 (酵母菌基因總數)
TOTAL_GENES = 6705

def home_view(request):
    """首頁視圖"""
    return render(request, 'home.html')


def load_go_feature_table(filepath):
    """
    讀取 GO 特徵表 CSV 檔案
    格式：Domains, number, ID
    回傳：dict {GO_ID: set(基因IDs)}
    """
    features = {}
    try:
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
    except FileNotFoundError:
        return None
    
    return features

def tool_view(request):
    """分析工具頁面"""
    context = {}
    
    if request.method == 'POST':
        # 1. 獲取使用者輸入的基因列表
        input_text = request.POST.get('gene_list', '')
        # 處理輸入：支援換行、逗號、空白分隔
        input_genes = set([x.strip() for x in input_text.replace('\n', ' ').replace(',', ' ').split() if x.strip()])
        
        if not input_genes:
            context['error'] = "請輸入至少一個基因 ID"
            return render(request, 'tool.html', context)

        # 2. 讀取特徵表
        try:
            domain_data = pd.read_excel("domain_test.xlsx")
        except Exception as e:
            context['error'] = f"找不到 domain_test.xlsx 或讀取錯誤: {str(e)}"
            return render(request, 'tool.html', context)

        results = []
        B = len(input_genes) # 輸入總數
        D = TOTAL_GENES      # 背景總數

        # 3. 遍歷所有 Domain 進行計算
        import math
        for index, row in domain_data.iterrows():
            domain_id = row["Domains"]
            # 假設 ID 欄位是以逗號分隔的基因字串
            domain_genes = set(str(row["ID"]).split(",")) 
            
            A = len(input_genes.intersection(domain_genes)) # 交集
            C = len(domain_genes) # 該 Domain 的基因數
            
            # Fisher's Exact Test
            # Enrichment (Greater)
            _, p_enrich = scipy.stats.fisher_exact([[A, C-A], [B-A, D-C-B+A]], alternative='greater')
            # Depletion (Less)
            _, p_deplete = scipy.stats.fisher_exact([[A, C-A], [B-A, D-C-B+A]], alternative='less')
            
            # 計算顯示用的數值
            observed_ratio = A / B if B > 0 else 0
            expected_ratio = C / D if D > 0 else 0
            fold_change = observed_ratio / expected_ratio if expected_ratio > 0 else 1
            # 轉換為 log2 fold change
            log2_fold_change = math.log2(fold_change) if fold_change > 0 else 0
            
            results.append({
                "domain": domain_id,
                "exp_str": f"{C}/{D} ({expected_ratio:.2%})",
                "obs_str": f"{A}/{B} ({observed_ratio:.2%})",
                "log2_fold_change": round(log2_fold_change, 2),
                "p_enrich_raw": p_enrich,
                "p_deplete_raw": p_deplete
            })

        # 4. 進行多重檢定校正 (FDR & Bonferroni)
        df = pd.DataFrame(results)
        
        # Enrichment 校正
        df['p_enrich_fdr'] = multipletests(df['p_enrich_raw'], method='fdr_bh')[1]
        df['p_enrich_bon'] = multipletests(df['p_enrich_raw'], method='bonferroni')[1]
        
        # Depletion 校正
        df['p_deplete_fdr'] = multipletests(df['p_deplete_raw'], method='fdr_bh')[1]
        df['p_deplete_bon'] = multipletests(df['p_deplete_raw'], method='bonferroni')[1]

        # 5. 轉換成 JSON 格式傳給前端
        full_data = df.to_dict(orient='records')
        context['result_json'] = json.dumps(full_data)
        context['input_count'] = B

    return render(request, 'tool.html', context)

def help_view(request):
    """幫助文檔頁面"""
    return render(request, 'help.html')

def developer_view(request):
    """開發者信息頁面"""
    return render(request, 'developer.html')

def tool2_view(request):
    """分析工具頁面 2 - 帶功能選擇 (包括 GO 項詞分析)"""
    context = {}
    
    if request.method == 'POST':
        # 1. 獲取使用者輸入的基因列表
        input_text = request.POST.get('gene_list', '')
        # 處理輸入：支援換行、逗號、空白分隔
        input_genes = set([x.strip() for x in input_text.replace('\n', ' ').replace(',', ' ').split() if x.strip()])
        
        if not input_genes:
            context['error'] = "請輸入至少一個基因 ID"
            return render(request, 'tool2.html', context)

        # 2. 獲取選擇的功能
        selected_features = request.POST.getlist('features')
        if not selected_features:
            context['error'] = "請至少選擇一個分析功能"
            return render(request, 'tool2.html', context)

        results_data = {}
        
        # ===== GO 項詞分析 =====
        go_types = ['go_mf', 'go_bp', 'go_cc']
        go_analysis_selected = [f for f in selected_features if f in go_types]
        
        if go_analysis_selected:
            # GO 特徵表檔案位置 (需要與 manage.py 同目錄)
            go_files = {
                'go_mf': 'GO_Term_domain_type_F_live_only.csv',  # MF = Molecular Function
                'go_bp': 'GO_Term_domain_type_P_live_only.csv',  # BP = Biological Process
                'go_cc': 'GO_Term_domain_type_C_live_only.csv'   # CC = Cellular Component
            }
            
            for go_type in go_analysis_selected:
                try:
                    feature_table = load_go_feature_table(go_files[go_type])
                    if feature_table is None:
                        context[f'warning_{go_type}'] = f"找不到 GO 特徵表檔案 {go_files[go_type]}"
                        continue
                    
                    # 執行 Fisher's Exact Test
                    results, B, D = perform_go_fisher_exact_test(input_genes, feature_table)
                    
                    # 多重檢定校正
                    results = apply_multiple_test_correction(results)
                    
                    # 轉為前端可用的格式
                    go_results = []
                    for r in results:
                        go_results.append({
                            'domain': r['go_id'],
                            'exp_str': r['exp_str'],
                            'obs_str': r['obs_str'],
                            'log2_fold_change': round(r['log2_fc'], 2),
                            'p_enrich_raw': r['p_enrich'],
                            'p_enrich_fdr': r['p_enrich_fdr'],
                            'p_enrich_bon': r['p_enrich_bonf'],
                            'p_deplete_raw': r['p_deplete'],
                            'p_deplete_fdr': r['p_deplete_fdr'],
                            'p_deplete_bon': r['p_deplete_bonf']
                        })
                    
                    results_data[f'go_{go_type}'] = {
                        'results': go_results,
                        'input_count': B,
                        'background_count': D
                    }
                    
                except Exception as e:
                    context[f'error_{go_type}'] = f"GO 分析出錯: {str(e)}"
        
        # ===== 原有的 Domain 分析 (保持相容) =====
        if 'protein_domain' in selected_features:
            try:
                domain_data = pd.read_excel("domain_test.xlsx")
            except Exception as e:
                context['error'] = f"找不到 domain_test.xlsx 或讀取錯誤: {str(e)}"
                return render(request, 'tool2.html', context)

            results = []
            B = len(input_genes) # 輸入總數
            D = TOTAL_GENES      # 背景總數

            # 遍歷所有 Domain 進行計算
            for index, row in domain_data.iterrows():
                domain_id = row["Domains"]
                # 假設 ID 欄位是以逗號分隔的基因字串
                domain_genes = set(str(row["ID"]).split(",")) 
                
                A = len(input_genes.intersection(domain_genes)) # 交集
                C = len(domain_genes) # 該 Domain 的基因數
                
                # Fisher's Exact Test
                # Enrichment (Greater)
                _, p_enrich = scipy.stats.fisher_exact([[A, C-A], [B-A, D-C-B+A]], alternative='greater')
                # Depletion (Less)
                _, p_deplete = scipy.stats.fisher_exact([[A, C-A], [B-A, D-C-B+A]], alternative='less')
                
                # 計算顯示用的數值
                observed_ratio = A / B if B > 0 else 0
                expected_ratio = C / D if D > 0 else 0
                fold_change = observed_ratio / expected_ratio if expected_ratio > 0 else 1
                # 轉換為 log2 fold change
                log2_fold_change = math.log2(fold_change) if fold_change > 0 else 0
                
                results.append({
                    "domain": domain_id,
                    "exp_str": f"{C}/{D} ({expected_ratio:.2%})",
                    "obs_str": f"{A}/{B} ({observed_ratio:.2%})",
                    "log2_fold_change": round(log2_fold_change, 2),
                    "p_enrich_raw": p_enrich,
                    "p_deplete_raw": p_deplete
                })

            # 進行多重檢定校正 (FDR & Bonferroni)
            df = pd.DataFrame(results)
            
            # Enrichment 校正
            df['p_enrich_fdr'] = multipletests(df['p_enrich_raw'], method='fdr_bh')[1]
            df['p_enrich_bon'] = multipletests(df['p_enrich_raw'], method='bonferroni')[1]
            
            # Depletion 校正
            df['p_deplete_fdr'] = multipletests(df['p_deplete_raw'], method='fdr_bh')[1]
            df['p_deplete_bon'] = multipletests(df['p_deplete_raw'], method='bonferroni')[1]

            # 轉換成 JSON 格式傳給前端
            full_data = df.to_dict(orient='records')
            results_data['protein_domain'] = {
                'results': full_data,
                'input_count': B
            }
        
        # 準備前端顯示
        if results_data:
            context['result_json'] = json.dumps(results_data)
            context['selected_features'] = selected_features
        else:
            context['error'] = "沒有可用的分析結果"

    return render(request, 'tool2.html', context)
