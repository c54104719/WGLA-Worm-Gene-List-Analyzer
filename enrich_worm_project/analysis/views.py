from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
import json
import pandas as pd
import scipy.stats
from statsmodels.stats.multitest import multipletests
import math
import csv
import io
from django.utils.encoding import escape_uri_path

# 導入 GO 富集分析工具函式
from .go_enrichment_utils import perform_go_fisher_exact_test, apply_multiple_test_correction, calculate_log10_p

# 設定全域變數
TOTAL_GENES = 6705
GO_TOTAL_GENES = 49164

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
                    results, B, D = perform_go_fisher_exact_test(
                        input_genes,
                        feature_table,
                        total_population=GO_TOTAL_GENES
                    )
                    
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
                    
                    # 調試：檢查結果是否正確
                    if go_results:
                        import sys
                        print(f"DEBUG: GO {go_type} results count: {len(go_results)}", file=sys.stderr)
                        print(f"DEBUG: First result keys: {list(go_results[0].keys())}", file=sys.stderr)
                    
                    results_data[go_type] = {
                        'results': go_results,
                        'input_count': B,
                        'background_count': D
                    }
                    
                except Exception as e:
                    context[f'error_{go_type}'] = f"GO 分析出錯: {str(e)}"
        
        # ===== 原有的 Domain 分析 (保持相容) =====
        if False and 'protein_domain' in selected_features:
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
            context['input_count'] = len(input_genes)
            # 儲存結果到會話以供下載
            request.session['analysis_results'] = results_data
            request.session['selected_features'] = selected_features
            # 儲存基因列表到會話以供證據頁面使用
            request.session['current_gene_list'] = list(input_genes)
        else:
            context['error'] = "沒有可用的分析結果"

    return render(request, 'tool2.html', context)


def download_analysis_results(request):
    """生成並下載分析結果 CSV 檔案"""
    analysis_type = request.GET.get('type')  # e.g., 'go_mf', 'go_bp', 'go_cc', 'protein_domain'
    result_format = request.GET.get('format')  # 'enrichment' or 'depletion'
    
    # 從會話獲取結果
    results_data = request.session.get('analysis_results', {})
    
    if not results_data or analysis_type not in results_data:
        return HttpResponse("分析結果不存在", status=404)
    
    # 取得該分析類型的結果
    analysis_result = results_data[analysis_type]
    results = analysis_result.get('results', [])
    
    if not results:
        return HttpResponse("沒有可用的結果", status=404)
    
    # 決定列名稱和要使用的 p-value 列
    if result_format == 'enrichment':
        p_value_col = 'p_enrich_raw'
        p_fdr_col = 'p_enrich_fdr'
        p_bon_col = 'p_enrich_bon'
        filename_suffix = 'enrichment'
    elif result_format == 'depletion':
        p_value_col = 'p_deplete_raw'
        p_fdr_col = 'p_deplete_fdr'
        p_bon_col = 'p_deplete_bon'
        filename_suffix = 'depletion'
    else:
        return HttpResponse("無效的格式參數", status=400)
    
    # 映射分析類型到顯示名稱
    type_names = {
        'go_mf': 'MolecularFunction',
        'go_bp': 'BiologicalProcess',
        'go_cc': 'CellularComponent',
        'protein_domain': 'ProteinDomain'
    }
    type_display = type_names.get(analysis_type, analysis_type)
    
    # 生成 CSV
    output = io.StringIO()
    writer = csv.writer(output)
    
    # 寫入表頭（與 f3 程序相同格式）
    writer.writerow([
        'Domain Name',
        'Expected Ratio',
        'Observed Ratio',
        'Fold Change',
        'P-value',
        'FDR P-value',
        'Bonferroni P-value',
        'log2(Fold Change)',
        '-log10(P-value)',
        '-log10(FDR P-value)',
        '-log10(Bonferroni P-value)'
    ])
    
    # 計算 -log10 值
    def calc_neg_log10(p_value):
        try:
            p = float(p_value)
            if p > 0:
                return -math.log10(p)
            else:
                return 0
        except:
            return 0
    
    # 寫入數據
    for row in results:
        p_raw = float(row.get(p_value_col, 1))
        p_fdr = float(row.get(p_fdr_col, 1))
        p_bon = float(row.get(p_bon_col, 1))
        log2_fc = float(row.get('log2_fold_change', 0))
        
        writer.writerow([
            row.get('domain', ''),
            row.get('exp_str', ''),
            row.get('obs_str', ''),
            f"{2 ** log2_fc:.14f}",  # 從 log2 fold change 計算 fold change
            f"{p_raw:.16e}",
            f"{p_fdr:.16e}",
            f"{p_bon:.16e}",
            f"{log2_fc:.16e}",
            f"{calc_neg_log10(p_raw):.16e}",
            f"{calc_neg_log10(p_fdr):.16e}",
            f"{calc_neg_log10(p_bon):.16e}"
        ])
    
    # 準備 HTTP 響應
    response = HttpResponse(output.getvalue(), content_type='text/csv; charset=utf-8')
    filename = f"GO_{type_display}_{filename_suffix}_result.csv"
    response['Content-Disposition'] = f'attachment; filename="{filename}"'
    
    return response


def evidence_view(request):
    """
    AJAX API - 返回特定 GO term 的證據數據（JSON格式）
    供 Modal 窗口在前端使用
    """
    from pathlib import Path
    
    try:
        # 從 GET 參數獲取 feature 和 term
        feature = request.GET.get('feature')
        term = request.GET.get('term')
        
        if not feature or not term:
            return JsonResponse({'error': '缺少 feature 或 term 參數'}, status=400)
        
        # 從會話獲取基因列表
        gene_list = request.session.get('current_gene_list', [])
        if not gene_list:
            return JsonResponse({'error': '會話中找不到基因列表，請先執行分析'}, status=400)
        
        gene_set = set(gene_list)
        
        # CSV 文件位置
        csv_dir = Path(Path(__file__).parent.parent.parent) / 'HW4PreprocessRawData'
        csv_file = csv_dir / 'f5_go_annotation_with_terms.csv'
        
        if not csv_file.exists():
            return JsonResponse({'error': f'找不到 CSV 文件: {csv_file}'}, status=500)
        
        # 讀取 CSV 並篩選數據
        evidence_list = []
        unique_genes = set()
        go_term_desc = None  # 儲存 GO term 描述
        
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # 檢查 GO ID 是否匹配
                if row.get('GO ID') == term:
                    # 第一次匹配時，取得 GO term 描述
                    if go_term_desc is None:
                        go_term_desc = row.get('GO_term', 'N/A')
                    
                    # 檢查基因是否在用戶輸入列表中
                    worm_id = row.get('WormBase ID', '')
                    symbol = row.get('Symbol', '')
                    
                    if worm_id in gene_set or symbol in gene_set:
                        evidence_list.append({
                            'gene_id': worm_id,
                            'gene_name': symbol,
                            'go_term': row.get('GO_term', ''),
                            'reference': row.get('Reference', ''),
                            'evidence_code': row.get('Evidence Code', ''),
                            'aspect': row.get('Aspect', '')
                        })
                        unique_genes.add(symbol if symbol else worm_id)
        
        # 返回 JSON 數據
        response_data = {
            'term_id': term,
            'feature': feature,
            'go_term': go_term_desc,
            'total_evidence_count': len(evidence_list),
            'unique_gene_count': len(unique_genes),
            'evidence_list': evidence_list
        }
        
        return JsonResponse(response_data)
        
    except Exception as e:
        import traceback
        return JsonResponse({'error': f'處理 CSV 文件時出錯: {str(e)}\n{traceback.format_exc()}'}, status=500)