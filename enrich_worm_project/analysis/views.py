import os
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
    """分析工具頁面 2 - 帶功能選擇 (包括 GO 項詞分析 + FGG)"""
    context = {}

    # 設定 FGG 檔案路徑
    fgg_file_path = os.path.join(os.path.dirname(__file__), '..', '..', 'HW42', 'd12_fgg_55feature_table.csv')
    fgg_options = []
    
    # 1. 讀取 FGG 選項供前端顯示 (GET/POST 都要做)
    if os.path.exists(fgg_file_path):
        try:
            fgg_df = pd.read_csv(fgg_file_path)
            # 根據你的截圖，欄位名稱是 'Gene list'
            for _, row in fgg_df.iterrows():
                fname = row.get('Gene list', '')
                if fname:
                    fgg_options.append({'id': fname, 'name': fname})
        except Exception as e:
            print(f"Error loading FGG file: {e}")
    else:
        print(f"Warning: FGG file not found at {fgg_file_path}")
        fgg_options = [{'id': 'CSR1IP', 'name': 'CSR1IP (Test)'}]

    context['fgg_options'] = fgg_options
    
    # 2. 處理分析請求 (POST)
    if request.method == 'POST':
        # (A) 獲取使用者輸入的基因列表
        input_text = request.POST.get('gene_list', '')
        input_genes = set([x.strip() for x in input_text.replace('\n', ' ').replace(',', ' ').split() if x.strip()])
        
        if not input_genes:
            context['error'] = "請輸入至少一個基因 ID"
            return render(request, 'tool2.html', context)

        # (B) 獲取選擇的功能 (features)
        selected_features = request.POST.getlist('features')
        if not selected_features:
            context['error'] = "請至少選擇一個分析功能"
            return render(request, 'tool2.html', context)

        results_data = {}
        
        # === 區塊 1: GO 項詞分析 (保持不變) ===
        go_types = ['go_mf', 'go_bp', 'go_cc']
        go_analysis_selected = [f for f in selected_features if f in go_types]
        
        if go_analysis_selected:
            go_files = {
                'go_mf': 'GO_Term_domain_type_F_live_only.csv',
                'go_bp': 'GO_Term_domain_type_P_live_only.csv',
                'go_cc': 'GO_Term_domain_type_C_live_only.csv'
            }
            for go_type in go_analysis_selected:
                try:
                    feature_table = load_go_feature_table(go_files[go_type])
                    if feature_table is None: continue
                    
                    results, B, D = perform_go_fisher_exact_test(input_genes, feature_table, total_population=GO_TOTAL_GENES)
                    results = apply_multiple_test_correction(results)
                    
                    # 轉為前端格式
                    formatted_results = []
                    for r in results:
                        formatted_results.append({
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
                    
                    results_data[go_type] = {
                        'results': formatted_results,
                        'input_count': B, 'background_count': D
                    }
                except Exception as e:
                    print(f"GO Analysis Error: {e}")

        # === 區塊 2: FGG 分析 (新增功能) ===
        # 找出使用者選了哪些 FGG (透過比對 fgg_options 的 id)
        selected_fgg_ids = set(f for f in selected_features if any(opt['id'] == f for opt in fgg_options))
        
        if selected_fgg_ids:
            try:
                # 1. 準備 Feature Table: { "FGG_Name": set("Gene1", "Gene2"...) }
                fgg_feature_table = {}
                
                # 重新讀取 DataFrame (如果上面已經讀過，這裡其實可以優化，但為了邏輯清晰先重讀)
                if os.path.exists(fgg_file_path):
                    fgg_df = pd.read_csv(fgg_file_path)
                    
                    for _, row in fgg_df.iterrows():
                        fname = row.get('Gene list', '') # 你的截圖欄位名
                        
                        # 只處理使用者勾選的 FGG
                        if fname in selected_fgg_ids:
                            # 你的截圖欄位名是 'Gene'，且用逗號分隔
                            gene_str = str(row.get('Gene', ''))
                            # 解析基因：逗號分隔，並去除前後空白
                            gene_set = set(g.strip() for g in gene_str.split(',') if g.strip())
                            fgg_feature_table[fname] = gene_set
                
                if fgg_feature_table:
                    # 2. 直接呼叫你的通用 Fisher Test 函式
                    # 背景基因數使用 GO_TOTAL_GENES (49164)
                    fgg_results_raw, B, D = perform_go_fisher_exact_test(
                        input_genes, 
                        fgg_feature_table, 
                        total_population=GO_TOTAL_GENES 
                    )
                    
                    # 3. 校正 P-value
                    fgg_results_corrected = apply_multiple_test_correction(fgg_results_raw)
                    
                    # 4. 格式化輸出
                    fgg_formatted_results = []
                    for r in fgg_results_corrected:
                        fgg_formatted_results.append({
                            'domain': r['go_id'], # 前端用 'domain' 當作顯示名稱
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
                    
                    # 5. 存入結果，Key 取名為 'fgg'
                    results_data['fgg'] = {
                        'results': fgg_formatted_results,
                        'input_count': B,
                        'background_count': D
                    }
                    
            except Exception as e:
                context['error'] = f"FGG 分析出錯: {str(e)}"
                print(f"FGG Error: {e}")

        # === 結束 ===
        
        if results_data:
            context['result_json'] = json.dumps(results_data)
            context['selected_features'] = selected_features
            context['input_count'] = len(input_genes)
            request.session['analysis_results'] = results_data
            request.session['current_gene_list'] = list(input_genes)
        else:
            if not context.get('error'):
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
    import os
    import csv
    from django.http import JsonResponse
    
    try:
        # 1. 從 GET 參數獲取 feature 和 term
        feature = request.GET.get('feature')
        term = request.GET.get('term') # 例如 "GO:0005942"
        
        if not feature or not term:
            return JsonResponse({'error': 'Missing feature or term parameter'}, status=400)
        
        # 2. 從會話獲取基因列表 (使用者輸入的那些基因)
        # 如果是測試階段 session 常常過期，這裡加個保險，如果拿不到就用空列表避免報錯
        gene_list = request.session.get('current_gene_list', [])
        gene_set = set(gene_list)
        
        # 3. 決定要讀取哪個 CSV 檔案
        current_dir = os.path.dirname(os.path.abspath(__file__))
        
        # === 關鍵修改：路徑設定 ===
        # 假設結構是 enrich_project/HW4PreprocessRawData/f5...
        # views.py 在 enrich_worm_project/analysis/views.py
        # 所以要往上兩層 (../../) 才能看到 HW4PreprocessRawData
        csv_path = os.path.join(current_dir, '..', '..', 'HW4PreprocessRawData', 'f5_go_annotation_with_terms.csv')
        csv_path = os.path.normpath(csv_path)

        # 除錯：印出路徑看看對不對
        print(f"DEBUG: Loading GO Evidence from: {csv_path}")
        
        if not os.path.exists(csv_path):
             print(f"Error: File not found at {csv_path}")
             return JsonResponse({'error': f'Evidence file not found at: {csv_path}'}, status=500)

        # 4. 讀取 CSV 並篩選數據
        evidence_list = []
        unique_genes = set()
        go_term_desc = "N/A" # 預設描述
        
        try:
            with open(csv_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                
                for row in reader:
                    # A. 檢查是否為目標 Term (GO ID)
                    if row.get('GO ID') == term:
                        
                        # 順便抓一下 GO Term 的描述文字
                        if go_term_desc == "N/A":
                            go_term_desc = row.get('GO_term', 'N/A')

                        # B. 檢查這行基因是否在使用者輸入的列表中 (只顯示交集)
                        wbid = row.get('WormBase ID', '').strip()
                        symbol = row.get('Symbol', '').strip()
                        
                        # 如果 WBID 或 Symbol 任一個有中，就加進去
                        if (wbid in gene_set) or (symbol in gene_set):
                            evidence_list.append({
                                'gene_name': symbol,   
                                'gene_id': wbid,   
                                'go_term': row.get('GO_term', ''),
                                'evidence_code': row.get('Evidence Code', ''),
                                'reference': row.get('Reference', '')
                            })
                            unique_genes.add(wbid)
        except Exception as e:
            print(f"Error reading CSV: {e}")
            return JsonResponse({'error': f'Error reading CSV file: {str(e)}'}, status=500)

        # 5. 回傳 JSON
        response_data = {
            'term_id': term,
            'go_term': go_term_desc,
            'feature': feature,
            'list_name': request.session.get('list_name', 'User Input List'),
            'input_gene_count': len(gene_set), 
            'associated_gene_count': len(unique_genes),
            'total_evidence_count': len(evidence_list),
            'evidence_list': evidence_list
        }
        
        return JsonResponse(response_data)
        
    except Exception as e:
        import traceback
        print(traceback.format_exc())
        return JsonResponse({'error': f"Backend Logic Error: {str(e)}"}, status=500)