import pandas as pd
import scipy.stats
from statsmodels.stats.multitest import multipletests

# 1. 設定基礎參數
# 酵母菌總基因數 (Total Population)
TOTAL_GENES = 6705 

# 讀取資料 (請確保 domain_test.xlsx 檔案存在且格式正確)
domain_data = pd.read_excel("domain_test.xlsx")

# 輸入的基因列表 (Input List)
input_list = [ 'YBR021W', 'YKL051W', 'YJR009C', 'YMR109W', 'YGR192C', 'YJL052W', 'YKR014C', 'YOL103W', 'YGR191W', 'YLR332W', 'YKL220C', 'YAL038W', 'YCR012W', 'YIL033C', 'YDR348C', 'YLR413W', 'YNL192W', 'YHR005C', 'YIR019C', 'YJL005W', 'YGR254W', 'YHR174W', 'YLL052C', 'YPR192W', 'YDR050C', 'YML132W', 'YER020W', 'YDL194W', 'YLR452C', 'YJL129C', 'YKL126W', 'YKL209C', 'YGL115W', 'YML116W', 'YLR081W', 'YDR122W', 'YLR096W', 'YFR012W', 'YDR040C', 'YGR152C', 'YGL167C', 'YJR152W', 'YPL036W', 'YER056C', 'YCR075C', 'YOR212W', 'YJR086W', 'YLR229C', 'YKR039W', 'YDR164C', 'YOR317W', 'YJL085W', 'YNL194C', 'YER091C', 'YLR303W', 'YMR307W', 'YNL183C', 'YLR342W', 'YFL014W', 'YHR135C', 'YNL154C', 'YIL009W', 'YMR011W', 'YNL084C', 'YLL043W', 'YPR165W', 'YCR098C', 'YCR037C', 'YCL048W', 'YGR032W', 'YCL073C', 'YCR010C', 'YFR022W', 'YCR021C', 'YJR100C', 'YJR125C', 'YMR246W', 'YFL005W', 'YGL045W', 'YPL204W', 'YGR068C', 'YDR497C', 'YGR121C', 'YDR522C', 'YAL030W', 'YDR090C', 'YDR039C', 'YOR018W', 'YDR153C', 'YHR094C', 'YDR345C', 'YHR092C', 'YMR306W', 'YPR159W', 'YLR206W', 'YDR011W', 'YKL203C', 'YER166W', 'YLR214W', 'YLL007C', 'YPL232W', 'YKR093W', 'YCR024C-A', 'YNR002C', 'YDR103W', 'YOR153W', 'YHL016C', 'YJR160C', 'YDL247W', 'YGL208W', 'YKR059W', 'YJL138C', 'YAL005C', 'YLL024C', 'YJR066W', 'YBL069W', 'YDL229W', 'YHR214W', 'YDR212W', 'YKL217W', 'YBR196C', 'YIL034C', 'YKR105C', 'YKR106W', 'YJR040W', 'YLR305C', 'YDR158W', 'YBR054W', 'YBR069C', 'YOR322C', 'YDL161W', 'YDR038C', 'YBR008C', 'YBL106C', 'YBL042C', 'YBR016W', 'YBR043C', 'YBR086C', 'YMR186W', 'YBR294W', 'YBR295W', 'YBR296C', 'YHR201C', 'YHR073W', 'YHL044W', 'YHL040C', 'YHR048W', 'YHR114W', 'YCL040W', 'YOL122C', 'YOL020W', 'YMR058W', 'YDR208W', 'YDR342C', 'YMR008C', 'YJL198W', 'YCR088W', 'YAR031W', 'YAR033W', 'YAL056W', 'YMR183C', 'YIL147C', 'YDR536W', 'YER155C', 'YER123W', 'YER060W', 'YER118C', 'YER120W', 'YER143W', 'YER145C', 'YDR461W', 'YER185W', 'YNL145W', 'YBL060W', 'YJL093C', 'YBL061C', 'YHR042W', 'YKL187C', 'YIL120W', 'YIL105C', 'YIL088C', 'YNL322C', 'YIL047C', 'YMR319C', 'YOR188W', 'YLR353W', 'YDR420W', 'YNL291C', 'YPL092W', 'YNL142W', 'YNL323W', 'YJL100W', 'YFL050C', 'YGR197C', 'YGR198W', 'YJL170C', 'YJL145W', 'YJL058C', 'YJR059W', 'YNL294C', 'YNL293W', 'YPR124W', 'YGR217W', 'YGR224W', 'YGR055W', 'YGL233W', 'YOR371C', 'YOR328W', 'YBL105C', 'YGR060W', 'YGR213C', 'YGR281W', 'YGL186C', 'YBL091C-A', 'YGL084C', 'YGL053W', 'YGL051W', 'YGR041W', 'YGR065C', 'YGR138C', 'YGR266W', 'YNR049C', 'YNR047W', 'YNR060W', 'YNR070W', 'YNL275W', 'YCR004C', 'YNL180C', 'YNL173C', 'YCR027C', 'YNL065W', 'YNL047C', 'YML052W', 'YOR008C', 'YCL058C', 'YNL093W', 'YDR276C', 'YCR017C', 'YKR055W', 'YCR028C', 'YDR463W', 'YPL058C', 'YOR047C', 'YJR005W', 'YDR459C', 'YBR023C', 'YML072C', 'YMR212C', 'YMR215W', 'YDR160W', 'YDR104C', 'YML006C', 'YMR192W', 'YMR017W', 'YMR243C', 'YDL222C', 'YLR120C', 'YGR014W', 'YMR063W', 'YOR316C', 'YMR251W-A', 'YKL094W', 'YKL007W', 'YLR237W', 'YLR138W', 'YIR006C', 'YLR422W', 'YPR156C', 'YLR414C', 'YJR054W', 'YPR194C', 'YPR201W', 'YDR309C', 'YER177W', 'YLL010C', 'YLL028W', 'YLR019W', 'YLR020C', 'YML047C', 'YOL019W', 'YOL078W', 'YOL130W', 'YOL158C', 'YOR011W', 'YOR049C', 'YOR071C', 'YOR306C', 'YOR378W', 'YOR381W', 'YOR390W', 'YPL274W', 'YPL279C', 'YPR032W', 'YOR104W', 'YAR050W', 'YDR033W', 'YER060W-A', 'YMR031C', 'YPR149W', 'YOR171C', 'YOR273C', 'YDL138W', 'YLR092W', 'YOL152W', 'YPL249C', 'YDR055W', 'YDR384C', 'YDL035C', 'YLL061W', 'YLR130C', 'YOL002C', 'YDL019C', 'YOL113W', 'YDL012C', 'YDR093W', 'YML125C', 'YGR031C-A', 'YOR129C', 'YOR030W', 'YDR129C', 'YEL047C', 'YKL035W', 'YKL046C', 'YDR099W', 'YKL105C', 'YBL007C', 'YLR058C', 'YKL196C', 'YKR100C', 'YNL090W', 'YDR414C', 'YPL158C', 'YBL099W', 'YAL014C', 'YAL022C', 'YAR042W', 'YBR264C', 'YBR129C', 'YLR025W', 'YOR275C', 'YMR120C', 'YLR109W', 'YBL085W', 'YBR299W', 'YHR096C', 'YNL257C', 'YDR343C', 'YER125W', 'YNL209W', 'YIR038C', 'YNL271C', 'YOR327C', 'YJL171C', 'YJL158C', 'YJR065C', 'YLR432W', 'YGL082W', 'YGR086C', 'YGR122W', 'YGR256W', 'YGR292W', 'YPL265W', 'YNL231C', 'YIL118W', 'YLL050C', 'YML128C', 'YHL007C', 'YPR075C', 'YDL124W', 'YDL223C', 'YOL109W', 'YPL004C', 'YDR032C', 'YOR161C', 'YDL135C', 'YOR086C', 'YBL029C-A', 'YOL009C', 'YPL240C', 'YGR131W', 'YML003W', 'YML002W', 'YBL075C', 'YLL053C', 'YGR241C', 'YHR126C', 'YHR161C', 'YJR001W', 'YIL140W', 'YAR066W', 'YCL025C', 'YKL064W', 'YNL128W', 'YDR506C', 'YDL192W', 'YOR353C', 'YER067W', 'YBR005W', 'YER093C', 'YDR210W', 'YIR039C', 'YJL212C', 'YNL006W', 'YFL054C', 'YFL051C', 'YFR029W', 'YGL160W', 'YJL156C', 'YNL012W', 'YNL087W', 'YLR443W', 'YLR047C', 'YGL108C', 'YKL060C', 'YNL283C', 'YNL279W', 'YNL033W', 'YNL019C', 'YJR058C', 'YPL056C', 'YDL185W', 'YGR167W', 'YMR238W', 'YMR034C', 'YLR194C', 'YLR219W', 'YDL137W', 'YLR343W', 'YMR244C-A', 'YOR348C', 'YDR373W', 'YER103W', 'YCR091W', 'YCR030C', 'YCR094W', 'YKR050W', 'YPR171W', 'YNR048W', 'YGL255W', 'YER008C', 'YGR143W', 'YBR207W', 'YIL043C', 'YHL047C', 'YHR050W', 'YHR155W', 'YAL026C', 'YEL065W', 'YGR009C', 'YIL121W', 'YIL048W', 'YFL047W', 'YMR279C', 'YMR266W', 'YMR006C', 'YML087C', 'YIR028W', 'YMR032W', 'YOR301W', 'YLR241W', 'YLR411W', 'YLL005C', 'YOL011W', 'YOR192C', 'YOR384W', 'YPL180W', 'YLR034C', 'YOL084W', 'YLR046C', 'YLL051C', 'YMR162C', 'YGR212W', 'YCL048W-A', 'YDR524C-B', 'YPR055W', 'YIL171W', 'YLR299W', 'YBR068C', 'YLL055W', 'YNL243W', 'YBL037W', 'YHR006W', 'YAR027W', 'YEL017C-A', 'YFL041W', 'YJR151C', 'YGR221C', 'YDR261C', 'YGR026W', 'YDR144C', 'YMR086W', 'YMR068W', 'YLR373C', 'YLR187W', 'YLR004C', 'YOL132W', 'YPL176C', 'YLR121C', 'YLR084C', 'YLR262C']
input_set = set(input_list) # 轉成 set 加速搜尋

# 2. 準備儲存結果的列表
results = []

# 3. 迴圈計算每個 Domain (計算 A, B, C, D)
for index, row in domain_data.iterrows():
    # 取得該 domain 的基因列表
    domain_genes = set(row["ID"].split(","))
    
    # 計算 Fisher's Exact Test 所需的 A, B, C, D
    # A (T): 交集數量 (Intersection)
    A = len(input_set.intersection(domain_genes))
    # B (S): 輸入基因總數 (Input Size)
    B = len(input_set)
    # C (G): 該 Domain 的基因總數 (Domain Size)
    C = len(domain_genes) # 也可以用 row["number"]
    # D (F): 生物體總基因數 (Total Population)
    D = TOTAL_GENES

    # 進行 Fisher's Exact Test
    # Enrichment (Greater)
    _, p_val_enrich = scipy.stats.fisher_exact([[A, C-A], [B-A, D-C-B+A]], alternative='greater')
    # Depletion (Less)
    _, p_val_deplete = scipy.stats.fisher_exact([[A, C-A], [B-A, D-C-B+A]], alternative='less')

    # 計算比率與倍數 (Fold Change)
    observed_ratio = A / B if B > 0 else 0
    expected_ratio = C / D if D > 0 else 0
    fold_change = observed_ratio / expected_ratio if expected_ratio > 0 else 0

    # 格式化字串 (如老師要求的: "18/6705 (0.27%)")
    exp_str = f"{C}/{D} ({expected_ratio:.2%})"
    obs_str = f"{A}/{B} ({observed_ratio:.2%})"

    results.append({
        "Domain_id": row["Domains"],
        "Expected_Ratio_Str": exp_str,
        "Observed_Ratio_Str": obs_str,
        "Fold_Change": fold_change,
        "P_Value_Enrich": p_val_enrich,
        "P_Value_Deplete": p_val_deplete,
        "A": A, "B": B, "C": C, "D": D  # 保留原始數字方便檢查
    })

# 4. 轉成 DataFrame
df_results = pd.DataFrame(results)

# 5. 多重檢定校正 (Multiple Testing Correction)
# 分別針對 Enrichment 和 Depletion 做校正
cut_off = 0.01

# --- 校正 Enrichment P-values ---
_, fdr_enrich, _, _ = multipletests(df_results["P_Value_Enrich"], alpha=cut_off, method="fdr_bh")
_, bon_enrich, _, _ = multipletests(df_results["P_Value_Enrich"], alpha=cut_off, method="bonferroni")
df_results["FDR_Enrich"] = fdr_enrich
df_results["Bonferroni_Enrich"] = bon_enrich

# --- 校正 Depletion P-values ---
_, fdr_deplete, _, _ = multipletests(df_results["P_Value_Deplete"], alpha=cut_off, method="fdr_bh")
_, bon_deplete, _, _ = multipletests(df_results["P_Value_Deplete"], alpha=cut_off, method="bonferroni")
df_results["FDR_Deplete"] = fdr_deplete
df_results["Bonferroni_Deplete"] = bon_deplete


# 6. 整理並輸出檔案 (根據老師要求)

# --- 檔案 1: Enrichment 底層全結果 (Base Result File) ---
cols_enrich = ["Domain_id", "Expected_Ratio_Str", "Observed_Ratio_Str", "Fold_Change", "P_Value_Enrich", "FDR_Enrich", "Bonferroni_Enrich"]
df_enrich_base = df_results[cols_enrich].copy()
df_enrich_base.sort_values("P_Value_Enrich", inplace=True)
df_enrich_base.to_excel("Enrichment_Base_Result.xlsx", index=False)

# --- 檔案 2: Depletion 底層全結果 (Base Result File) ---
cols_deplete = ["Domain_id", "Expected_Ratio_Str", "Observed_Ratio_Str", "Fold_Change", "P_Value_Deplete", "FDR_Deplete", "Bonferroni_Deplete"]
df_deplete_base = df_results[cols_deplete].copy()
df_deplete_base.sort_values("P_Value_Deplete", inplace=True)
df_deplete_base.to_excel("Depletion_Base_Result.xlsx", index=False)

# --- 檔案 3: 篩選後的 Enrichment 結果 (Filtered Result) ---
# 這裡範例使用 FDR < 0.01 篩選
df_enrich_filtered = df_enrich_base[df_enrich_base["FDR_Enrich"] < cut_off]
df_enrich_filtered.to_excel("Enrichment_Filtered_FDR0.01.xlsx", index=False)

# --- 檔案 4: 篩選後的 Depletion 結果 (Filtered Result) ---
df_deplete_filtered = df_deplete_base[df_deplete_base["FDR_Deplete"] < cut_off]
df_deplete_filtered.to_excel("Depletion_Filtered_FDR0.01.xlsx", index=False)


# 7. 顯示簡單結果
print(f"Total Domains Processed: {len(df_results)}")
print(f"Enrichment Significant (FDR < {cut_off}): {len(df_enrich_filtered)}")
print(f"Depletion Significant (FDR < {cut_off}): {len(df_deplete_filtered)}")
print("\n--- Enrichment Top 5 ---")
print(df_enrich_base.head(5))