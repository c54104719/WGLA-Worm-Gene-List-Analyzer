import csv
import random

# 讀取 GO feature table
go_genes = {}
with open('../enrich_worm_project/GO_Term_domain_type_C_live_only.csv', 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        go_id = row['Domains'].strip()
        genes = set(row['ID'].strip('"').split(','))
        genes = {g.strip() for g in genes if g.strip()}
        go_genes[go_id] = genes

# 找最大的 GO term（最多基因）
largest_go = max(go_genes.items(), key=lambda x: len(x[1]))
print(f"Largest GO term: {largest_go[0]} with {len(largest_go[1])} genes")

# 收集所有基因
all_genes = set()
for genes in go_genes.values():
    all_genes.update(genes)

print(f"Total genes: {len(all_genes)}")

# 找第二大、第三大的 GO terms 來排除
sorted_go = sorted(go_genes.items(), key=lambda x: len(x[1]), reverse=True)
print("\nTop 10 largest GO terms:")
for i, (go_id, genes) in enumerate(sorted_go[:10]):
    print(f"  {i+1}. {go_id}: {len(genes)} genes")

# 選擇要排除的 GO terms（選前 3 個最大的）
exclude_gos = [sorted_go[0][0], sorted_go[1][0], sorted_go[2][0]]
print(f"\nExcluding genes from: {exclude_gos}")

# 找不在這些 GO terms 中的基因
excluded_genes = set()
for go_id in exclude_gos:
    excluded_genes.update(go_genes[go_id])

available_genes = all_genes - excluded_genes
print(f"Available genes (not in excluded GO terms): {len(available_genes)}")

# 隨機選 32 個
if len(available_genes) >= 32:
    random.seed(100)
    test_genes = sorted(random.sample(list(available_genes), 32))
    
    print("\nGenerated test gene set (should show depletion for excluded GO terms):")
    for gene in test_genes:
        print(gene)
    
    # 儲存
    with open('test_genes_depletion_specific.txt', 'w') as f:
        for gene in test_genes:
            f.write(gene + '\n')
else:
    print("Not enough genes available!")
