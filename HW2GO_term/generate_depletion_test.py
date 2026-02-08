import csv
import random

# 讀取所有基因
all_genes = set()
with open('../enrich_worm_project/GO_Term_domain_type_C_live_only.csv', 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        genes = row['ID'].strip('"').split(',')
        for gene in genes:
            all_genes.add(gene.strip())

print(f'Total genes in database: {len(all_genes)}')

# 隨機選 32 個基因
random.seed(42)  # 固定種子以便重複
test_genes = sorted(random.sample(list(all_genes), 32))

# 儲存到文件
with open('test_genes_depletion.txt', 'w') as f:
    for gene in test_genes:
        f.write(gene + '\n')

print('Generated test gene set:')
for gene in test_genes:
    print(gene)
