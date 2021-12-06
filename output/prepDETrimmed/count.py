TAIR = open("TAIR10_gene_count_matrix.csv")
ARA = open("ARA11_gene_count_matrix.csv")

def genes(file):
    gene = {}
    
    skip_first = False
    for line in file:
        if not skip_first:
            skip_first = True
            continue
        stuff = line.split(",")
        gene[stuff[0]] = stuff[1:]
    
    return gene

TAIR_genes = genes(TAIR)
ARA_genes = genes(ARA)

print(f"TAIR: {len(TAIR_genes)}")
print(f"ARA: {len(ARA_genes)}")



