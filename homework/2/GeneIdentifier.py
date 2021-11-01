mod = 43507093  # address distance between given FASTA file and gene file
gene_file = open("Chr21.txt")
genes = []
for line in gene_file:
    genes.append(line.strip().split('\t'))


def identify(isl_end):
    ret_genes = []
    for gene in genes:
        start = int(gene[0])
        if 500 >= start - (isl_end + mod) >= 0:  # find gene within 500 bp of scaled island end
            # print("start:", start, "isl_end:", isl_end, "isl_end+mod:", isl_end+mod, "diff:", start - (isl_end+mod))
            ret_genes.append(gene)

    return ret_genes

# if __name__ == '__main__':
#     identify("Chr21.txt", "islands.txt")
