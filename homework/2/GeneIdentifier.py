mod = 43507093  # address distance between given FASTA file and gene file
gene_file = open("Chr21.txt")
genes = []
for line in gene_file:
    genes.append(line.strip().split('\t'))


def identify(isl_start, isl_end):
    ret_genes = []
    for gene in genes:
        start = int(gene[0])
        end = int(gene[1])
        if 500 >= start - (isl_end + mod) >= 0 and gene[2] == '+':  # find gene within 500 bp of scaled island end
            print("start:", start, "isl_end:", isl_end, "isl_end+mod:", isl_end+mod, "diff:", start - (isl_end+mod), "dir:", gene[2])
            ret_genes.append(gene)
        if 500 >= end - (isl_start + mod) >= 0 and gene[2] == '-':
            print("end:", end, "isl_start:", isl_start, "isl_start+mod:", isl_start+mod, "diff:", end - (isl_start+mod), "dir:", gene[2])
            ret_genes.append(gene)

    return ret_genes
