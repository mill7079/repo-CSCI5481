mod = 43507093  # address distance between given FASTA file and gene file
gene_file = open("Chr21.txt")
genes = []
for line in gene_file:
    genes.append(line.strip().split('\t'))


# def identify(island_file):
#     # g_file = open(gene_file)
#     # genes = []
#     # for line in g_file:
#     #     genes.append(line.strip().split('\t'))
#     # # print(genes)
#
#     i_file = open(island_file)
#     islands = []
#     i_file.readline()
#     for line in i_file:
#         bpsplit = line.strip().split(' bp ')
#         island = bpsplit[1][1:-1].split(' - ')
#         isl = [int(island[0]), int(island[1])]
#         print(isl)
#         islands.append(isl)
#     # print(islands)
#
#     mod = 43507093  # address distance between given FASTA file and gene file
#     for gene in genes:
#         start = int(gene[0])
#         for isl in islands:
#             if start - (isl[1]+mod) <= 500:
#                 print("CpG island length", isl[1]-isl[0], " nucleotides:", isl[0]+mod, '-', isl[1]+mod, " Gene:", gene[3], gene[2])
#                 break


def identify(isl_end):
    ret_genes = []
    for gene in genes:
        start = int(gene[0])
        if (isl_end+mod) - start <= 500:  # find gene within 500 bp of scaled island end
            ret_genes.append(gene)

    return ret_genes

# if __name__ == '__main__':
#     identify("Chr21.txt", "islands.txt")
