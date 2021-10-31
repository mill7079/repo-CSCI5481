def identify(gene_file, island_file):
    g_file = open(gene_file)
    genes = []
    for line in g_file:
        genes.append(line.split('\t'))
    print(genes)


if __name__ == '__main__':
    identify("Chr21.txt", )