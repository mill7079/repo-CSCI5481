from GeneIdentifier import identify
from HMM import *
from sys import argv

if __name__ == '__main__':

    if len(argv) != 2:
        print("wrong number of arguments")
        exit(-1)

    _, seqfile = argv
    hmm = HMM()

    file = open(seqfile)
    file.readline()
    seq = ''
    for line in file:  # concatenate parts of fasta file
        seq += line.strip()

    # find optimal path using viterbi algorithm
    path = hmm.viterbi(seq)

    cg_count = 0  # length of CpG island
    start_idx = 0  # start index of CpG island
    island_count = 0  # number of CpG islands found so far
    empty_count = 0

    # print("CpG islands.txt longer than 200 bp:")
    # for i in range(1, len(path)):
    #     state = path[i]
    #     if state == '':
    #         empty_count += 1
    #     else:
    #         if state[1] == '-' and cg_count > 0:
    #             if cg_count > 200:
    #                 island_count += 1
    #                 print("CpG Island " + str(island_count) + ":", cg_count, "bp (" + str(start_idx), '-', str(i) + ')')
    #             cg_count = 0
    #         elif state[1] == '+':
    #             if cg_count == 0:
    #                 start_idx = i
    #             cg_count += 1

    islands = []
    with_gene = 0
    
    # find islands and potential genes
    for i in range(1, len(path)):
        state = path[i]
        if state == '':  # meant for handling errors; I don't think this actually needs to be here anymore
            empty_count += 1
        else:
            if state[1] == '-' and cg_count > 0:
                if cg_count > 200:
                    island_count += 1
                    genes = identify(i)
                    if len(genes) > 0:  # island has downstream gene
                        with_gene += 1

                    island = [island_count, start_idx, i, cg_count, genes]  # island num, start, end, length, genes
                    islands.append(island)
                cg_count = 0
            elif state[1] == '+':
                if cg_count == 0:
                    start_idx = i
                cg_count += 1

    print("Total CpG islands found:", str(island_count)+';', with_gene, "out of", island_count,
          "islands are followed by a coding region")

    # print formatted islands
    for island in islands:
        print("CpG Island " + str(island[0]) + ":", island[3], "bp (" + str(start_idx), '-', str(i) + '),', end='')
        if len(island[4]) > 0:
            for gene in island[4]:
                print(",", gene[3], "(" + gene[2] + ")", end='')
        else:
            print("No gene", end='')
        print()
