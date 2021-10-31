from HMM import *
from sys import argv

if __name__ == '__main__':

    # print(hmm.emit('A+', 'A'))
    # print(hmm.emit('T+', 'A'))
    # hmm.viterbi('ATACGACA')
    # hmm.viterbi("CGCCG")

    if len(argv) != 2:
        print("wrong number of arguments")
        exit(-1)

    _, seqfile = argv
    hmm = HMM()

    file = open(seqfile)
    file.readline()
    seq = ''
    for line in file:
        seq += line.strip()

    path = hmm.viterbi(seq)

    cg_count = 0
    start_idx = 0
    island_count = 0
    empty_count = 0

    print("CpG islands.txt longer than 200 bp:")
    for i in range(1, len(path)):
        if i % 100000 == 0:
            print("step")

        state = path[i]
        if state == '':
            empty_count += 1
        else:
            if state[1] == '-' and cg_count > 0:
                if cg_count > 200:
                    island_count += 1
                    print("CpG Island " + str(island_count) + ":", cg_count, "bp (" + str(start_idx), '-', str(i) + ')')
                cg_count = 0
            elif state[1] == '+':
                if cg_count == 0:
                    start_idx = i
                cg_count += 1
