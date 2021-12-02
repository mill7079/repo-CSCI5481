from _collections import deque

# end of file character appended to every sequence
eof = '$'


# helper for reverseBwt
def lf(index, bwt, first):
    char = bwt[index]
    offset = 0
    for i in range(0, index):
        if bwt[i] == char:
            offset += 1

    new_index = first.index(char) + offset
    return new_index


# learning from hw2 mistakes and using a deque
def bwt(seq):
    table = deque()

    for i in range(0, len(seq)):
        table.append(seq[i:] + seq[:i])

    sortedTable = sorted(table)
    bwtSeq = ''
    for row in sortedTable:
        bwtSeq += row[-1]

    return bwtSeq


def reverseBwt(bwtSeq):
    # recreate first column of table
    first = sorted(bwtSeq)
    seq = bwtSeq[0]
    i = 0

    while bwtSeq[i] != eof:
        i = lf(i, bwtSeq, first)
        seq = bwtSeq[i] + seq

    return seq[1:]


if __name__ == '__main__':
    print(bwt('acaacg$'))
    print(reverseBwt(bwt('acaacg$')))
    # print(bwt('abc$'))
    # print(reverseBwt(bwt('abc$')))
