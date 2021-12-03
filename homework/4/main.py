from _collections import deque

# end of file character appended to every sequence
eof = "$"


# helper for reverseBwt
# def lf(index, bwt, first):
#     char = bwt[index]
#     offset = 0
#     for i in range(0, index):
#         if bwt[i] == char:
#             offset += 1
#
#     new_index = first.index(char) + offset
#     return new_index


# theoretically faster lookup
def lf(index, bwt, offsets):
    char = bwt[index]
    offset = 0
    for i in range(0, index):
        if bwt[i] == char:
            offset += 1

    return offsets[char] + offset


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


# def reverse_bwt(bwtSeq):
#     # recreate first column of table
#     first = sorted(bwtSeq)
#     seq = bwtSeq[0]
#     i = 0
#
#     offsets = dict()  # stores offset of each new character
#     for j in range(0, len(first)):
#         if first[j] not in offsets:
#             offsets[first[j]] = j
#
#     while bwtSeq[i] != eof:
#         # i = lf(i, bwtSeq, first)
#         i = lf(i, bwtSeq, offsets)
#         seq = bwtSeq[i] + seq
#
#     return seq[1:]

def reverse_bwt(bwtSeq):
    # recreate first column of table
    first = sorted(bwtSeq)
    seq = bwtSeq[0]
    i = 0

    # calculate offsets of new characters in first column
    offsets = dict()  # stores offset of each new character
    for j in range(0, len(first)):
        if first[j] not in offsets:
            offsets[first[j]] = j

    # precompute lf-mapping to speed up calculations
    # because trying to reverse bwtHomoSapiens takes AGES
    # even with the theoretically faster lookup
    lf = [0] * len(bwtSeq)
    for j in range(0, len(bwtSeq)):
        lf[j] = offsets[bwtSeq[j]]
        offsets[bwtSeq[j]] += 1  # need to update offset value as you go to prevent re-searching the array

    while bwtSeq[i] != eof:
        # i = lf(i, bwtSeq, offsets)
        i = lf[i]
        seq = bwtSeq[i] + seq

    return seq[1:]


if __name__ == '__main__':
    # testing
    # print(bwt('acaacg$'))
    # print(reverse_bwt(bwt('acaacg$')))
    # print(bwt('abc$'))
    # print(reverse_bwt(bwt('abc$')))


    # first function: bwt
    bwt_in = open("files/sample1.fa")
    bwt_in.readline()
    seq = bwt_in.read().replace('\n', '') + eof
    temp_seq = bwt(seq)

    bwt_out = open("files/bwtSample1.fa", 'w')
    bwt_out.write(">sample1\n")

    while len(temp_seq) > 0:
        if len(temp_seq) < 70:
            bwt_out.write(temp_seq + '\n')
            temp_seq = ''
        else:
            bwt_out.write(temp_seq[0:70] + '\n')
            temp_seq = temp_seq[70:]

    bwt_out.close()


    # second function: reverse_bwt
    # rev_in = open("files/bwtSample1.fa")
    rev_in = open("files/bwtHomoSapiens.fa")
    rev_in.readline()
    seq = rev_in.read().replace('\n', '')
    rev_seq = reverse_bwt(seq)

    # rev_out = open("files/reversedSample1.fa", "w")
    rev_out = open("files/rBwtSample2.fa", "w")
    rev_out.write(">sample2\n")

    while len(rev_seq) > 0:
        rev_out.write(rev_seq[0:70] + '\n')
        rev_seq = rev_seq[70:]

    rev_out.close()
