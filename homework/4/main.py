from _collections import deque

# end of file character appended to every sequence
eof = "$"

# print actual locations in chromosome
chr_offset = 25000000


# helper to parse fastq file into arrays of reads
def parse_fastq(file):
    fastq = open(file)
    line_count = 0
    reads = []
    head = ""
    for line in fastq:
        if line_count % 4 == 0:
            head = line.strip().split("/")[0]
        elif line_count % 4 == 1:
            reads.append((head, line.strip()))
        line_count += 1

    return reads

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
# def lf(index, bwt, offsets):
#     char = bwt[index]
#     offset = 0
#     for i in range(0, index):
#         if bwt[i] == char:
#             offset += 1
#
#     return offsets[char] + offset

# helper for reverse_bwt, fm indexing
# finds offsets of new characters in first column
def calc_offsets(col):
    offsets = dict()  # stores offset of each new character
    for j in range(0, len(col)):
        if col[j] not in offsets:
            offsets[col[j]] = j

    return offsets


# helper for reverse_bwt, fm indexing
# precomputes lf array
def calc_lf(bwt_seq, offsets):
    lf = [0] * len(bwt_seq)
    for j in range(0, len(bwt_seq)):
        lf[j] = offsets[bwt_seq[j]]
        offsets[bwt_seq[j]] += 1  # need to update offset value as you go to prevent re-searching the array

    return lf


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
    sa = [0] * len(bwtSeq)  # suffix array
    count = len(bwtSeq) - 1
    sa[0] = count

    # calculate offsets of new characters in first column
    # offsets = dict()  # stores offset of each new character
    # for j in range(0, len(first)):
    #     if first[j] not in offsets:
    #         offsets[first[j]] = j
    offsets = calc_offsets(first)

    # precompute lf-mapping to speed up calculations
    # because trying to reverse bwtHomoSapiens takes AGES
    # even with the theoretically faster lookup
    # lf = [0] * len(bwtSeq)
    # for j in range(0, len(bwtSeq)):
    #     lf[j] = offsets[bwtSeq[j]]
    #     offsets[bwtSeq[j]] += 1  # need to update offset value as you go to prevent re-searching the array
    lf = calc_lf(bwtSeq, offsets)

    while bwtSeq[i] != eof:
        # i = lf(i, bwtSeq, offsets)
        count -= 1
        i = lf[i]
        seq = bwtSeq[i] + seq
        sa[i] = count


    # return seq[1:]
    return seq[1:], sa


# full fm index; keeps entire suffix array in memory
def full_index_fm(bwt_seq, all_reads):
    seq, sa = reverse_bwt(bwt_seq)  # original sequence and suffix array
    first = sorted(bwt_seq)  # first column of matrix

    offsets = calc_offsets(first)
    lf = calc_lf(bwt_seq, offsets)

    # calculate array of row indices
    init_rows = []
    for i in range(0, len(bwt_seq)):
        init_rows.append(i)

    matches = dict()

    for reads in all_reads:
        progress_count = 0
        for read in reads:  # align each short read
            if progress_count % 10 == 0:
                print("progress:", progress_count, "out of", len(reads), "reads")
            rev_read = read[1][::-1]
            rows = init_rows  # contains row indices for currently valid rows

            # while len(rows) > 1:
            while len(rev_read) > 0:
                char = rev_read[0]
                rev_read = rev_read[1:]
                temp_rows = []
                for row in rows:  # find valid rows (row ends with next char)
                    if bwt_seq[row] == char:
                        temp_rows.append(row)

                # map valid rows to rows in first
                rows = []
                for row in temp_rows:
                    rows.append(lf[row])

                if len(rows) < 1:
                    # print("no matching rows")
                    break

            if len(rows) > 0:
                # print("matching:", rows)
                indices = []
                for row in rows:
                    indices.append(sa[row] + chr_offset)

                if read[0] in matches:
                    matches[read[0]].append(indices)
                else:
                    matches[read[0]] = [indices]

            progress_count += 1

    return matches


# hybrid fm index; keeps partial suffix array in memory
def hybrid_fm(bwt_seq, all_reads):
    seq, full_sa = reverse_bwt(bwt_seq)  # original sequence and suffix array
    # only store portion of SA (but it's auto-calculated in reverse_bwt so....pretend this happened earlier)
    sa = [0] * (len(full_sa) // 3 + 1)
    for i in range(0, len(sa)):  # shrink by 2/3
        sa[i] = full_sa[i*3]

    first = sorted(bwt_seq)  # first column of matrix

    offsets = calc_offsets(first)
    lf = calc_lf(bwt_seq, offsets)

    # calculate array of row indices
    init_rows = []
    for i in range(0, len(bwt_seq)):
        init_rows.append(i)

    matches = dict()

    for reads in all_reads:
        progress_count = 0
        for read in reads:  # align each short read
            if progress_count % 10 == 0:
                print("progress hybrid:", progress_count, "out of", len(reads), "reads")
            rev_read = read[1][::-1]
            rows = init_rows  # contains row indices for currently valid rows

            # while len(rows) > 1:
            while len(rev_read) > 0:
                char = rev_read[0]
                rev_read = rev_read[1:]
                temp_rows = []
                for row in rows:  # find valid rows (row ends with next char)
                    if bwt_seq[row] == char:
                        temp_rows.append(row)

                # map valid rows to rows in first
                rows = []
                for row in temp_rows:
                    rows.append(lf[row])

                if len(rows) < 1:
                    # print("no matching rows")
                    break

            if len(rows) > 0:
                # print("matching:", rows)
                indices = []
                for row in rows:
                    r = row
                    steps = 0
                    while r % 3 != 0:
                        # print(r)
                        r = lf[r]
                        steps += 1
                    indices.append(sa[r // 3] + steps + chr_offset)

                if read[0] in matches:
                    matches[read[0]].append(indices)
                else:
                    matches[read[0]] = [indices]

            progress_count += 1

    return matches


if __name__ == '__main__':
    # testing
    # print(bwt('acaacg$'))
    # print(reverse_bwt(bwt('acaacg$')))
    # print(bwt('abc$'))
    # print(reverse_bwt(bwt('abc$')))


    # first function: bwt
    # bwt_in = open("files/sample1.fa")
    # bwt_in.readline()
    # seq = bwt_in.read().replace('\n', '') + eof
    # temp_seq = bwt(seq)
    #
    # bwt_out = open("files/bwtSample1.fa", 'w')
    # bwt_out.write(">sample1\n")
    #
    # while len(temp_seq) > 0:
    #     if len(temp_seq) < 70:
    #         bwt_out.write(temp_seq + '\n')
    #         temp_seq = ''
    #     else:
    #         bwt_out.write(temp_seq[0:70] + '\n')
    #         temp_seq = temp_seq[70:]
    #
    # bwt_out.close()


    # second function: reverse_bwt
    # # rev_in = open("files/bwtSample1.fa")
    # rev_in = open("files/bwtHomoSapiens.fa")
    # rev_in.readline()
    # seq = rev_in.read().replace('\n', '')
    # rev_seq = reverse_bwt(seq)[0]
    #
    # # rev_out = open("files/reversedSample1.fa", "w")
    # rev_out = open("files/rBwtSample2.fa", "w")
    # rev_out.write(">sample2\n")
    #
    # while len(rev_seq) > 0:
    #     rev_out.write(rev_seq[0:70] + '\n')
    #     rev_seq = rev_seq[70:]
    #
    # rev_out.close()


    # third function: full index
    reads = [("head", "aba"), ("head2", "bba")]
    reads2 = [("head", "aba"), ("head2", "bba")]
    chr_bwt_seq = "abba$aa"
    # parse fastq files for reads
    # reads = parse_fastq("files/SRR089545_1.fq")
    # reads2 = parse_fastq("files/SRR089545_2.fq")

    # extract chr Y BWT
    # chry_bwt = open("files/bwtChrYnew(25M-26M).fa")
    # chry_bwt.readline()
    # chr_bwt_seq = chry_bwt.read().replace("\n", "")

    # run fm indexing
    matches = full_index_fm(chr_bwt_seq, [reads, reads2])

    # write matches to file
    full_out = open("files/mapping_fullindexFM.txt", "w")
    full_out.write(str(len(matches)) + " out of " + str(len(reads) + len(reads2)) + " short read pairs map on chrY:20M~40M:\n")
    for match in matches:
        full_out.write(match + " " + str(matches[match][0]))
        try:
            full_out.write(" " + str(matches[match][1]))
        except IndexError:
            pass
        full_out.write("\n")

    full_out.close()


    # fourth function: hybrid fm
    matches = hybrid_fm(chr_bwt_seq, [reads, reads2])

    hybrid_out = open("files/mapping_FMhybrid.txt", "w")
    hybrid_out.write(str(len(matches)) + " out of " + str(len(reads) + len(reads2)) + " short read pairs map on chrY:20M~40M:\n")
    for match in matches:
        hybrid_out.write(match + " " + str(matches[match][0]))
        try:
            hybrid_out.write(" " + str(matches[match][1]))
        except IndexError:
            pass
        hybrid_out.write("\n")

    hybrid_out.close()
