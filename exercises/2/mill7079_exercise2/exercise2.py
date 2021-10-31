def read(genome, start, direction):
    stops = ['TGA', 'TAG', 'TAA']  # start is ATG
    orfs = []
    frame = ''

    if direction < 0:
        genome = genome[::-1]

    start_index = 0

    for i in range(start, len(genome), 3):
        codon = genome[i:i+3]
        if codon == 'ATG' and frame == '':
            frame = codon
            start_index = i
        elif codon in stops and frame != '':
            frame += codon
            orfs.append((frame, direction, start_index, i+2))
            frame = ''
        elif frame != '':
            frame += codon

    return orfs


def print_orf(orf):
    print('(' + str(orf[2]), ':', str(orf[3]) + ')', 'Strand:', '+' if orf[1] > 0 else '-', 'Length:', len(orf[0]))


if __name__ == '__main__':
    file = open("GCA_002859165.1_ASM285916v1_genomic.fna")
    file.readline()
    genome = file.read().replace('\n', '')
    orfs = []
    for i in range(0, 6):  # handles both strands
        arr = read(genome, i % 3, (i%2)*2-1)
        for orf in arr:
            orfs.append(orf)

    print(len(genome))
    print(len(genome)/3)

    # print("ORFs longer than 60: ")
    # for orf in orfs:
    #     if (len(orf[0])/3) >= 60:
    #         print(orf, len(orf[0]), len(orf[0])/3)
    #
    # print('\n****************************')
    # print("ORFs longer than 100: ")
    # for orf in orfs:
    #     if (len(orf[0])/3) >= 100:
    #         print(orf, len(orf[0]), len(orf[0])/3)
    #
    # print('\n****************************')
    # print("ORFs longer than 200: ")
    # for orf in orfs:
    #     if (len(orf[0])/3) >= 200:
    #         print(orf, len(orf[0]), len(orf[0])/3)

    print("ORFs longer than 60: ")
    for orf in orfs:
        if len(orf[0]) >= 60:
            print_orf(orf)

    print('\n****************************')
    print("ORFs longer than 100: ")
    for orf in orfs:
        if len(orf[0]) >= 100:
            print_orf(orf)

    print('\n****************************')
    print("ORFs longer than 300: ")
    for orf in orfs:
        if len(orf[0]) >= 300:
            print_orf(orf)
            if orf[1] < 0:
                print("reversed indices:")
                print(len(genome) - orf[2], len(genome) - orf[3])

    count = 0
    print('\n****************************')
    print("ORFs longer than 600: ")
    for orf in orfs:
        if len(orf[0]) >= 600:
            count += 1
            print(str(count) + ': ', end='')
            print_orf(orf)

    #find counts
    counts = {30: 0, 75: 0, 150: 0, 300: 0, 600: 0}
    for orf in orfs:
        if len(orf[0]) >= 600:
            counts[600] += 1
        if len(orf[0]) >= 300:
            counts[300] += 1
        if len(orf[0]) >= 150:
            counts[150] += 1
        if len(orf[0]) >= 75:
            counts[75] += 1
        if len(orf[0]) >= 30:
            counts[30] += 1

    for count in counts:
        print(count, ':', counts[count])
