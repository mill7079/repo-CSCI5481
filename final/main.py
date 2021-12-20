# chaos time :)
import argparse


# calculate parsimony score between two individuals
def parsimony(s1, s2):
    score = 0
    for l in s2.upper():
        if l not in s1.upper():
            score += 1
    return score


def parse_file(file, num):
    chrm = open(file)
    headers = chrm.readline()
    data = dict()
    for l in chrm:
        line = l.split(" ")
        if line[0] not in data:
            data[line[0]] = []

        for i in range(11, num + 11):
            data[line[0]].append(line[i])

    return data


if __name__ == '__main__':
    # print(parsimony("tt", "tt"))
    # print(parsimony("tt", "ta"))
    # print(parsimony("tt", "aa"))
    # print(parsimony("at", "ta"))

    # cool shit
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--num_inds", default=10, type=int)  # individuals per data file
    parser.add_argument("-f", "--files", nargs='+')  # data files
    args = parser.parse_args()
    # print(args.num_inds)
    # print(args.files)

    snps = None
    for file in args.files:
        file_data = parse_file(file, args.num_inds)

        # combine data from all files
        if snps is None:
            snps = file_data
        else:
            for snp in file_data:
                if snp in snps:
                    snps[snp].extend(file_data[snp])
                else:
                    snps[snp] = file_data[snp]

    print(snps)




