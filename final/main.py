# chaos time :)
import argparse
from node import Node
from math import inf


# calculate parsimony score between two individuals
# formula given by professor in assignment video - count number of mismatches
def parsimony(s1, s2):
    score = 0
    for l in s2.upper():
        if l not in s1.upper():
            score += 1
    return score


# method wrapper in case I figure out a different distance method
# but I don't know what else I'd use so...it's just passing the args to parsimony
# takes in two nodes, passes sequences to parsimony
def dist(n1, n2):
    return parsimony(n1.seq, n2.seq)


# parses a single file
# returns a dictionary: key = SNP name, value = array of individual values
def parse_file(file, num):
    try:  # can't parse anything if there's no file to open
        chrm = open(file)
    except FileNotFoundError:
        print("file", file, "not found.")
        return None

    headers = chrm.readline()  # not used but still need to read the file anyway
    data = dict()

    for l in chrm:
        line = l.split(" ")  # the files look space-separated so this should work
        if line[0] not in data:  # if new SNP, add to dictionary
            data[line[0]] = []

        for i in range(11, num + 11):  # append individual values for later processing
            data[line[0]].append(line[i])

    return data


# combines data from files into one dictionary for tree building
# args = list of files passed in through command line
def combine_files(args):
    data = None
    for file in args.files:
        file_data = parse_file(file, args.num_inds)
        if file_data is None:  # skip file if was not found in parse_file
            continue

        # combine data from all files
        if data is None:  # first file, nothing to merge
            data = file_data
        else:  # merge lists
            for snp in file_data:
                if snp in data:  # if the SNP name is present in both files...
                    data[snp].extend(file_data[snp])  # ...then append the individuals from the new population onto the list
                else:  # new SNP
                    data[snp] = file_data[snp]

    return data


# create starting tree for each snp
# star shape - all leaf nodes (individuals' values) connected to central node
def create_star_tree(snp):
    tree = []

    # create leaves
    for seq in snp:
        tree.append(Node(seq))

    # connect leaves to center node
    center = Node()
    for node in tree:
        center.connect(node)

    # center node is always last in list
    tree.append(center)
    return tree


# finds nxn Q matrix used in neighbor joining algorithm
# matrix is a 2D dictionary: keys = nodes, values = either list of distances for that node or distance value
def find_q_matrix(unpaired, D):
    # for n in D:
    #     print(n, D[n])

    q = dict()
    # print(unpaired)
    for node in unpaired:
        q[node] = dict()
        for node2 in unpaired:
            # Q(i,j) = (n-2)*D(i,j) - sum(d(i,k)) - sum(d(j,k))
            q[node][node2] = (len(unpaired) - 2) * D[node][node2]
            sum_1 = 0
            sum_2 = 0
            for nodek in unpaired:
                sum_1 += D[node][nodek]
                sum_2 += D[node2][nodek]
            q[node][node2] -= sum_1
            q[node][node2] -= sum_2

    return q


# perform neighbor joining algorithm
# tree is an array of nodes, initially connected in star tree pattern
def join(tree):
    # find all unpaired nodes (avoid using central node because it doesn't have distance yet)
    unpaired = []
    for node in tree:
        if node.seq != "":
            unpaired.append(node)

    # find initial distance matrix - 2D dictionary like Q matrix with nodes as keys
    D = dict()
    for node in unpaired:
        D[node] = dict()
        for node2 in unpaired:
            D[node][node2] = dist(node, node2)  # should be 0 if equal

    # neighbor joining algorithm - iterate through steps until tree structure is complete
    # number of iterations = n-3 where n is number of (leaf) nodes
    for _ in range(0, len(tree) - 3):
        # find Q matrix
        Q = find_q_matrix(unpaired, D)

        # find pair of distinct taxa for which Q(i,j) is min
        # mi = [0, 0]  # min index tracker
        # for x in range(0, len(Q)):
        #     for y in range(0, len(Q[x])):
        #         if x != y and Q[x][y] < Q[mi[0]][mi[1]]:
        #             mi[0] = x
        #             mi[1] = y
        m_nodes = [unpaired[0], unpaired[1]]
        for node in unpaired:
            for node2 in unpaired:
                if node is not node2 and Q[node][node2] < Q[m_nodes[0]][m_nodes[1]]:
                    m_nodes = [node, node2]

        # n1 = tree[mi[0]]
        # n2 = tree[mi[1]]
        # set references for ease of typing
        n1 = m_nodes[0]
        n2 = m_nodes[1]
        parent = None

        # find shared node between two nodes - should only be one, i think
        for node in n1.connections:
            if node in n2.connections:
                parent = node
                break

        # join taxa to new node (disconnect and reconnect)
        # redo connections
        new_parent = Node()

        # disconnect from old parent, reconnect to new parent (empty node), connect that new node to old parent
        n1.disconnect(parent)
        n2.disconnect(parent)
        new_parent.connect(n1)
        new_parent.connect(n2)
        parent.connect(new_parent)

        # remove nodes from list of unpaired, add new node
        unpaired.remove(n1)
        unpaired.remove(n2)
        unpaired.append(new_parent)

        # find distance from each of paired taxa to new node - do I actually need these...?
        # d_n1 = 0.5 * D[mi[0]][mi[1]]
        # d_n1 = 0.5 * D[n1][n2]
        # sum_f = 0
        # sum_g = 0
        # # for k in range(0, len(tree)):
        # for node in unpaired:
        #     # sum_f += D[mi[0]][k]
        #     sum_f += D[n1][node]
        #     # sum_g += D[mi[1]][k]
        #     sum_g += D[n2][node]
        #
        # d_n1 += (1 / (2 * (len(tree) - 2))) * (sum_f - sum_g)
        # # d_n2 = D[mi[0]][mi[1]] - d_n1
        # d_n2 = D[n1][n2] - d_n1

        # find distance to new node from all other taxa - redo distance matrix
        # first add new values...
        D[new_parent] = dict()
        for node in unpaired:
            if node != new_parent:
                D[new_parent][node] = 0.5 * (D[n1][node] + D[n2][node] - D[n1][n2])
                D[node][new_parent] = D[new_parent][node]

        D[new_parent][new_parent] = 0

        # ...then remove old values from matrix
        # have to remove both row for each node and the nodes' entries from the other rows
        # since I have the entire matrix defined, not just above the diagonal
        D.pop(n1)
        D.pop(n2)
        for row in D:
            if n1 in D[row]:
                D[row].pop(n1)
            if n2 in D[row]:
                D[row].pop(n2)

    # print("tree created. ", len(unpaired))


# print the tree - debugging
# visual representation of tree - shows which nodes are paired, though still kind of hard to understand entire structure
def print_tree(node, tabs):
    node.visited = True
    # tabs += '\t'
    tabs += '---'
    acc = ''
    for c in node.connections:
        if not c.visited:
            acc += '\n' + tabs + print_tree(c, tabs)
            c.visited = True

    return str(node) + acc


if __name__ == '__main__':
    # print(parsimony("tt", "tt"))
    # print(parsimony("tt", "ta"))
    # print(parsimony("tt", "aa"))
    # print(parsimony("at", "ta"))

    # cool shit
    # parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--num_inds", default=10, type=int)  # individuals per data file
    parser.add_argument("-f", "--files", nargs='+')  # data files
    args = parser.parse_args()
    # print(args.num_inds)
    # print(args.files)

    # process file data
    snps = combine_files(args)
    # print(snps)

    # create trees out of data
    trees = dict()
    for snp in snps:
        # trees[snp] = create_leaves(snps[snp])
        trees[snp] = create_star_tree(snps[snp])
    # print(trees)

    # use neighbor joining algorithm to join the trees (in place)
    for tree in trees:
        join(trees[tree])

    # # debugging
    # for tree in trees:
    #     print(print_tree(trees[tree][-1], ""))
    #     break



