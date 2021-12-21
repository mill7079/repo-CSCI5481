# chaos time :)
import argparse
import copy
from node import Node, alph
from math import inf, log


# calculate parsimony score between two individuals
# formula given by professor in assignment video - count number of mismatches
def parsimony(s1, s2):
    score = 0
    for l in s2.upper():
        if l not in s1.upper():
            score += 1
    return score


# Jukes-Cantor model based distance
def dist(n1, n2):
    # print(parsimony(n1.seq, n2.seq))
    p = parsimony(n1.seq, n2.seq) / len(n1.seq)
    # if p == 0.0:
    #     return 0
    # print(1-p)
    # return -log(1 - p)

    return (-19/20) * log(1 - ((19/20) * p))


# parses a single file
# returns a dictionary: key = SNP name, value = array of individual values
def parse_file(file, num):
    try:  # can't parse anything if there's no file to open
        chrm = open(file)
    except FileNotFoundError:
        print("file", file, "not found.")
        return None

    headers = chrm.readline()  # not used but still need to read the line anyway to get it out of the way
    data = dict()

    # each line is a new SNP
    for l in chrm:
        line = l.split(" ")  # the files look space-separated so this should work
        if line[0] not in data:  # if new SNP, add to dictionary
            data[line[0]] = []

        for i in range(11, num + 11):  # append individual values for later processing
            if i < len(line):  # avoids errors if user requested more values than the file contains
                data[line[0]].append(line[i])
                if line[i][0] != line[i][1]:  # testing if a file contains a mismatched pair - none of the 11 files did
                    print(line[i])
            else:  # if i >= len(line)
                break

    chrm.close()
    return data


# combines data from files into one dictionary for tree building
# args = list of files passed in through command line
def combine_files(args):
    data = None
    for file in args.files:
        file_data = parse_file(file, args.num_inds)  # create name/values dictionary for each file
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
            # formula from Wikipedia because it's easier to search than the textbook:
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
    all_nodes = []  # start list of all nodes in tree; returned for later use in Sankoff
    count = 0  # number of total nodes; used for labeling internal nodes for ease of results checking
    for node in tree:
        all_nodes.append(node)

    # find all unpaired nodes (avoid using central node for pairing because it doesn't have distance yet)
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

        # find shared node between two nodes - should only be one, i think, otherwise there's a cycle
        for node in n1.connections:
            if node in n2.connections:
                parent = node
                break

        # join taxa to new node (disconnect and reconnect)
        # redo connections
        # new_parent = Node()
        new_parent = Node(str(count))
        count += 1

        # disconnect from old parent, reconnect to new parent (empty node), connect that new node to old parent
        n1.disconnect(parent)
        n2.disconnect(parent)
        new_parent.connect(n1)
        new_parent.connect(n2)
        parent.connect(new_parent)

        # remove old nodes from list of unpaired, add new node
        unpaired.remove(n1)
        unpaired.remove(n2)
        unpaired.append(new_parent)
        all_nodes.append(new_parent)

        # set parent/children relationship for use in sankoff/nni
        n1.parent = new_parent
        n2.parent = new_parent
        new_parent.children = [n1, n2]

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
        # first add new values (both full row and in other rows)...
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
    return all_nodes


# sankoff scoring function
def score(parent, children):
    c1 = children[0]
    c2 = children[1]
    for i in range(0, len(parent.scores)):  # for each letter in parent...
        min_c1 = inf
        min_c2 = inf
        for j in range(0, len(alph)):  # for each letter in child, find minimum score per child
            c1_score = c1.scores[j] + parsimony(alph[i], alph[j])
            c2_score = c2.scores[j] + parsimony(alph[i], alph[j])

            if c1_score < min_c1:
                min_c1 = c1_score
            if c2_score < min_c2:
                min_c2 = c2_score

        # print(min_c1, min_c2)
        parent.scores[i] = min_c1 + min_c2
        parent.scored = True


# returns number of unscored nodes
# for starting sankoff - mostly just finding leaves
def unscored_nodes(nodes):
    count = 0
    for node in nodes:
        if not node.scored:
            count += 1
    return count


# sankoff implementation
# same iffy implementation as HW3...
def sankoff(nodes):
    nodes_left = unscored_nodes(nodes)
    while nodes_left > 0:  # while some nodes have not yet been scored
        for node in nodes:
            if not node.scored:  # if parent has not been computed...
                scored = True
                for child in node.children:  # ...check if children have been
                    scored = scored and child.scored

                if scored:  # if both are true then score the parent using the children
                    score(node, node.children)
                    nodes_left -= 1
                    if nodes_left == 0:  # scored last node, so return parsimony score
                        return min(node.scores[0:len(alph)])


# nearest neighbor interchange implementation
# nodes is just a list of all nodes in the tree in no particular order
# so this is likely not super efficient
# num is the -i parameter, only used for printing the trees
# copies are used to avoid screwing with other data, though this might be overboard
def nni(nodes, num):
    nodes_copy = copy.deepcopy(nodes)
    # for node in nodes_copy:
    for i in range(0, len(nodes_copy)):
        node = nodes_copy[i]
        num_children = 0
        for child in node.children:  # can exchange subtrees if node has 4 grandchildren - overly simple case
            num_children += len(child.children)

        # not really sure how to do this iteratively..
        if num_children == 4:
            # create a copy of the tree for each possibility
            trees = [copy.deepcopy(nodes), copy.deepcopy(nodes), copy.deepcopy(nodes)]
            sankoffs = [-1, -1, -1]

            # original layout
            sankoffs[0] = sankoff(trees[0])
            # print("Tree 0:\n", print_tree(trees[0][num], ""))

            # first alternate layout - exchange right subtrees
            node = trees[1][i]
            c1 = node.children[0]
            c2 = node.children[1]
            temp = c1.children[1]
            c1.children[1] = c2.children[1]
            c2.children[1] = temp

            for j in range(0, len(c1.children)):
                c1.children[j].parent = c1
                c2.children[j].parent = c2
            sankoffs[1] = sankoff(trees[1])
            # print("Tree 1:\n", print_tree(trees[1][num], ""))
            # print(trees[1])

            # second alternate layout - exchange right of one with left of the other
            node = trees[2][i]
            c1 = node.children[0]
            c2 = node.children[1]
            temp = c1.children[1]
            c1.children[1] = c2.children[0]
            c2.children[0] = temp

            for j in range(0, len(c1.children)):
                c1.children[j].parent = c1
                c2.children[j].parent = c2
            sankoffs[2] = sankoff(trees[2])
            # print("Tree 2:\n", print_tree(trees[2][num], ""))

            # print(sankoffs)
            # use the optimal tree for the rest of the algorithm
            # though in practice this...doesn't really affect anything
            nodes_copy = trees[sankoffs.index(min(sankoffs))]

    # return sankoff(nodes_copy)


# print the tree - debugging
# visual representation of tree - shows which nodes are paired, though still kind of hard to understand entire structure
def print_tree(node, tabs=""):
    node.visited = True
    tabs += '---'
    acc = ''
    for c in node.children:
        if not c.visited:
            acc += '\n' + tabs + print_tree(c, tabs)
            c.visited = True
    return str(node) + acc


if __name__ == '__main__':
    # parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--num_inds", default=10, type=int)  # individuals per data file
    parser.add_argument("-f", "--files", nargs='+')  # data files
    parser.add_argument("-t", "--num_trees", default=1, type=int)  # number of trees to see in output
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
    full_trees = dict()
    full_copies = dict()
    count = 0
    print("joining trees")
    for tree in trees:
        full_trees[tree] = join(trees[tree])
        trees[tree][-1].children = trees[tree][-1].connections
        full_copies[tree] = copy.deepcopy(full_trees[tree])
        count += 1
        print(count, "out of", len(trees), "joined")

    # Sankoff algorithm
    print('\n\nSankoff/Neighbor Joining results:\n')
    i = 0
    for tree in trees:
        # print(len(full_trees[tree]))
        print("Parsimony score:", sankoff(full_trees[tree]))
        print("Tree:\n" + print_tree(trees[tree][-1], ""), '\n')
        # print("Copy:\n" + print_tree(full_copies[tree][4], ""))
        i += 1
        if i >= args.num_trees:
            break

    # Nearest Neighbor Interchange
    print('\n\nNearest Neighbor Interchange results:\n')
    i = 0
    for tree in trees:
        # run nearest neighbor interchange on tree copies
        # print('\n\n', print_tree(trees[tree][-1], ""))
        nni(full_copies[tree], args.num_inds)
        print("NNI Parsimony score:", sankoff(full_copies[tree]))
        print("NNI Tree:\n" + print_tree(full_copies[tree][args.num_inds*len(args.files)], ""), '\n')
        i += 1
        if i >= args.num_trees:
            break
