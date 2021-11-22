from Bio import Phylo
import copy
from sankoff import Node, dist, alph
from math import inf
#
# sequences = []
# trees = dict()
#
#
# def read_seqs(align_file):
#     # files = glob.glob("seqs/*.fasta")
#     # for file in files:
#     #     f = open(file)
#     #     line = f.readline()
#     #     sequences.append((line.split(' ')[0], f.read().replace('\n', '')))
#
#
#
# if __name__ == '__main__':
#     tree = Phylo.read("tree.dnd", "newick")
#     print(tree)
#     # print(tree.find_clades("silva|CR378665|83603|85179|"))
#     # Phylo.draw(tree)
#     read_seqs()
#     seq = sequences[0][1]
#     for i in range(0, len(seq)):
#         new_tree = copy.deepcopy(tree)
#         leaves = new_tree.get_terminals()
#         tree_dict = dict()
#         for j in range(0, len(leaves)):
#             if len(seq) != len(sequences[j][1]):
#                 print(len(seq), len(sequences[j][1]))
#             tree_dict[leaves[j]] = Node(leaf_char=sequences[j][1][i])
#         trees[new_tree] = tree_dict
#     # for seq in sequences:
#     #     new_tree = copy.deepcopy(tree)
#     #     print(tree.find_clades(seq[1:]))
#
#     # for i in range(0, len(tree.get_terminals())):
#     #     new_tree = copy.deepcopy(tree)
#     #
#     #     trees.append(new_tree)
#


# read alignment file to concatenate aligned sequences if needed
# assumes structure like the one returned from clustal omega
def read_alignment(align_file):
    sequences = dict()

    file = open(align_file)
    file.readline()
    file.readline()
    file.readline()
    for line in file:
        split = line.split(' ')
        if split[0] != '' and split[0] != '\n':
            if split[0] not in sequences:
                sequences[split[0]] = split[-1].strip()
            else:
                sequences[split[0]] += split[-1].strip()

    file.close()
    return sequences


# creates the individual trees for each sequence position
def create_trees(init_tree, sequences):
    trees = dict()
    seq_len = 0
    for seq in sequences:
        seq_len = len(sequences[seq])
        break

    for i in range(0, seq_len):
        new_tree = copy.deepcopy(init_tree)
        trees[new_tree] = create_nodes(new_tree, sequences, i)

    return trees


# creates the dictionary of nodes for each tree
def create_nodes(tree, sequences, position):
    leaves = tree.get_terminals()
    others = tree.get_nonterminals()
    nodes = dict()

    for node in others:
        nodes[node] = Node()

    for leaf in leaves:
        nodes[leaf] = Node(leaf_char=sequences[leaf.name][position])

    return nodes


# taken from https://biopython.org/wiki/Phylo_cookbook
# unused
# def get_parent(tree, child_clade):
#     node_path = tree.get_path(child_clade)
#     return node_path[-2]


# mostly taken from sankoff.py
# finds score matrix of parent given its children
def score(parent, left, right):
    for i in range(0, len(alph)):  # for each letter in parent
        min_left = inf
        min_right = inf
        trace = [0, 0]
        for j in range(0, len(alph)):  # for each letter in leaf  --  finding mins
            left_score = left.scores[j] + dist(alph[i], alph[j])
            right_score = right.scores[j] + dist(alph[i], alph[j])

            if left_score < min_left:
                min_left = left_score
                trace[0] = left.scores[j]
            if right_score < min_right:
                min_right = right_score
                trace[1] = right.scores[j]

        # set score matrix of parent node and traceback refs
        parent.scores[i] = min_left + min_right
        parent.scores[i + len(alph)] = (trace[0], trace[1])
        parent.char = 'checked'


# yikes
# super inefficient implementation of sankoff
# couldn't figure out how to do the tree traversal any better
def sankoff(tree, nodes):
    nodes_left = len(nodes) - len(tree.get_terminals())
    while nodes_left > 0:
        # print(nodes_left)
        # for node in nodes:
        #     n1 = nodes[node]
        #     for node2 in nodes:
        #         n2 = nodes[node2]
        #         if n1.char != 'z' and n2.char != 'z':  # if nodes are leaves or already computed
        #             if get_parent(tree, node) == get_parent(tree, node2):  # if nodes share same parent
        #                 # compute score for that parent
        for node in nodes:
            parent = nodes[node]
            if parent.char == 'z':  # if parent has not been computed
                scored = True
                children = []
                # for child in parent.clade:  # if children have been computed
                for child in node:
                    # scored &= child.char != 'z'
                    scored &= nodes[child].char != 'z'
                    # children.append(child)
                    children.append(nodes[child])

                if scored:
                    # score(parent, nodes[children[0]], nodes[children[1]])
                    score(parent, children[0], children[1])
                    nodes_left -= 1
                    if nodes_left == 0:  # just scored parent - find parsimony score
                        # parsimony_score += min(parent.scores[0:len(alph)])
                        return min(parent.scores[0:len(alph)])


# trace scores back, setting name of clade to char for now
def trace_back(root, nodes):
    # root.char = alph[root.scores.index(min(root.scores[0:len(alph)]))]  # set min char on root
    root_node = nodes[root]
    root_node.char = alph[root_node.scores.index(min(root_node.scores[0:len(alph)]))]
    root.name = root_node.char

    for child in root:
        trace_back(child, nodes)


# merges chars of two trees
def merge(tree1, tree2):
    tree1.name += tree2.name

    children1 = []
    for child in tree1:
        children1.append(child)

    children2 = []
    for child in tree2:
        children2.append(child)

    for i in range(0, len(children1)):
        merge(children1[i], children2[i])


# sums all trees to get final sequences
def sum_trees(trees):
    final_tree = None
    for tree in trees:
        if final_tree is None:
            final_tree = copy.deepcopy(tree)
        else:
            merge(final_tree.clade, tree.clade)

    return final_tree


if __name__ == '__main__':
    # sequences = read_alignment('seqs/clustal.aln')  # parse alignment file into dict - fasta_name:seq
    # tree = Phylo.read("tree.dnd", "newick")  # create tree from alignment results (tree from Clustal Omega)
    sequences = read_alignment('toy_align.txt')
    tree = Phylo.read('toy.dnd', 'newick')
    # print(tree)

    # run sankoff
    trees = create_trees(tree, sequences)
    parsimony_score = 0
    for tree in trees:
        parsimony_score += sankoff(tree, trees[tree])
        # trace_back(tree, trees[tree])
        trace_back(tree.clade, trees[tree])
        # print(tree)
        # break

    pars_tree = sum_trees(trees)
    print("final tree:")
    print(pars_tree)
    # Phylo.draw(pars_tree)
    print("score:")
    print(parsimony_score)
