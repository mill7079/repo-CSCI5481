from math import inf

alph = ['A', 'U', 'G', 'C', '-']

# indices correspond with alph for easy matching
scoring_matrix = [[0, 3, 4, 9, 8],
                  [3, 0, 2, 4, 8],
                  [4, 2, 0, 4, 8],
                  [9, 4, 4, 0, 8],
                  [8, 8, 8, 8, 8]]


class Node:
    def __init__(self, leaf_char=None, left=None, right=None):
        self.scores = [inf, inf, inf, inf, inf, (inf, inf), (inf, inf), (inf, inf), (inf, inf), (inf, inf)]
        if leaf_char is not None:
            self.scores[alph.index(leaf_char.upper())] = 0
            self.char = leaf_char

        self.left = left
        self.right = right
        self.char = 'z'

    def set_score(self, char, score):
        self.scores[alph.index(char.upper())] = score

    def set_min_char(self):
        self.char = alph[self.scores.index(min(self.scores))]

    def __str__(self):
        return self.char + ' ' + str(self.scores) + '\n'

    def __repr__(self):
        return self.char + ' ' + str(self.scores) + '\n'


def length_check(arr):
    if len(arr) < 1:
        print("empty array")
        return False

    val = len(arr[0])
    for seq in arr:
        if len(seq) != val:
            return False
    return True


def dist(a, b):
    return scoring_matrix[alph.index(a)][alph.index(b)]


# creates trees consisting only of leaves to start with
def create_trees(seqs):
    trees = []
    for i in range(0, len(seqs[0])):
        tree = []
        for j in range(0, len(seqs)):
            tree.append(Node(leaf_char=seqs[j][i]))
        trees.append(tree)

    # print(trees)
    return trees


def trace_back(root):
    root.char = alph[root.scores.index(min(root.scores[0:len(alph)]))]  # set min char on root

    if root.left is not None:
        trace_back(root.left)
    if root.right is not None:
        trace_back(root.right)


# performs sankoff algorithm on given list of sequences
def sankoff(seqs):
    if not length_check(seqs):
        print("bad lengths")
        return None

    # each tree deals with one position from the sequence
    # so if len(trees) = 2, len(seq) = 2, 1st tree is for 1st letter, 2nd for 2nd
    trees = create_trees(seqs)  # create trees from leaves

    # for each tree, perform alg
    for tree in trees:
        # just for now - join each pair of seqs for next step up
        for i in range(0, len(tree), 2):
            if i + 1 >= len(tree):
                print("exit loop")
                break

            parent = Node(left=tree[i], right=tree[i+1])
            for j in range(0, len(alph)):  # for each parent letter
                min_left = inf
                min_right = inf
                trace = [0, 0]
                for x in range(0, len(alph)):  # for each character in alph in leaf - finding mins
                    left_score = tree[i].scores[x] + dist(alph[j], alph[x])
                    right_score = tree[i+1].scores[x] + dist(alph[j], alph[x])
                    if left_score < min_left:
                        min_left = left_score
                        trace[0] = tree[i].scores[x]

                    if right_score < min_right:
                        min_right = right_score
                        trace[1] = tree[i+1].scores[x]
                    # min_left = min(min_left, tree[0].scores[x] + dist(alph[j], alph[x]))
                    # min_right = min(min_right, tree[1].scores[x] + dist(alph[j], alph[x]))

                parent.scores[j] = min_left + min_right
                parent.scores[j + len(alph)] = (trace[0], trace[1])
            tree.append(parent)

    # TODO trace back down each tree
        # set min char on root
        # use stored trace values and left/right pointers to set chars on leaves
    for tree in trees:
        trace_back(tree[-1])
        # root = tree[-1]
        # root.char = alph[root.scores.index(min(root[0:len(alph)]))]  # set min char on root
        #
        # left = root.left
        # right = root.right
        # parent = root
        #
        # while left is not None and right is not None:  # set chars for other nodes

    # TODO combine all trees into single tree

    final_tree = [''] * len(trees[0])
    for tree in trees:
        for i in range(0, len(tree)):
            # final_tree[i] += tree[i]
            print(tree[i])
            final_tree[i] += tree[i].char
    return final_tree


if __name__ == '__main__':
    print("main")
    # sequences = ['AU', 'GU', 'CG', 'CA']
    sequences = ['AUUCGUGAUU', 'AUUGAA-AUU', 'GUCCUCGGUU', 'GA-CACGAUC']

    # sankoff(sequences)
    print(sankoff(sequences))
