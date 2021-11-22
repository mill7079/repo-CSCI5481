from tree_sankoff import *

if __name__ == '__main__':
    sequences = read_alignment('seqs/clustal.aln')  # parse alignment file into dict - fasta_name:seq
    tree = Phylo.read("tree.dnd", "newick")  # create tree from alignment results (tree from Clustal Omega)
    # sequences = read_alignment('toy_align.txt')
    # tree = Phylo.read('toy.dnd', 'newick')

    # run sankoff
    trees = create_trees(tree, sequences)
    parsimony_score = 0
    for tree in trees:
        new_score = sankoff(tree, trees[tree])
        if new_score < inf:
            parsimony_score += new_score
            print(parsimony_score, new_score)
        trace_back(tree.clade, trees[tree])

    pars_tree = sum_trees(trees)
    print("final tree:")
    print(pars_tree)
    # Phylo.draw(pars_tree)
    print("score:")
    print(parsimony_score)
