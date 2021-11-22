Allison Miller
mill7079
Homework 3

Note: The included .zip file contains all the submission files, in case that is
easier to grade.

To run:
Problem 2:
python3 tree_sankoff.py
- automatically runs toy example, using toy.dnd as the tree file and toy_align as the alignment file

Problem 3:
python3 p_main.py
- automatically runs actual example, using tree.dnd as the tree file and clustal.aln as the alignment file


Files included for each problem:

Problem 1:
    - msa.pdf - multiple sequence alignment results
    - tree.txt - pairwise distances with tree structure
    - tree.pdf - image file of tree

Problem 2:
    - tree_sankoff.py - implementation of sankoff algorithm using BioPython trees, with main method to run as listed above
        - (submission items 1 and 2 in the writeup for this problem are the same file)
    - toy_example.pdf - results from toy example
    - sankoff.py - old implementation of sankoff algorithm; required for Node class access
    - toy.dnd - formatted tree file of toy example for BioPython
    - toy_align.txt - alignment file for toy example

Problem 3:
    - tree_sankoff.py - implementation of sankoff algorithm using BioPython trees (also submitted for problem 2)
    - p_main.py - main method to run actual example
    - tree.dnd - formatted tree file of alignment for BioPython
    - clustal.aln - multiple sequence alignment results, needed for sankoff
    - results.pdf - total parsimony score and inferred most likely sequence of internal nodes


Issues/Notes:

I had a hard time understanding how to use Sankoff on sequences and not just single characters for a while; as a result, the first
    attempt (sankoff.py) is limited to mostly running the toy example, and was unsuitable for problem 3. I rewrote the algorithm and
    the surrounding support code to use BioPython trees (tree_sankoff.py) but it took a long time to figure out how to use the trees,
    and I was unable to complete this portion in time for the deadline. I believe the solution is mostly correct; however, there is
    at least one tree somewhere in the mix that returned infinity as the minimum score for the root, causing the parsimony score to be
    infinity, so something's not quite right. I have no idea how to debug the issue and it's late enough as it is, so hopefully it's not
    screwing too many things up - I just ignored that value for the final score computation.

For problem 2, my intermediate results don't quite match the given solution in the writeup for the toy example despite the final result
    matching perfectly (i.e. for the two intermediate results, indices 3 and 5 (from 0) of the left-hand sequence are
    different from the given solution, and indices 1 and 2 have the same issue with the right-hand sequence); I calculated some of the
    values by hand, and it appears to just be an issue with which base is chosen in the event of tied minimum scores. For example, for
    position 3 in the left intermediate sequence, the solution lists C, while my output lists G; when calculating the scores by hand,
    both G and C contributed the minimum score of 4 at that stage, so either of those characters should work for the analysis despite
    not matching the solution.

Additionally, for problem 1, I had trouble figuring out how to get just the pairwise distances, as the flags on the tool I was using
    (the downloaded version of Clustal Omega) were apparently deprecated and thus ignored when I tried running the program with them;
    however, I think the distances are included with the tree structure, so I've included that file in place of the pairwise distances.
    I also wasn't entirely sure what the submission instructions for problem 3 were asking for, so hopefully the results file is correct.