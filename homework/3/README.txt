Allison Miller
mill7079
Homework 3

Note: The included .zip file contains all the submission files, in case that is
easier to grade.

To run:
Problem 2:
The toy example is coded into the main method. To run:
    python3 sankoff.py


Files included for each problem:

Problem 1:
    - msa.pdf - multiple sequence alignment results
    - tree.txt - pairwise distances with tree structure
    - tree.pdf - image file of tree

Problem 2:
    - sankoff.py - implementation of sankoff algorithm for toy examples, with main method to run on given example
    - toy_output.pdf - results from toy example

Problem 3:
    - not successfully completed
    - tree_sankoff.py - attempt at re-doing sankoff algorithm with BioPython trees

Issues/Notes:

It took me a rather embarrassingly long time to figure out how to use Sankoff on sequences and not just single characters; I couldn't
    figure out how the BioPython trees were relevant to the toy examples, so I ended up trying to create my own tree structure...which
    works, but only for trees with a number of leaves that is an exact power of 2 (aka, the toy example). I realized way too late
    that it wouldn't translate well to using actual sequences and complicated pre-existing tree structures, and as such, had to
    completely re-implement the algorithm for Problem 3.

For problem 3, I attempted to figure out how to use the BioPython trees, but the documentation is all over the place, and I couldn't
    figure out how to either traverse the tree effectively or store the information I needed at each node; hence, problem 3 currently
    isn't implemented (though you can see the attempt in tree_sankoff.py).

For problem 2, my code works for the toy example, though it does appear to make different choices than the provided solution when it
    comes across ties for lowest score (i.e. for the two intermediate results, indices 3 and 5 (from 0) of the left-hand sequence are
    different from the given solution, and indices 1 and 2 have the same issue with the right-hand sequence. I calculated the values by
    hand, and the algorithm appears to be working correctly, just with different tie-breaker choices (i.e. for position 3 in the left
    intermediate sequence, the solution lists C, while my solution lists G; when calculating the scores by hand, both G and C contributed
    the minimum score of 4 at that stage, so either of those characters should work).

Additionally, for problem 1, I had trouble figuring out how to get just the pairwise distances, as the flags on the tool I was using
    (the downloaded version of Clustal Omega) were apparently deprecated and thus ignored when I tried running the program with them;
    however, I think the distances are included with the tree structure, so I've included that file in place of the pairwise distances.

Sorry about the mess...