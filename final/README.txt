Allison Miller
Final Project - Option 3 (Parsimony Analysis)
5278509

To run:
The program includes three command-line flags, one for the number of individuals to read from each data file (-i or --num_inds),
    one to indicate the start of the list of files (-f or --files), and one for the number of trees to be printed in the output.
    The flags can be included in any order (though the files should all be listed together). Examples:

To run an analysis on two files with 20 individuals from each file and see 3 trees as output:
    python3 main.py -i 20 -t 3 -f data/genotypes_chrM_ASW_phase3.2_consensus.b36_fwd.txt data/genotypes_chrM_JPT_phase3.2_consensus.b36_fwd.txt

To run an analysis on three files with 50 individuals from each file (with the default one result tree):
    python3 main.py --num_inds 50 --files data/genotypes_chrM_ASW_phase3.2_consensus.b36_fwd.txt data/genotypes_chrM_JPT_phase3.2_consensus.b36_fwd.txt data/genotypes_chrM_CEU_phase3.2_consensus.b36_fwd.txt

Each file should be a path to the file (i.e. I had the files all in a directory named 'data' while testing, hence the paths). Number
    of individuals should be an integer.

Most code is contained in main.py (i.e. functions for joining the trees and performing Sankoff and such); node.py just contains
    the definition for the Node class used to create the trees, as well as the array of possible values.


Output:

The output consists of three parts: progress printing while the trees are being created using the neighbor joining algorithm, so you
    know about how much longer the code will run; the normal parsimony scores and trees; and parsimony scores and trees after using
    the Nearest Neighbor Interchange algorithm. The number of score/tree pairs printed for each section depends on the command line
    flag -t (defaults to 1).


Testing ended up being severely restricted due to ongoing battery issues with my computer. (Additionally, I wasn't entirely sure
    what we were supposed to write in the report, so hopefully it's not too far off...)


Resources:
https://en.wikipedia.org/wiki/Neighbor_joining#The_algorithm
    - Description and example walkthrough (with pictures!) of the neighbor joining algorithm
https://en.wikipedia.org/wiki/Tree_rearrangement
    - Picture of nearest neighbor interchange made the concept click


