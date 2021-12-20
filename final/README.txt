Allison Miller
Final Project - Option 3 (Parsimony Analysis)
5278509

To run:
The program includes two command-line flags, one for the number of individuals to read from each data file (-i or --num_inds), and
    one to indicate the start of the list of files (-f or --files). Examples:

To run an analysis on two files with 20 individuals from each file:
    python3 main.py -i 20 -f data/genotypes_chrM_ASW_phase3.2_consensus.b36_fwd.txt data/genotypes_chrM_JPT_phase3.2_consensus.b36_fwd.txt

To run an analysis on three files with 50 individuals from each file:
    python3 main.py --num_inds 50 --files data/genotypes_chrM_ASW_phase3.2_consensus.b36_fwd.txt data/genotypes_chrM_JPT_phase3.2_consensus.b36_fwd.txt data/genotypes_chrM_CEU_phase3.2_consensus.b36_fwd.txt

Each file should be a path to the file (i.e. I had the files all in a directory named 'data' while testing, hence the paths). Number
    of individuals should be an integer.

Issues:

Took longer than expected, as I had to completely re-implement the Sankoff algorithm. My homework 3 implementation used a different
    tree structure, since I used an online tool for neighbor joining, and I couldn't figure out how I would go about creating that
    particular tree structure from scratch. The neighbor-joining algorithm seems a little off at times, but I can't figure out why
    that might be.
