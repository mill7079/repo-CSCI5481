Allison Miller
mill7079
Homework 4

All functions listed in the homework writeup are contained in main.py.


To run:

python3 main.py


The main method automatically runs the functions on the file names suggested in the homework writeup;
    however, these files aren't included in the submission as per the instructions, so if you want to
    run the code you'll need to place those files in the files directory (as all the paths look for files
    in that directory), and if you want to test new files, you'll need to change the filenames in the code.

Additionally, for ease of grading, I've left the source file out in the main project folder instead of storing
    it in its own folder; the IDE I use starts looking for files in the main project directory anyway, so it
    saves some hassle of fixing relative paths before submission.

All files mentioned in the rubric are included and named appropriately; there's an additional file called
    reversedSample1.fa that's the reverse_bwt of sample1.fa for testing purposes.


Issues:

The indices in mapping_FMhybrid.txt are incorrect, but I couldn't resolve the issue (something with how the
    indices relate between the suffix array and rows).

Additionally, I was unsure how to implement the checkpoint version of FM indexing with the way in which I'm
    performing the LF lookup - I precalculate all of the lookup values before the algorithm begins (because
    it's slow to an ungodly degree otherwise), so each LF is already constant time (there's no counting going
    on in the call so including checkpoints doesn't do anything). I've included the algorithm with the
    checkpoints calculated but it functions the same as the original version.