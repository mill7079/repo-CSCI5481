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


3.4

The runtime for 3.1 and 3.2 should be the same, though the space usage of 3.2 is lower than that of 3.1 as
    the program is only storing a subset of the suffix array values. In an ideal implementation, then, the
    time usage of 3.3 is less than that of both 3.1 and 3.2, as the checkpoints make the LF lookup constant
    time. The space complexity of 3.3 is the same as that of 3.2 (assuming the algorithm in 3.3 is based off
    3.2) as it's still using the subset of the suffix array; adding the checkpoints presumably consumes a
    constant amount of space, so the space usage shouldn't increase due to that.


Issues:

I was unsure how to implement the checkpoint version of FM indexing with the way in which I'm
    performing the LF lookup - I precalculate all of the lookup values before the algorithm begins (because
    the first two indexing algorithms and the reverse bwt are slow to an ungodly degree otherwise), so each
    LF call is already constant time (there's no counting going on in the call so including checkpoints doesn't
    do anything). I've included the algorithm with the checkpoints calculated but it functions the same as the
    original version.


Resources:

https://www.cs.jhu.edu/~langmea/resources/lecture_notes/bwt_and_fm_index.pdf
- Super helpful in helping me understand what I wasn't understanding about FM index from the class slides (i.e. the row lookups)