Allison Miller
mill7079

How to Run:
python3 NWA.py seq1.fa seq2.fa
  -- runs the basic implementation, with -3 for mismatches and -2 for gaps
python3 NWA_B62.py seq1.fa seq2.fa
  -- runs the implementation using the BLOSUM62 matrix and -5 for gaps
python3 NWA_anchor.py seq1.fa seq2.fa matches.txt
  -- runs the anchored version with -5 for gaps

For the anchored algorithm, seq1 should always be the human version
 to line up with the matches.txt file.

I don't know how to run these without first including the python3 command, so
unfortunately it differs slightly from the submission instructions.

Issues/Notes:
The code is structured rather poorly, so my apologies if you read it. I had
been intending on having all the code in just one file as I misunderstood the
submission instructions; hence, all of the actual code is in NWA.py and the
other two files simply contain main methods for running each respective version
of the algorithm.

Additionally, I had originally been removing the newlines when reading the
sequences in from the .fa files; however, I did not notice until late that
whether they are included or not changes the result. I have commented out the
code to leave them in on lines 81 and 89; if the newlines are supposed to be
included, uncomment lines 81 and 89 and comment out lines 77-80 and 85-88. The
alignment scores in the alignment_results.txt file are the results from the
alignment without the newlines.
