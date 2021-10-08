import sys
import urllib.request

# used for backtracing alignment; holds backtraced node and direction of node
class Node:
    def __init__(self):
        self.score = 0;
        self.back_ptr = None;
        self.dir = "*";

    def trace(self):
        return self.back_ptr;

    def link(self, other, direction):
        self.back_ptr = other;
        self.dir = direction;

    def set_score(self, score):
        self.score = score;

    def __str__(self):
        return str(self.score) + str(self.dir) + " ";
    

# holds BLOSUM62 matrix
class Blosum:
    def __init__(self):
        # this is kinda gross but I kept getting certificate errors if I tried reading
            # directly from the website. my apologies
        matrix = "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  * \n\
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 \n\
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 \n\
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 \n\
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 \n\
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 \n\
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 \n\
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 \n\
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 \n\
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 \n\
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 \n\
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 \n\
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 \n\
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 \n\
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 \n\
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 \n\
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 \n\
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 \n\
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 \n\
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 \n\
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 \n\
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 \n\
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 \n\
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 \n\
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1" 

        lines = matrix.split("\n");
        letters = lines[0].strip().split();
        numbers = [];
        for i in range (1, len(lines)):
             numbers.append(lines[i].strip().split()[1:]);

        self.letters = letters;
        self.numbers = numbers;
        
    # look up score from matrix
    def score(self, xi, yj):
        i = self.letters.index(xi.upper());
        j = self.letters.index(yj.upper());
        return int(self.numbers[i][j]);
    

# class for actual sequence alignment
class Alignment:
    def __init__(self, file1, file2, gap, mismatch, matches):
        self.seq1 = open(file1);
        self.seq1.readline();
        self.seq1 = self.seq1.readline().strip();

        self.seq2 = open(file2);
        self.seq2.readline();
        self.seq2 = self.seq2.readline().strip();

##        print(self.seq1);
##        print(self.seq2);

        self.blosum = Blosum();
        self.match = 1;
        self.mismatch = mismatch;
        self.gap = gap;
        
        self.scores = [[0] * (len(self.seq1)+1) for i in range(len(self.seq2)+1)];
        for i in range (0, len(self.scores)):  # not sure how to do this more elegantly
            for j in range (0, len(self.scores[i])):
                self.scores[i][j] = Node();
                
        for i in range (0, len(self.scores)):  # initialize score matrix
            self.scores[i][0].set_score(i * self.gap);

        for j in range (0, len(self.scores[0])):
            self.scores[0][j].set_score(j * self.gap);

        self.matches = matches;
        self.align_seq1 = "";
        self.align_seq2 = "";

    # basic scoring function
    def similar(self, xi, yj):
        if xi == yj:
            return self.match;  # match
        return self.mismatch;  # mismatch

    # nodelist is list of max calc score prevnodes with directions
    def max_node(self, nodelist):
        num = nodelist[0][0].score
        node = nodelist[0]
        for pair in nodelist:
            if (pair[0].score > num):
                num = pair[0].score
                node = pair;

        return pair;
        
    # calculate score for a node and link the previous node for alignment
    def find_max(self, r, c, score_func):
        case1 = self.scores[r-1][c-1].score + score_func(self.seq1[c-1], self.seq2[r-1]);
        case2 = self.scores[r-1][c].score + self.gap;
        case3 = self.scores[r][c-1].score + self.gap;
        sc = max(case1, case2, case3);

        # add all to list with their directions
        nodelist = [];
        if sc == case1:
            nodelist.append([self.scores[r-1][c-1],"D"]);
        if sc == case2:
            nodelist.append([self.scores[r-1][c], "U"]);
        if sc == case3:
            nodelist.append([self.scores[r][c-1], "L"]);

        # if only one max calculated score, use that for backtrace
        # otherwise find max of original values of max calc score nodes
        if (len(nodelist) == 1):
            self.scores[r][c].link(nodelist[0][0], nodelist[0][1]);
        else:
            mnode = self.max_node(nodelist)
            self.scores[r][c].link(mnode[0], mnode[1]);

        return sc;

    # print score table in a reasonable format
    def print_scores(self):
        for row in self.scores:
            for node in row:
                print(node, end="");
            print("\n");

    # implement algorithm and print the score table
    def needleman_wunsch(self, score_func):
        if (score_func == "basic"):
            score_func = self.similar;
        else:
            score_func = self.blosum.score;
            
        for r in range (1, len(self.scores)):  # r -> s2
            for c in range (1, len(self.scores[r])):  # c -> s1
                self.scores[r][c].set_score(self.find_max(r, c, score_func));

        self.print_scores();

    # print the final alignment
    def print_alignment(self):
        s1 = len(self.seq1) - 1;
        s2 = len(self.seq2) - 1;
        nseq1 = "";
        nseq2 = "";
        temp = self.scores[-1][-1];

        # backtrace
        while (temp != None):
            if (temp.dir == "D"):
                nseq1 = self.seq1[s1] + nseq1;
                nseq2 = self.seq2[s2] + nseq2;
                s1-=1; 
                s2-=1;
            elif temp.dir == "U":
                nseq1 = "-" + nseq1;
                nseq2 = self.seq2[s2] + nseq2;
                s2-=1;
            elif temp.dir == "L":
                nseq1 = self.seq1[s1] + nseq1;
                nseq2 = "-" + nseq2; 
                s1-=1;
                
            temp = temp.trace();

        print(nseq1);
        print(nseq2);
        return nseq1, nseq2;

##    def anchored_nwa(self, score_func):
##        file = open(self.matches);
####        line = file.readline();
##        segments = [];
####        while (line != "" and line != "\n"):
####            split = line.split();
####            print(split);
####            line = file.readline();
##        for line in file:
##            split = line.split();
##            print(int(split[1]));
##            nseq1 = self.seq1[int(split[0]) - 1 : int(split[1])];
##            print(nseq1)
##
##        print(segments)


# main
if __name__ == '__main__':
    if (len(sys.argv) == 3):
        _, seq1, seq2 = sys.argv
        matchfile = None;
    elif (len(sys.argv) == 4):
        _, seq1, seq2, matchfile = sys.argv
    else:
        print("Wrong number of arguments.")
        exit(-1);
        
    print(seq1)
    print(seq2)

    a = Alignment(seq1, seq2, -5, -3, matchfile)
##    if (matchfile == None):
##        a.needleman_wunsch("basic");
##    else:
##        a.anchored_nwa("basic");
    a.needleman_wunsch("blosum");
        
    a.print_alignment();
    
