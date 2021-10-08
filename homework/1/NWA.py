#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3

import sys

class Node:  # used for backtracing alignment
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
        file = open("blosum.txt");
        letters = file.readline();
        letters = letters.split();
        numbers = [];
        line= file.readline();
        while line != "":
            numbers.append(line.split()[1:]);
            line = file.readline();

        self.letters = letters;
        self.numbers = numbers;

    def score(self, xi, yj):
        i = self.letters.index(xi.upper());
        j = self.letters.index(yj.upper());
        return int(self.numbers[i][j]);
    

# class for actual sequence alignment
class Alignment:
    def __init__(self, file1, file2, matches=None):
        self.seq1 = open(file1);
        self.seq1.readline();
        self.seq1 = self.seq1.readline().strip();

        self.seq2 = open(file2);
        self.seq2.readline();
        self.seq2 = self.seq2.readline().strip();

        print(self.seq1);
        print(self.seq2);

        self.blosum = Blosum();
        self.match = 1;
##        self.mismatch = -3;
##        self.gap = -2;
        self.mismatch = -1;
        self.gap = -1;
        
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

    def similar(self, xi, yj):
        if xi == yj:
            return self.match;  # match
        return self.mismatch;  # mismatch

    # returns list of max values given three
##    def max_values(self, a, b, c):
##        ret = [];
##        m = max(a, b, c);
##        if (m == a):
##            ret.append(a);
##        if (m == b):
##            ret.append(b);
##        if (m == c):
##            ret.append(c);
##
##        return ret;

    # nodelist is list of max calc score prevnodes with directions
    def max_node(self, nodelist):
        num = nodelist[0][0].score
        node = nodelist[0]
        for pair in nodelist:
            if (pair[0].score > num):
                num = pair[0].score
                node = pair;

        return pair;
        

    def find_max(self, r, c, score_func):
        case1 = self.scores[r-1][c-1].score + score_func(self.seq1[c-1], self.seq2[r-1]);
        case2 = self.scores[r-1][c].score + self.gap;
        case3 = self.scores[r][c-1].score + self.gap;
        sc = max(case1, case2, case3);
##        numvals = len(self.max_values(case1, case2, case3))

        # god this is gross my apologies
##        if numvals == 1:
##            if sc == case1:
##                self.scores[r][c].link(self.scores[r-1][c-1], "D");
##            elif sc == case2:
##                self.scores[r][c].link(self.scores[r-1][c], "U");
##            else:
##                self.scores[r][c].link(self.scores[r][c-1], "L");
##        else:
##            nodelist = [];
##            if sc == case1:
##                nodelist.append([self.scores[r-1][c-1],"D"]);
##            if sc == case2:
##                nodelist.append([self.scores[r-1][c], "U"]);
##            if sc == case3:
##                nodelist.append([self.scores[r][c-1], "L"]);
##
##            mnode = self.max_node(nodelist)
##            self.scores[r][c].link(mnode[0], mnode[1]);
        nodelist = [];
        if sc == case1:
            nodelist.append([self.scores[r-1][c-1],"D"]);
        if sc == case2:
            nodelist.append([self.scores[r-1][c], "U"]);
        if sc == case3:
            nodelist.append([self.scores[r][c-1], "L"]);

        if (len(nodelist) == 1):
            self.scores[r][c].link(nodelist[0][0], nodelist[0][1]);
        else:
            mnode = self.max_node(nodelist)
            self.scores[r][c].link(mnode[0], mnode[1]);

        return sc;

    def print_scores(self):
        for row in self.scores:
            for node in row:
                print(node, end="");
            print("\n");

    def needleman_wunsch(self, score_func):
        if (score_func == "basic"):
            score_func = self.similar;
        else:
            score_func = self.blosum.score;
            
        for r in range (1, len(self.scores)):  # r -> s2
            for c in range (1, len(self.scores[r])):  # c -> s1
                self.scores[r][c].set_score(self.find_max(r, c, score_func));

        self.print_scores();

    def print_alignment(self):
        s1 = len(self.seq1) - 1;
        s2 = len(self.seq2) - 1;
        nseq1 = "";
        nseq2 = "";
        temp = self.scores[-1][-1];

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
        
# print scores
##def print_scores(matrix):
##    for row in matrix:
##        for node in row:
##            print(node, end="");
##        print("\n",);
    

# returns score matrix
##def needleman_wunsch(seq1, seq2, match, gap, mismatch):
##    def s(xi, yj):
##        if xi == yj:
##            return match;
##        return mismatch;
##
##    def find_max(matrix, r, c):
##        # find score
##        case1 = matrix[r-1][c-1].score + s(seq1[c-1], seq2[r-1]);
##        case2 = matrix[r-1][c].score + gap;
##        case3 = matrix[r][c-1].score+ gap;
##        sc = max(case1, case2, case3);
##
##        # set link
##        if sc == case1:
##            matrix[r][c].link(matrix[r-1][c-1], "D");
##        elif sc == case2:
##            matrix[r][c].link(matrix[r-1][c], "U");
##        else:
##            matrix[r][c].link(matrix[r][c-1], "L");
##
##        return sc;
##        
##    # first is cols, second is rows - seq1 across top (cols), seq2 down side (rows)
##    score = [[0] * (len(seq1)+1) for i in range(len(seq2)+1)]
##    for i in range (0, len(score)):  # not sure how to do this more elegantly
##        for j in range (0, len(score[0])):
##            score[i][j] = Node();
##    
##    # initialize
##    for i in range (0, len(score[0])):
##        score[0][i].set_score(i * gap);
##
##    for j in range (0, len(score)):
##        score[j][0].set_score(j * gap);
##
##    # main iteration
##    for r in range (1, len(score)):  # r -> s2
##        for c in range (1, len(score[r])):  # c -> s1
##            score[r][c].set_score(find_max(score, r, c));
##
##    print_scores(score);
##
##    # find alignment
##    s1 = len(seq1) - 1;
##    s2 = len(seq2) - 1;
##    nseq1 = "";
##    nseq2 = "";
##    temp = score[-1][-1];
##
##    while (temp != None):
##        print(temp);
##        if (temp.dir == "D"):
##            nseq1 = seq1[s1] + nseq1;
##            nseq2 = seq2[s2] + nseq2;
##            s1-=1; 
##            s2-=1;
##        elif temp.dir == "U":
##            nseq1 = "-" + nseq1;
##            nseq2 = seq2[s2] + nseq2;
##            s2-=1;
##        elif temp.dir == "L":
##            nseq1 = seq1[s1] + nseq1;
##            nseq2 = "-" + nseq2; 
##            s1-=1;
##            
##        temp = temp.trace();
##
##    print(nseq1);
##    print(nseq2);
##    print(len(nseq1));
##    print(len(nseq2));
##    print(len(seq1));
##    print(len(seq2));
    
##x = Node();
##y = Node();
##z = Node();
##
##x.set_score(1);
##y.set_score(2);
##z.set_score(3);
##
##y.link(x);
##z.link(y);
##
##temp = z;
##while (temp != None):
##    print(temp.score);
##    temp = temp.trace();

        
##d2 = [[1,2],[2,3],[3,4]]
##print_scores(d2);

##match = 1;
##gap = -2;
##mismatch = -3;

if __name__ == '__main__':
    _, seq1, seq2 = sys.argv
    print(seq1)
    print(seq2)

    a = Alignment(seq1, seq2)
    a.needleman_wunsch("basic");
    a.print_alignment();

##    print(a.max_values(1, 2, 3));
##    print(a.max_values(1, 3, 3));
##    print(a.max_values(3,3,3));
#0 11 3
##    b = Blosum();
##    print(b.score("A","V"));
##    print(b.score("W","W"));
##    print(b.score("B","N"));

##    seq1 = open(seq1);
##    seq2 = open(seq2);
##    seq1.readline();
##    seq2.readline();
##
##    needleman_wunsch(seq1.read(), seq2.read(), match, gap, mismatch);
##    needleman_wunsch("aaagaattca", "aaatca", match, -1, -1);
    
