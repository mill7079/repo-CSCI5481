import sys
sys.path.append(".");
from NWA import Alignment

if __name__ == '__main__':
    if (len(sys.argv) != 4):
        print("Too few arguments.")
        exit(1);
        
    _, file1, file2, matchfile = sys.argv;
    
    sys.path.append(".");
    align = Alignment(file1, file2, -5, -1, matchfile);
    s1, s2, score = align.anchored_nwa("blosum");
##    result = align.print_alignment();

    print(s1)
    print(s2)
    print("Score:", score);
##    print(result[0])
##    print(result[1])
