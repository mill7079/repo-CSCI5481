import sys
sys.path.append(".");
from NWA import Alignment

if __name__ == '__main__':
    _, file1, file2 = sys.argv;
    
    sys.path.append(".");
    align = Alignment(file1, file2, -5, -1, None);
    align.needleman_wunsch("blosum");
    align.print_alignment();


