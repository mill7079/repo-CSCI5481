from HMM import *

if __name__ == '__main__':
    hmm = HMM()
    print(hmm.emit('A+', 'A'))
    print(hmm.emit('T+', 'A'))
    hmm.viterbi('ATACGACA')
