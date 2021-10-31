from math import log, inf
from collections import deque


def print_m(matrix):
    for row in matrix:
        print(row)


class HMM:
    def __init__(self):
        p = 0.98
        q = 0.999

        self.states = ['A+', 'C+', 'G+', 'T+', 'A-', 'C-', 'G-', 'T-']
        self.trans_matrix = [[0.180*p, 0.274*p, 0.426*p, 0.120*p, (1-p)/4, (1-p)/4, (1-p)/4, (1-p)/4],
                             [0.171*p, 0.368*p, 0.274*p, 0.188*p, (1-p)/4, (1-p)/4, (1-p)/4, (1-p)/4],
                             [0.161*p, 0.339*p, 0.375*p, 0.125*p, (1-p)/4, (1-p)/4, (1-p)/4, (1-p)/4],
                             [0.079*p, 0.355*p, 0.384*p, 0.182*p, (1-p)/4, (1-p)/4, (1-p)/4, (1-p)/4],
                             [(1-q)/4, (1-q)/4, (1-q)/4, (1-q)/4, 0.300*q, 0.205*q, 0.285*q, 0.210*q],
                             [(1-q)/4, (1-q)/4, (1-q)/4, (1-q)/4, 0.322*q, 0.298*q, 0.078*q, 0.302*q],
                             [(1-q)/4, (1-q)/4, (1-q)/4, (1-q)/4, 0.248*q, 0.246*q, 0.298*q, 0.208*q],
                             [(1-q)/4, (1-q)/4, (1-q)/4, (1-q)/4, 0.177*q, 0.239*q, 0.292*q, 0.292*q]]
        # print(1/4)
        # csv = open('transition.csv')
        # csv.readline()
        # for line in csv:
        #     trans_matrix.append(line[3:].split(','))
        # print(self.trans_matrix)

    def transit(self, sk, sl):
        return self.trans_matrix[self.states.index(sk)][self.states.index(sl)]

    def emit(self, state, obs):
        # print(state, obs)
        return 1 if (state[0:1] == obs or obs == 'N') else 0

    def viterbi(self, seq):
        seq = seq.upper()
        v = [[0.0]*len(seq) for _ in range(len(self.states))]  # states = rows
        ptr = [['']*len(seq) for _ in range(len(self.states))]
        # initialization
        for i in range(0, len(self.states)):
            v[i][0] = log(1/len(self.states))

        # print(v)
        # iteration
        for i in range(1, len(seq)):
            # print(seq[i])
            for j in range(0, len(self.states)):
                argmax = ''
                max_prob = -inf
                emission = self.emit(self.states[j], seq[i])
                for k in range(0, len(self.states)):
                    # prob = v[k][i-1] * self.transit(self.states[k], self.states[j])
                    prob = v[k][i-1] + log(self.transit(self.states[k], self.states[j]))
                    # prob = log(v[k][i-1]) + log(self.transit(self.states[k], self.states[j]))
                    # print(k, i)
                    # if emission > 0:
                    #     print("prev:", v[k][i-1], "log transit:", log(self.transit(self.states[k], self.states[j])), "prob:", prob)
                    if prob > max_prob:
                        max_prob = prob
                        argmax = self.states[k]
                v[j][i] = max_prob + (log(emission) if emission > 0 else -inf)
                ptr[j][i] = argmax

        # print_m(v)
        # print_m(ptr)

        # termination - retrace
        max_end_prob = -inf
        idx = -1
        for i in range(0, len(v)):
            if v[i][-1] > max_end_prob:
                max_end_prob = v[i][-1]
                idx = i

        print(self.states[idx], max_end_prob)
        # path = [ptr[idx][-1]]
        # print(path)
        # print("starting path")
        # for i in range(len(seq) - 2, -1, -1):
        #     if i % 250000 == 0:
        #         print("step in path")
        #     idx = self.states.index(path[0])
        #     path.insert(0, ptr[idx][i])

        path = deque()
        path.append(ptr[idx][-1])
        print(path)
        print("starting path")
        for i in range(len(seq) - 2, -1, -1):
            if i % 250000 == 0:
                print("step in path")
            idx = self.states.index(path[0])
            path.insert(0, ptr[idx][i])

        # print(path)
        print("done with path")
        print(len(path))
        print(len(seq))
        return path

