from math import inf

# all possibilities for sankoff scoring
alph = ['AA', 'TT', 'GG', 'CC', 'NN']


class Node:
    def __init__(self, sequence=""):
        self.seq = sequence.upper()  # store sequence in uppercase to avoid conflicts
        self.visited = False  # used for printing the tree
        self.connections = []  # all nodes that connect to this node in tree, set during join

        # useful things for sankoff
        self.scored = False  # used for scoring during Sankoff
        self.scores = [inf, inf, inf, inf, inf]

        self.parent = None  # explicit refs to parent and child nodes
        self.children = []

        if self.seq != "":  # leaf node, so set score accordingly
            self.scores[alph.index(self.seq)] = 0
            self.scored = True

    # connect two nodes - double linking
    def connect(self, other):
        if other not in self.connections:
            self.connections.append(other)

        if self not in other.connections:
            other.connections.append(self)

    # disconnect two nodes
    def disconnect(self, other):
        if other in self.connections:
            self.connections.remove(other)

        if self in other.connections:
            other.connections.remove(self)

    def __str__(self):
        return "<Node: " + self.seq + ", Connections: " + str(len(self.connections)) #+ \
               #", num kids " + str(len(self.children)) + ">"

    def __repr__(self):
        return self.__str__()
