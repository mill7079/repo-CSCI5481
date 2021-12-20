class Node:
    def __init__(self, sequence=""):
        self.seq = sequence
        self.visited = False
        self.connections = []

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
        return "Node: " + self.seq + ", Connections: " + str(len(self.connections))

    def __repr__(self):
        return self.__str__()
