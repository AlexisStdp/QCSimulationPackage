### First draft of a Quantum Circuit object

import numpy as np


def kron(*args):
    ## multiple kronecker product
    qb = np.array([[1.0]])
    for q in args:
        qb = np.kron(qb, q)
    return qb

def n_kron(n, vector):
    ## n kronecker product with itself
    ret = np.array([[1.0]])
    for _ in range(n):
        ret = np.kron(ret, vector)
    return ret

def dot(*args):
    # multiple dot products
    qb = 1
    for q in args:
        qb = np.dot(qb, q)
    return qb

def gate_operator(O, i, n):
    # O is the matrix representation of the gate
    # It should be at the i-th place
    # There are n qubits in total
    I = np.eye(2)
    return kron(n_kron(i, I), O, n_kron(n-i-1, I))

def gate_multiple_operator(O, args, n):
    # O is the matrix representation of the gate
    # It should be at the places indicated by the elements of args
    # There are n qubits in total
    I = np.eye(2)
    ret = np.array([[1.0]])
    for i in range(n):
        if i+1 in args: # because of qubit notation
            ret = np.kron(ret, O)
        else:
            ret = np.kron(ret, I)
    return ret

def prepare_projector(P, i, n):
    # P is the projector (generally of the form [[1,0],[0,0]]
    # Measurement on the i-th qubit
    # n qubits in the system
    I = np.eye(2)
    return kron(n_kron(n-i-1, I), P, n_kron(i, I))


class QCircuit:
    ### states are updated after every gate addition
    def __init__(self, number_of_qubits = 1):
        self.number_of_qubits = number_of_qubits
        self.state_zero = np.array([[1.0],[0.0]])
        # self.state_one = np.array([[0.0],[1.0]])
        self.initial_state = n_kron(self.number_of_qubits, self.state_zero)
        self.state = self.initial_state
        self.track_gates1 = [[] for _ in range(self.number_of_qubits)]
        self.track_gates2 = [[] for _ in range(self.number_of_qubits)]
        self.basis = [('{:0'+str(self.number_of_qubits)+'b}').format(i) for i in range(2**self.number_of_qubits)]
        # self.I = np.eye(2)
        
    def one_X(self, i = 1): # first qubit per default
        # add X gate to i-th qubit
        i -= 1
        self.X = np.array([[0.0, 1.0], [1.0, 0.0]])
        self.state = np.dot(gate_operator(self.X, i, self.number_of_qubits), self.state)
        for j in range(self.number_of_qubits):
            if j == i:
                self.track_gates1[j].append("X")
                self.track_gates2[j].append("X")
            else:
                self.track_gates2[j].append("-")

    def X(self, *args):
        # add X gate to multiple qubits
        if len(args) == 0:
            args = [1]
        self.X = np.array([[0.0, 1.0], [1.0, 0.0]])
        self.state = np.dot(gate_multiple_operator(self.X, args, self.number_of_qubits), self.state)
        for j in range(self.number_of_qubits):
            if j+1 in args:
                self.track_gates1[j].append("X")
                self.track_gates2[j].append("X")
            else:
                self.track_gates2[j].append("-")

    def H(self, *args):
        # add H gate to multiple qubits
        if len(args) == 0:
            args = [1]
        self.H = 1.0 / 2**.5 *  np.array([[1, 1], [1, -1]])
        self.state = np.dot(gate_multiple_operator(self.H, args, self.number_of_qubits), self.state)
        for j in range(self.number_of_qubits):
            if j+1 in args:
                self.track_gates1[j].append("H")
                self.track_gates2[j].append("H")
            else:
                self.track_gates2[j].append("-")

    def CNOT(self, control=1, target=2):
        # add CNOT gate w.r.t. control and target (both should be valid qubits)
        # for now, the control and the target have to be next to each other
        if abs(control - target) > 1:
            print("Warning, the control and target should be next to eachother, nothing added")
        elif control == target:
            print("Warning, the control and target should be different, nothing added")
        elif control < target:
            self.CNOT = np.array([[1.0, 0.0, 0.0, 0.0],
                                  [0.0, 1.0, 0.0, 0.0],
                                  [0.0, 0.0, 0.0, 1.0],
                                  [0.0, 0.0, 1.0, 0.0],])
            self.state = np.dot(gate_operator(self.CNOT, control-1, self.number_of_qubits-1), self.state)
        else:
            self.CNOT = np.array([[0.0, 1.0, 0.0, 0.0],
                                  [1.0, 0.0, 0.0, 0.0],
                                  [0.0, 0.0, 1.0, 0.0],
                                  [0.0, 0.0, 0.0, 1.0],])
            self.state = np.dot(gate_operator(self.CNOT, target-1, self.number_of_qubits-1), self.state)
        if abs(control - target) == 1:
            for j in range(self.number_of_qubits):
                if j+1 == control:
                    self.track_gates1[j].append("ctrl")
                    self.track_gates2[j].append("ctrl")
                elif j+1 == target:
                    self.track_gates1[j].append("CNOT")
                    self.track_gates2[j].append("CNOT")
                else:
                    self.track_gates2[j].append("----")

    def measure(self, i = 1):
        # add measurement gate at i-th qubit
        i -= 1
        self.P = np.dot(self.state_zero, self.state_zero.T)
        prob = dot(np.conjugate(self.state).T,gate_operator(self.P,i,self.number_of_qubits),self.state)
        if np.random.rand() < prob:
            self.state = np.dot(gate_operator(self.P,i,self.number_of_qubits),self.state) / np.sqrt(prob)
        elif prob > 0.0:
            self.state_one = np.array([[0.0],[1.0]])
            self.P1 = np.dot(self.state_one, self.state_one.T)
            self.state = np.dot(gate_operator(self.P1,i,self.number_of_qubits),self.state) / np.sqrt(prob)
        for j in range(self.number_of_qubits):
            if j == i:
                self.track_gates1[j].append("M")
                self.track_gates2[j].append("M")
            else:
                self.track_gates2[j].append("-")

    def draw(self):
        print("First Drawing")
        for _,q in enumerate(self.track_gates1):
            ret = "|0> --- " + " --- ".join(q)
            print(ret)
        print("Second Drawing")
        for _,q in enumerate(self.track_gates2):
            ret = "|0> --- " + " --- ".join(q)
            print(ret)

    def reinitialize(self):
        # carefull for the previous states will be lost
        self.state = self.initial_state
        self.track_gates1 = [[] for _ in range(self.number_of_qubits)]
        self.track_gates2 = [[] for _ in range(self.number_of_qubits)]

    def dirac(self):
        # returns a nice description of the state of the system
        equation = "|Psi> = "
        for i,s in enumerate(self.state):
            equation += str(s[0]) + '|' + self.basis[i] + '> + '
        print(equation[:-2])
        equation = "|Psi> = "
        for i,s in enumerate(self.state):
            if s > 0.0:
                equation += str(s[0]) + '|' + self.basis[i] + '> + '
        print(equation[:-2])

## Introducing QCChain

class QCChain:
    ### states are updated after every simulation call, allows for more flexibility, and introduces a "running time"
    def __init__(self, number_of_qubits = 1):
        self.number_of_qubits = number_of_qubits
        self.state_zero = np.array([[1.0],[0.0]])
        # self.state_one = np.array([[0.0],[1.0]])
        self.initial_state = n_kron(self.number_of_qubits, self.state_zero)
        self.state = self.initial_state
        self.track_gates1 = [[] for _ in range(self.number_of_qubits)]
        self.track_gates2 = [[] for _ in range(self.number_of_qubits)]
        self.basis = [('{:0'+str(self.number_of_qubits)+'b}').format(i) for i in range(2**self.number_of_qubits)]
        # self.I = np.eye(2)
        self.count = {basis:0 for basis in self.basis}

    def X(self, *args):
        # add X gate to multiple qubits, default to first qubit
        if len(args) == 0:
            args = [1]
        for j in range(self.number_of_qubits):
            if j+1 in args:
                self.track_gates1[j].append("X")
                self.track_gates2[j].append("X")
            else:
                self.track_gates2[j].append("-")

    def H(self, *args):
        # add H gate to multiple qubits
        if len(args) == 0:
            args = [1]
        for j in range(self.number_of_qubits):
            if j+1 in args:
                self.track_gates1[j].append("H")
                self.track_gates2[j].append("H")
            else:
                self.track_gates2[j].append("-")

    def CNOT(self, control=1, target=2):
        # add CNOT gate w.r.t. control and target (both should be valid qubits)
        # for now, the control and the target have to be next to each other
        if abs(control - target) > 1:
            print("Warning, the control and target should be next to eachother, nothing added")
        elif control == target:
            print("Warning, the control and target should be different, nothing added")
        else:
            for j in range(self.number_of_qubits):
                if j+1 == control:
                    self.track_gates1[j].append("ctrl")
                    self.track_gates2[j].append("ctrl")
                elif j+1 == target:
                    self.track_gates1[j].append("CNOT")
                    self.track_gates2[j].append("CNOT")
                else:
                    self.track_gates2[j].append("----")

    def measure(self, *args):
        # add measurement gate at i-th qubit
        if len(args) == 0:
            args = [1]
        for j in range(self.number_of_qubits):
            if j+1 in args:
                self.track_gates1[j].append("M")
                self.track_gates2[j].append("M")
            else:
                self.track_gates2[j].append("-")

    def draw(self):
        print("First Drawing")
        for _,q in enumerate(self.track_gates1):
            ret = "|0> --- " + " --- ".join(q)
            print(ret)
        print("Second Drawing")
        for _,q in enumerate(self.track_gates2):
            ret = "|0> --- " + " --- ".join(q)
            print(ret)

    def reinitialize(self):
        # carefull for the previous states will be lost
        self.state = self.initial_state
        self.track_gates1 = [[] for _ in range(self.number_of_qubits)]
        self.track_gates2 = [[] for _ in range(self.number_of_qubits)]

    def add(self, gates=[['X']], qubit=[1], place=[0]):
        # special method that adds a gate or several gate to specified place.
        # Example: q.add([['X','H'],['X'],['H']],[1,5,6],[-1,0,1])
        #       this will add two gates X and H to the first qubit before the last gate,
        #       X to the fifth qubit after all the other added gates,
        #       H to the sixth qubit at first place (so before all the other).
        for j in range(self.number_of_qubits):
            if j+1 in qubit:
                i = qubit.index(j+1)
                if place[i] == 0:
                    for gate in gates[i]:
                        self.track_gates1[j].append(gate)
                        self.track_gates2[j].append(gate)
                if place[i] > 0:
                    for gate in gates[i]:
                        self.track_gates1[j].insert(place[i]-1,gate)
                        self.track_gates2[j].insert(place[i]-1,gate)
                if place[i] < 0:
                    for gate in gates[i]:
                        self.track_gates1[j].insert(place[i],gate)
                        self.track_gates2[j].insert(place[i],gate)
            else:
                self.track_gates2[j].append("-")

    def delq(self,qubit=[1],place=[0]):
        # allows to delete gates at specified places
        for j in range(self.number_of_qubits):
            if j+1 in qubit:
                i = qubit.index(j+1)
                if place[i] == 0:
                    del self.track_gates1[j][-1]
                    del self.track_gates2[j][-1]
                else:
                    del self.track_gates1[j][place[i]-1]
                    del self.track_gates2[j][place[i]-1]

    def simulate(self):
        # simulate the circuit! uses the second tracker
        self.state = self.initial_state
        for j,_ in enumerate(self.track_gates2[0]):
            queue = [self.track_gates2[e][j] for e in range(self.number_of_qubits)]
            app = []
            for i,g in enumerate(queue):
                if g not in ['-','ctrl','CNOT','----']:
                    app.append(i+1)
                    c = g
                elif g == 'ctrl':
                    control = i
                elif g == 'CNOT':
                    target = i
                    c = g
            if c == 'X':
                self.X = np.array([[0.0, 1.0], [1.0, 0.0]])
                self.state = np.dot(gate_multiple_operator(self.X, app, self.number_of_qubits), self.state)
            elif c == 'H':
                self.H = 1.0 / 2**.5 *  np.array([[1, 1], [1, -1]])
                self.state = np.dot(gate_multiple_operator(self.H, app, self.number_of_qubits), self.state)
            elif c == 'M':
                for i in app:
                    self.P = np.dot(self.state_zero, self.state_zero.T)
                    prob = dot(np.conjugate(self.state).T,gate_operator(self.P,i-1,self.number_of_qubits),self.state)
                    if np.random.rand() < prob:
                        self.state = np.dot(gate_operator(self.P,i-1,self.number_of_qubits),self.state) / np.sqrt(prob)
                    elif prob > 0.0:
                        self.state_one = np.array([[0.0],[1.0]])
                        self.P1 = np.dot(self.state_one, self.state_one.T)
                        self.state = np.dot(gate_operator(self.P1,i-1,self.number_of_qubits),self.state) / np.sqrt(prob)
            elif c == 'CNOT':
                if control < target:
                    self.CNOT = np.array([[1.0, 0.0, 0.0, 0.0],
                                          [0.0, 1.0, 0.0, 0.0],
                                          [0.0, 0.0, 0.0, 1.0],
                                          [0.0, 0.0, 1.0, 0.0],])
                    self.state = np.dot(gate_operator(self.CNOT, control, self.number_of_qubits-1), self.state)
                else:
                    self.CNOT = np.array([[0.0, 1.0, 0.0, 0.0],
                                          [1.0, 0.0, 0.0, 0.0],
                                          [0.0, 0.0, 1.0, 0.0],
                                          [0.0, 0.0, 0.0, 1.0],])
                    self.state = np.dot(gate_operator(self.CNOT, target, self.number_of_qubits-1), self.state)

    def run(self, shots=1):
        # allows to simulate multiple times
        self.results = []
        self.count = {basis:0 for basis in self.basis}
        for i in range(shots):
            self.simulate()
            self.results.append([i,self.state])
            for i,s in enumerate(self.state):
                if s == [1]:
                    self.count[self.basis[i]] += 1
    
    def dirac(self):
        # returns a nice description of the state of the system
        equation = "|Psi> = "
        for i,s in enumerate(self.state):
            equation += str(s[0]) + '|' + self.basis[i] + '> + '
        print(equation[:-2])
        equation = "|Psi> = "
        for i,s in enumerate(self.state):
            if s > 0.0:
                equation += str(s[0]) + '|' + self.basis[i] + '> + '
        print(equation[:-2])

    def plot(self):
        # uses matplotlib
        import matplotlib.pyplot as plt
        plt.bar(list(self.count.keys()),self.count.values())
        plt.title('Results of the quantum experiment')
        plt.xlabel('Basis states')
        plt.ylabel('Count')
        plt.show()

""" !EPR pair!
qc = QCircuit(2)
qc.H(1)
print(qc.state)
qc.CNOT(1,2)
print(qc.state)
qc.measure(1)
qc.measure(2)
print(qc.state)
qc.draw()
print(qc.number_of_qubits)
print(qc.track_gates1)
print(qc.track_gates2)
"""
