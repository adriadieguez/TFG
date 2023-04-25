# Required packages
from Bio import Phylo
from io import StringIO
from Bio import SeqIO
from sympy import symbols, Eq, solve
from scipy.linalg import expm
import numpy as np
import sympy as sp
import sys

tree = sys.argv[1]
case = 1

num_lengths = len(sys.argv)
lengths = []
if (sys.argv[2].startswith("L")):
    for i in range(2, num_lengths):
        lengths.append(int(sys.argv[i][1:]))
        case = 2
else:
    t = int(sys.argv[2])
    L = int(sys.argv[3])


def generate_alignment(length, distribution):
    """
    Generates an alignment of length `length` using a given `distribution`.
    """
    seq = ''
    for i in range(length):
        # Generate a random sample from the multinomial distribution
        nucleotide = np.random.choice(['A', 'G', 'C', 'T'], p=distribution)
        seq += nucleotide
    return seq

class Edge:
    def __init__(self, edge, transition_matrix=None):
        self.edge = edge
        self.transition_matrix = transition_matrix

class MM:
    def __init__(self, source, target, matrix):
        self.source = source
        self.target = target
        self.matrix = matrix

def get_matrix_from_dict(d):
    """
    Return a matrix from a dictionary.
    """
    Q2 = np.zeros((4,4))
    coefficients = list(d.values())
    for i in range(4):
        for j in range(4):
            Q2[i,j] = coefficients[i*4+j]
    return Q2

def alpha(new_distribution, Q, i, k):
    return min(1, (new_distribution[k]*Q[k,i])/(new_distribution[i]*Q[i,k]))

def get_M2(new_distribution,d2, l):
    """
    Metropolis - Hastings Markov chain Monte Carlo implementation to get M2
    """
    P = np.zeros((4,4))
    iter = True
    while iter:
        Q = np.zeros((4,4))
        i=0
        while i<4:
            dir = np.ones(4)
            dir[i] = (50/np.sqrt(l))*np.exp(-l)
            R = np.random.dirichlet(dir)
            if R[i] > 0.4:
                Q[i,:] = R
                i = i + 1
        for i in range(4):
            for j in range(4):
                if i == j:
                    sum = 0
                    for k in range(4):
                        if k != i:
                            sum += (Q[i,k] * (1 - alpha(new_distribution,Q,i,k)))
                    P[i,j] = Q[i,i] + sum
                else:
                    P[i,j] = Q[i,j]*alpha(new_distribution,Q,i,j)
        assert (np.abs(np.sum(new_distribution - np.matmul(new_distribution,P)))) < 10**-6
        # Adjust the matrix diagonalising
        vaps, _ = np.linalg.eig(P)
        vaps = sorted(vaps, reverse=True)
        A = symbols('A')
        eq = Eq(-d2+(((1-A)*vaps[1]+A)*((1-A)*vaps[2]+A)*((1-A)*vaps[3]+A)),0)
        sol = solve(eq, A)
        # We only want the real solution between 0 and 1
        res = 0
        for s in sol:
            if s.is_real and s > 0 and s < 1:
                res = s
                res = np.float64(res)
                P = (1-res)*P + res*np.identity(4)
                iter = False
                break
            elif s.is_complex:
                b = np.imag(s)
                a = sp.re(s)
                if np.abs(b) < 10**-20 and a > 0 and a < 1:
                    res = sp.re(s)
                    res = np.float64(res)
                    P = (1-res)*P + res*np.identity(4)
                    iter = False
                    break
    return P

def generate_random_matrix(distribution, l):
    """
    Returns the markov chain matrix M = M1 · M2 given a branch length
    """
    res = 1
    # Compute M1
    while res >= 1:
        M1 = np.zeros((4,4))
        i=0
        while i<4:
            dir = np.ones(4)
            dir[i] = (50/np.sqrt(l))*np.exp(-l)
            R = np.random.dirichlet(dir)
            if R[i] > 0.4:
                M1[i,:] = R
                i = i + 1

        new_distribution = np.matmul(distribution,M1)
        D = np.diag(distribution)
        D_ = np.diag(new_distribution)
        res = np.exp(-l)*np.sqrt(np.linalg.det(D_))/np.sqrt(np.linalg.det(D))
        detM1 = np.linalg.det(M1)
        if detM1 > np.exp(-l)*np.sqrt(np.linalg.det(D_))/np.sqrt(np.linalg.det(D)):
            pass
        else:
            res = 1
    d2 = np.exp(-l)*np.sqrt(np.linalg.det(D_))/(detM1*np.sqrt(np.linalg.det(D)))
    # Obtain M2
    M2 = get_M2(new_distribution,d2,l)


    detM2 = np.linalg.det(M2)
    assert(np.abs(detM2 - d2) < 10**-6)
    M = np.matmul(M1,M2)
    return M

def generate_sequences(M, seq):
    """
    Given a sequence and the transition matrix it returns the next state sequence.
    """
    new_seq = ""
    for s in seq:
        if s == "A":
            new_seq += np.random.choice(['A', 'G', 'C', 'T'], p=M[0,:])
        elif s == "G":
            new_seq += np.random.choice(['A', 'G', 'C', 'T'], p=M[1,:])
        elif s == "C":
            new_seq += np.random.choice(['A', 'G', 'C', 'T'], p=M[2,:])
        else:
            new_seq += np.random.choice(['A', 'G', 'C', 'T'], p=M[3,:])
    return new_seq


def matrix_generation(tree, length, lengths):
    """
    Markov chain propagation. We call the prior methods for every tree edge.
    """
    node_distribution = dict()
    B = False
    while not B:
        R = np.random.dirichlet([1, 1, 1, 1])
        Res = all(ele > 0.2 and ele < 0.3 for ele in R)
        if Res == True:
            B = True

    node_distribution["Root"] = R

    path_t = tree
    tree_file = open(path_t, "r")
    tree = tree_file.read()
    tree = Phylo.read(StringIO(tree), "newick")
    # Change nodes names
    for idx, clade in enumerate(tree.get_nonterminals()):
        clade.name = "Node_" + str(idx) if idx > 0 else "Root"
    # Change leaves names
    for idx, clade in enumerate(tree.get_terminals()):
        clade.name = "Leaf_" + clade.name
    #Phylo.draw(tree)
    net = Phylo.to_networkx(tree) # to graph

    iter = 0
    edges = []
    for edge in net.edges():
        l = edge[1].branch_length
        new_edge = Edge(edge, generate_random_matrix(node_distribution[edge[0].name], l))
        edges.append(new_edge)
        node_distribution[edge[1].name] = np.matmul(node_distribution[edge[0].name],new_edge.transition_matrix)
        for i in range(4):
            assert(np.sum(new_edge.transition_matrix[i,:])<1.000000001 and np.sum(new_edge.transition_matrix[i,:])>0.999999999)
        iter += 1

    assert(iter == len(net.edges()))

    # We save all the transition matrices
    real_matrices = []
    for e in edges:
        real_matrices.append(MM(e.edge[0].name, e.edge[1].name, e.transition_matrix))
        # Check that matrices are DLC
        for i in range(4):
            assert(e.transition_matrix[i,i] == max(e.transition_matrix[i,:]))

    file_names = []
    if case == 1:


        # Once the transition matrices are computed, we generate t alignments of length L using the same matrices.
        for align in range(t):

            node_sequence = dict()
            node_sequence["Root"] = generate_alignment(L, node_distribution["Root"])
            i = 0
            for edge in net.edges():
                node_sequence[edge[1].name] = generate_sequences(edges[i].transition_matrix, node_sequence[edge[0].name])
                i += 1

            leaves_seq = {k: v for k, v in node_sequence.items() if k.startswith('L')}
            sequences_in_leaves = list(leaves_seq.values())
            keys_for_sequences = list(leaves_seq.keys())
            iter = 0
            file_name = str(len(sequences_in_leaves))+ "_leaves_" + str(len(sequences_in_leaves[0])) + "length_sequences_num" + str(align+1) +".fasta"
            file_names.append(file_name)
            file = open(file_name, "w")
            for seq in sequences_in_leaves:
                file.write(">Seq" + str(keys_for_sequences[iter]) + "\n" + seq + "\n")
                iter += 1
            file.close()

    else:
        for l in lengths:
            node_sequence = dict()
            node_sequence["Root"] = generate_alignment(l, node_distribution["Root"])
            i = 0
            for edge in net.edges():
                node_sequence[edge[1].name] = generate_sequences(edges[i].transition_matrix, node_sequence[edge[0].name])
                i += 1

            leaves_seq = {k: v for k, v in node_sequence.items() if k.startswith('L')}
            sequences_in_leaves = list(leaves_seq.values())
            keys_for_sequences = list(leaves_seq.keys())
            iter = 0
            file_name = str(len(sequences_in_leaves))+ "_leaves_" + str(len(sequences_in_leaves[0])) + "length_sequences.fasta"
            file_names.append(file_name)
            file = open(file_name, "w")
            for seq in sequences_in_leaves:
                file.write(">Seq" + str(keys_for_sequences[iter]) + "\n" + seq + "\n")
                iter += 1
            file.close()


    return real_matrices, file_names

if (case == 1):
    matrix_generation(tree, L, None)
else:
    matrix_generation(tree, None, lengths)
