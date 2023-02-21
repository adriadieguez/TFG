# Import modules
import sys
import pylab
import networkx, pylab
from Bio import Phylo
from io import StringIO

# Read tree
path = sys.argv[1]
tree_file = open(path, "r")
tree = tree_file.read()
tree = Phylo.read(StringIO(tree), "newick")

# Change nodes names
def tabulate_names(tree):
    for idx, clade in enumerate(tree.get_nonterminals()):
        if clade.name:
            clade.name = clade.name
        else:
            clade.name = "Internal_" + str(idx+1)
tabulate_names(tree)

# Tree to graph
net = Phylo.to_networkx(tree)

# Display graph properties
# Nodes
print(" ")
print("---------------------------------------------------------------------------------------------------")
print(" ")
print("Nodes:")
print(" ")
for node in net.nodes():
    print(node)

# Edges
print(" ")
print("---------------------------------------------------------------------------------------------------")
print(" ")
print("Edges:")
print(" ")
for edge in net.edges():
    print(edge)

# Degrees
print(" ")
print("---------------------------------------------------------------------------------------------------")
print(" ")
print("Degrees:")
print(" ")
degrees = [(node,val) for (node, val) in net.degree()]
for node_degree in degrees:
    print(node_degree)
print(" ")
print("---------------------------------------------------------------------------------------------------")

# Draw tree
Phylo.draw(tree)
