From newick file to graph:
	`python3 newick_to_graph.py exemple_arbre.txt`

From FASTA to frequecies:
	`python3 align_to_freq.py Alignment.fasta`

`alignments_4_leaves` contains files with alignments for 4-leaves-trees inside each file. Some of them have length 1000 and others 10000. File **a15_b101** may have some issues to consider (i.e. division by 0...), so some problems could arise. In case we want alignments for 3-leaves-trees, just erase one of them.

https://github.com/TomMakesThings/Phylogenetic-Tree-Inference/blob/main/Phylogenetics.ipynb
