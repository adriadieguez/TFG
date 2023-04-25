From newick file to graph:
	`python3 newick_to_graph.py exemple_arbre.txt`

From FASTA to frequecies:
	`python3 align_to_freq.py Alignment.fasta`

`alignments_4_leaves` contains files with alignments for 4-leaves-trees inside each file. Some of them have length 1000 and others 10000. File **a15_b101** may have some issues to consider (i.e. division by 0...), so some problems could arise. In case we want alignments for 3-leaves-trees, just erase one of them.

https://github.com/TomMakesThings/Phylogenetic-Tree-Inference/blob/main/Phylogenetics.ipynb

-----------------------HOW TO GENERATE ALIGNMENTS?------------------------------------------------

Option 1: Generate t FASTA files of sequences of length L given a tree in a newick format:

	For example, t = 5 and L = 1000;  `python3 GenGM.py tree_4L.txt 5 1000` 

Option 2: Generate FASTA files of sequences of given lengths L1...Ld. **Note that in this case, lenghts sequences are preceeded by an L**

	For example, L1 = 500, L2 = 1000 and L3 = 10000;  `python3 GenGM.py tree_4L.txt L500 L1000 L10000` 
