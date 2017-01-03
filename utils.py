from Bio import Phylo
from Bio import AlignIO
import math
import numpy

def read_tree(tree_file):
	"""Reads a newick tree file and returns a biopython Tree object"""
	tree = Phylo.read(tree_file, 'newick')
	tree.root_at_midpoint()
	return tree

def tree_nodes_by_name(tree):
	"""Returns a dictionary of tree nodes that can be indexed by their respective names"""	
	names = {}
	for clade in tree.find_clades():
		if clade.name:
			if clade.name in names:
				raise ValueError("Duplicate key: %s" % clade.name)
			names[clade.name] = clade
	return names

def sequence_names_by_index(alignment):
	seq_names = []
	for seq in alignment:
		seq_names.append(seq.name)
	return seq_names
		
def indices_by_sequence_names(alignment):
	indices = {}
	i = 0
	for seq in alignment:
		indices[seq.name] = i
		i += 1
	return indices

def add_IDs_to_nodes(tree):
    names = {}
    n_clades = 0;
    for idx, clade in enumerate(tree.find_clades()):
        clade.ID = idx;
        names[idx] = clade.name
        if n_clades < idx:
            n_clades = idx;
    return names, n_clades + 1

if __name__ == "__main__":
	print("-----------------------------------------------")
	print("Loading a tree:")
	tree = read_tree('../data/beta_hemoglobin_tree.nwk')
        names = add_IDs_to_nodes(tree)

	print(tree)
        print(names)

	print("\n-----------------------------------------------")
	print("Iterating over immediate child nodes of a node (the parent):")
	root_node = tree.get_nonterminals()[0] # get root node
        print root_node[1][1].is_terminal()
	for child_node in root_node:
		print child_node
	
	print("\n-----------------------------------------------")
	print("Calculate the total branch lengths between any two nodes:")
	root_node = tree.get_nonterminals()[0] # get root node
	for child_node in root_node:
		print("t(%s, %s) = %f" % (root_node, child_node, tree.distance(root_node, child_node)))
		
	print("\n-----------------------------------------------")
	print("Getting a tree node by its name:")
	nodes_by_name = tree_nodes_by_name(tree)
	didelphis_node = nodes_by_name["Didelphis"]
	print(didelphis_node)

	print("\n-----------------------------------------------")
	print("Testing whether a node is a leaf node or a non-terminal node:")
	print("Is %s a leaf node? %s" % (didelphis_node, didelphis_node.is_terminal()))
	print("Is the root node a leaf node? %s" % (root_node.is_terminal()))

	print("\n-----------------------------------------------")
	print("Loading an alignment:")
	alignment = AlignIO.read('data/beta_hemoglobin_short.fasta', "fasta")
	print(alignment)

	print("\n-----------------------------------------------")
	print("Obtaining a column (a site) from the alignment by column index:")
	print("Column %d: %s" % (0, alignment[:, 0]))
	
	print("\n-----------------------------------------------")
	print("Retrieving the name of a node:")
	print("name=%s" % didelphis_node.name)
	
	print("\n-----------------------------------------------")
	print("Obtaining the nucleotide character assigned to a particular node by sequence name:")
	sequence_name = "Homo_sapiens"
	seq_indices = indices_by_sequence_names(alignment)
	seq_index = seq_indices[sequence_name]
	for i in range(0, 5):
		column = alignment[:, i]
		print("Nucleotide for %s in column %d: %s" % (sequence_name, i, column[seq_index]))
	

	print("\n-----------------------------------------------")
	print("Calculating the logarithm of the sum of exponentials (requires numpy):")
	print("log(0.25 + 0.75) = log[e^log(0.25) + e^log(0.75)] = %f" % numpy.logaddexp(math.log(0.25), math.log(0.75)))
	print("(this is useful for summing likelihoods in a log-space in a numerically stable way)")
