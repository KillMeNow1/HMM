from __future__ import division

import sys
sys.path.append("biopython-1.63/")

import math

from Bio import Phylo
from Bio import AlignIO
from cStringIO import StringIO


import utils

# Defining the letters that we work with and their position

alphabet = "ACGT"
alphabet_index={}
alphabet_index["A"] = 0
alphabet_index["C"] = 1
alphabet_index["G"] = 2
alphabet_index["T"] = 3

if __name__ == "__main__":
	tree = utils.read_tree('data/beta_hemoglobin_tree.nwk')
	alignment = AlignIO.read('data/beta_hemoglobin_short.fasta', "fasta")
	IDs = utils.add_IDs_to_nodes(tree)
	names = utils.tree_nodes_by_name(tree)
	print("The tree is %s" % (tree))
	print("The alingment is %s" % (alignment))
	print(IDs)
	
		#def jukesCantor (ratePair, t):
			#mu=1
			#v = 3*t*mu/4 
			#if a == b:
			#		prob = 0.25 + 0.75*math.exp(-4*v/3)
			#		print("P(%s->%s | mu=1, t=2) = %s" % (a, b, prob))
			#	else:
			#		prob = 0.25 - 0.25*math.exp(-4*v/3)
			#		print("P(%s->%s | mu=1, t=2) = %s" % (a, b, prob))
	
	
	

# running the function here
	
#jukesCantor(ratePair)

#	for i in range(0, len(alignment[0])):
#		column_i =  alignment[:, i]
#		likelihood = "?"
#		print("Col %d: %s, likelihood = %s" % (i, column_i, likelihood))
