from __future__ import division

import sys
sys.path.append("biopython-1.63/")

from Bio import Phylo
from Bio import AlignIO

import utils

alphabet = "ACGT"
alphabet_index={}
alphabet_index["A"] = 0
alphabet_index["C"] = 1
alphabet_index["G"] = 2
alphabet_index["T"] = 3
	
if __name__ == "__main__":
	tree = utils.read_tree('data/beta_hemoglobin_tree.nwk')
	alignment = AlignIO.read('data/beta_hemoglobin_short.fasta', "fasta")
	
	"""S2.1: Compute all combinations of substitution rate probabilities under 
	the Jukes and Cantor model with mu=1 and t=2"""
	mu = 1
	for a in alphabet:
		for b in alphabet:
			prob = "?"
			print("P(%s->%s | mu=1, t=2) = %s" % (a, b, prob))
	
	print("\n")
	
	"""S3.2: Your first task is to implement Felsentein's algorithm to compute 
	the likelihood of each individual site in the alignment, given the tree and 
	the Jukes and Cantor model with mu = 1"""
	for i in range(0, len(alignment[0])):
		column_i =  alignment[:, i]
		likelihood = "?"
		print("Col %d: %s, likelihood = %s" % (i, column_i, likelihood))

		
	"""S3.4: Now calculate the likelihood of the entire alignment using the 
	same parameters as in S3.2"""

	
	"""S4.3: Your second task is to implement the HMM and compute the likelihood
	of the alignment using the parameters below"""	
	rate_cats = [0.5, 1.0, 2.0]
	freqs = [1.0/len(rate_cats) for rate_cat in rate_cats] # 1/3, 1/3, 1/3
	transition_probs = []  # you are required to implement this calculation
	autocorrelation = 0.75
	
	"""Sanity test: Perform a sanity test of your implementation in D using the parameters
	below. You should get the same likelihood as the one you calculated in B
	(the same model with no rate variation and no autocorrelation)"""
	rate_cats = [1.0, 1.0, 1.0]
	freqs = [1.0/len(rate_cats) for rate_cat in rate_cats] # 1/3, 1/3, 1/3
	transition_probs = [] # you are required to implement this calculation
	autocorrelation = 0

	
	"""S5.1: Your final task is to find the most probable assignment	of
	rates to categories"""
	rate_cats = [0.5, 1.0, 2.0]
	freqs = [1.0/len(rate_cats) for rate_cat in rate_cats] # 1/3, 1/3, 1/3
	transition_probs = []  # you are required to implement this calculation
	autocorrelation = 0.75


