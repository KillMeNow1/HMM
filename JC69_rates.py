from __future__ import division

import sys
sys.path.append("biopython-1.63/")

import math

from Bio import Phylo
from Bio import AlignIO

import utils

alphabet = "ACGT"
alphabet_index={}
alphabet_index["A"] = 0
alphabet_index["C"] = 1
alphabet_index["G"] = 2
alphabet_index["T"] = 3

ratePair = list(raw_input("Please enter a basepair: "))

"""S2.1: Compute all combinations of substitution rate probabilities under 
	the Jukes and Cantor model with mu=1 and t=2"""
def jukesCantor (ratePair):
	t=2
	mu=1
	v = 3*t*mu/4 
	if len(ratePair) != 2:
		print "Please enter a pair of bases!"
	else:
		a = ratePair[0]
		b = ratePair[1]
		if (a in alphabet and b in alphabet):
			if a == b:
				prob = 0.25 + 0.75*math.exp(-4*v/3)
				print("P(%s->%s | mu=1, t=2) = %s" % (a, b, prob))
			else:
				prob = 0.25 - 0.25*math.exp(-4*v/3)
				print("P(%s->%s | mu=1, t=2) = %s" % (a, b, prob))
		else:
			print "Please enter nucleotides basepair only."
	
	

	
	
jukesCantor(ratePair)

	
