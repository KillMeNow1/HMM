# Need to calculate the algorithm with logarithms - more accurate results
# N.B. Python uses log as the syntax for ln (natural log)

from __future__ import division

import sys
sys.path.append("biopython-1.63/")

import math

from Bio import Phylo
from Bio import AlignIO
from cStringIO import StringIO

import numpy as np

import utils


# Calculating the jukes and cantor probabilities of base mutation

def jukesCantor (i, j, t):
    mu=1
    v = 3*t*mu/4 
    if i == j:
        prob = 0.25 + 0.75*math.exp(-4*v/3)
        return prob
    else:
        prob = 0.25 - 0.25*math.exp(-4*v/3)
        return prob
        
# Function containing felsensteins algorithm

def Felsenstein (node, a, alignment, cache, seq_indices, col):  
    key = (node, a)                                
    if key in cache:
        return cache[key]
    elif node.is_terminal():
        seq_index = seq_indices[node.name]
        if alignment[seq_index, col] == a:
            return np.log(1)
        else:
            return np.NINF          # Python doesn't like log(0) so need to use -ve infinity instead
    else:
        total = np.NINF
        for b in ['A', 'C', 'G', 'T']:
            for c in ['A', 'C', 'G', 'T']:
                trans_prob_left = np.log(jukesCantor (a, b, tree.distance(node, node[0])))
                recursive_left = Felsenstein (node[0], b, alignment, cache, seq_indices, col)
                # This does not need to be logged as the result of felsenstein is already a log              
                trans_prob_right = np.log(jukesCantor (a, c, tree.distance(node, node[1])))
                recursive_right = Felsenstein (node[1], c, alignment, cache, seq_indices, col)
                total = np.logaddexp(total, trans_prob_left + recursive_left + trans_prob_right + recursive_right)
                # This function can only accept two variables, so perform the sum in the function
        cache[key] = total        
        return total

# Obtaining/Reading the necessary files

if __name__ == "__main__":
    tree = utils.read_tree('data/beta_hemoglobin_tree.nwk')
    alignment = AlignIO.read('data/beta_hemoglobin_short.fasta', "fasta")
    names = utils.tree_nodes_by_name(tree)

# Setting up the empty variables
        
    seq_indices = utils.indices_by_sequence_names(alignment)
    seq_index = ""
    cache = {}
    col = 0
    root = tree.get_nonterminals()[0]
    node = root   
    total = np.NINF

# Implementing the algorithm for all bases at the node
    
    for a in  ['A', 'C', 'G', 'T']:
        total = np.logaddexp(total, Felsenstein (node, a, alignment, cache, seq_indices, col) + np.log(0.25))
      
    print total

