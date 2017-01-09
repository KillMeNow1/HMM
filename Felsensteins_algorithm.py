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
    key = (node, a)                                     # creating the dictionary entry
    if key in cache:                                    # checking if entry has been calculated already
        return cache[key]
    elif node.is_terminal():                            # formula for leaf-nodes
        seq_index = seq_indices[node.name]
        if alignment[seq_index, col] == a:
            return 1
        else:
            return 0
    else:                                               # formula for ancestral-nodes
        total = 0
        for b in ['A', 'C', 'G', 'T']:
            for c in ['A', 'C', 'G', 'T']:
                trans_prob_left = jukesCantor (a, b, tree.distance(node, node[0]))              
                # time variable taken from length of branch in the tree
                recursive_left = Felsenstein (node[0], b, alignment, cache, seq_indices, col)   
                # calling the function within the function - recurssion: node[0] = lefthand child
                trans_prob_right = jukesCantor (a, c, tree.distance(node, node[1]))
                # node[1] = righthand child
                recursive_right = Felsenstein (node[1], c, alignment, cache, seq_indices, col)
                total += trans_prob_left * recursive_left * trans_prob_right * recursive_right
        cache[key] = total                              # saving the result to the dictionary - prevent repeated calculations 
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
    root = tree.get_nonterminals()[0]       # selects the root node
    node = root   
    total = 0

# Implementing the algorithm for all bases at the node
    
    for a in  ['A', 'C', 'G', 'T']:
        total += Felsenstein (node, a, alignment, cache, seq_indices, col) * 0.25
      
    print total

