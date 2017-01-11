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

def jukesCantor (i, j, mu, t):
    v = 3*t*mu/4
    if i == j:
        prob = 0.25 + 0.75*math.exp(-4*v/3)
        return prob
    else:
        prob = 0.25 - 0.25*math.exp(-4*v/3)
        return prob

def TransitionProbs (i, j, lam, freq):
    if i == j:
        prob = lam + (1 - lam)*freq
        return prob
    else:
        prob = (1 - lam)*freq
        return prob

def hmmlikelihood (likelihoods, cache, col, r):
    key = (r, col)
    if key in cache:
        return cache[key]
    elif col == 119:
        return likelihoods[r,col]
    else:
        lik = likelihoods[r, col]
        s = np.NINF
        for k in xrange(0,3):
            s = np.logaddexp(s, np.log(TransitionProbs (r, k, 0.75, 1.0/3.0)) + hmmlikelihood(likelihoods, cache, col+1, k))
        lik = lik + s
        cache[key] = lik
        return lik

def viterbi (likelihoods, cache, cache2, col, r):
    key = (r, col)
    if key in cache2:
        return cache2[key]
    elif col == 119:
        return (-1,likelihoods[r,col])
    else:
        max_index = 0
        temp = np.NINF
        for k in xrange(0,3):
            v = np.log(TransitionProbs (r, k, 0.75, 1.0/3.0)) + viterbi (likelihoods, cache, cache2, col+1, k)[1]
            if v > temp:
                max_index = k
                temp = v
        temp = likelihoods[r,col] + temp
        cache2[key] = (max_index,temp)
        return (max_index,temp)

# Function containing felsensteins algorithm

def Felsenstein (node, a, rate, alignment, cache, seq_indices, col):
    key = (node, a)
    if key in cache:
        return cache[key]
    elif node.is_terminal():
        seq_index = seq_indices[node.name]
        if alignment[seq_index, col] == a:
            return np.log(1)
        else:
            return np.NINF
    else:
        total = np.NINF
        for b in ['A', 'C', 'G', 'T']:
            for c in ['A', 'C', 'G', 'T']:
                trans_prob_left = np.log(jukesCantor (a, b, rate, tree.distance(node, node[0])))
                recursive_left = Felsenstein (node[0], b, rate, alignment, cache, seq_indices, col)
                trans_prob_right = np.log(jukesCantor (a, c, rate, tree.distance(node, node[1])))
                recursive_right = Felsenstein (node[1], c, rate, alignment, cache, seq_indices, col)
                total = np.logaddexp(total, trans_prob_left + recursive_left + trans_prob_right + recursive_right)
        cache[key] = total
        return total

# Obtaining/Reading the necessary files

if __name__ == "__main__":
    tree = utils.read_tree('data/beta_hemoglobin_tree.nwk')
    alignment = AlignIO.read('data/beta_hemoglobin_short.fasta', "fasta")
    IDs = utils.add_IDs_to_nodes(tree)
    names = utils.tree_nodes_by_name(tree)

# Setting up the empty variables

    seq_indices = utils.indices_by_sequence_names(alignment)
    seq_index = ""
    col = 0
    root = tree.get_nonterminals()[0]
    node = root
    overall_log_likelihood = 0
    rate_cats = [0.5, 1.0, 2.0]


# Implementing the algorithm for all bases at the node for all columns
    likelihoods = np.zeros((3, 120))
    for r in xrange(0,3):
        rate = rate_cats[r]
        overall_log_likelihood = 0
        for col in range(0, 120):
            cache = {}
            total = np.NINF
            for a in  ['A', 'C', 'G', 'T']:
                total = np.logaddexp(total, Felsenstein (node, a, rate, alignment, cache, seq_indices, col) + np.log(0.25))
            likelihoods[r, col] = total
            overall_log_likelihood = overall_log_likelihood + total
        #print overall_log_likelihood
    #print likelihoods
    cache = {}
    print hmmlikelihood (likelihoods, cache, 0, 0)
    print hmmlikelihood (likelihoods, cache, 0, 1)
    print hmmlikelihood (likelihoods, cache, 0, 2)

    for r in xrange(0,3):
        cache2 = {}
        for col in xrange(119,-1,-1):
            print col, r, viterbi (likelihoods, cache, cache2, col, r), likelihoods[r,col]

    if ((cache2[0, 0] > cache2 [1, 0]) and (cache2[0, 0] > cache2 [2, 0])):
        print cache2[0,0]
    elif cache2[1,0] > cache2[2,0]:
        print cache2[1,0]
    else:
        print cache2[2,0]


