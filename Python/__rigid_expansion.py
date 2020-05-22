# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:55:08 2018
@author: Ricardo Chávez Cáliz - RCC

Computational experimentation to generate rigid expansions of a subgraph
"""
import imageio
import networkx as nx

from __aux import power_set
from __aux_graphic import visual_rigid_exp
from numpy import log

def single_verification(G,S):
    """
    Given a graph G and a subset of vertices S, it says if the subset determine
    a vertex uniquely
    Input: Graph G (dictionary), S set
    Output: Set
    """
    i = 0
    S = list(S)
    if len(S) == 0:
        return set()
    else:
        intersection = set(G.neighbors(S[i]))
        while len(intersection) != 0 and i < len(S)-1:
            i += 1
            intersection = set(intersection).intersection(G.neighbors(S[i]))
        return intersection

def is_this_uniquely_det(G,A,v):
    """
    Given a graph G and a subset of vertices A, it says if the subset A 
    determine the vertex v uniquely
    Input: Graph G (dictionary), A set, v int
    Output: int (0 or 1) 0 stands for failure and 1 for success 
    """
    Int = list(single_verification(G,A))
     
    if (len(Int)!=1):
        return 0
    elif (v==Int[0]):
        return 1
    else:
        return 0

def iterative_expansion(G, A, G_r, A_r, p, visual=True, optimized=True):
    """
    Auxiliar method which give a single Expansion in the algorithm rigid
    expansion. Given a graph G and a subset of vertices A, returns the set 
    of vertices that can be uniquely determined using A
    Input: Graph G (dictionary), set A to expand, set R for optmization, 
           float p as parameter, boolean options
    Output: Set, expanded A
    """
    i = 0
    n = len(G.nodes)
    k = len(A)

    if(visual):
        filenames = []
        print("We begin with A = "+ str(A.union(A_r)))
        visual_rigid_exp(G, A, G_r, A_r, filenames,i)
  
    if k == 0:
        return A.union(A_r)

    #Optimization 2
    if optimized and k*log(2) > (log(n-k) + (k*p)*log(2)):
        faster_with_complement = True
        expand = single_expansion_complement
    else:
        faster_with_complement = False
        expand = single_expansion

    #Iterative
    A,are_there_new_ones = expand(G, A, A_r)
    while are_there_new_ones:
        # Change expansion method if convinient. Optimization 2
        if optimized and not faster_with_complement:
            k = len(A)
            if k*log(2) > (log(n-k) + (k*p)*log(2)):
                faster_with_complement = True
                expand = single_expansion_complement
  
        if(visual):
            i+=1
            print("Expansion "+ str(i) + ", with A="+ str(A.union(A_r)))
            visual_rigid_exp(G, A, G_r, A_r, filenames,i)

        A,are_there_new_ones = expand(G, A, A_r, optimized)

    #Generate gif
    if(visual):
        print("Generating gif... ")
        images = []
        for filename in filenames:
            images.append(imageio.imread(filename))
            imageio.mimsave('Figures/rigid_expansion_gif/rigid_expansion.gif', images, duration=1)

    return A.union(A_r)

def sparse_graph_optimization(G,A):
    """
    ---------------------------------------------------------------------------
    Optimization.
    ---------------------------------------------------------------------------
    Isolated vertices. This vertices don't have an effect to rigid expansions.
    A-Leaves. Leaves in A, should be removed and petioles (neighbors of leaves) 
    should be added A.
    ---------------------------------------------------------------------------
    Given a graph G and a subset of vertices A, returns:
    1. G' - Graph without isolated vertices and A-leaves.
    2. A' - The usefull set A; without leaves or isolated points but with 
    petioles added.
    3. R - Removed vertices ie G = G' \cup R.
    Input: Graph G (dictionary), set A
    Output: (G', A', G_removed, A_removed) - (dictionary, set, set)
    """

    isolates = set(v for v in nx.isolates(G))
    A_isolates = isolates.intersection(A)
    G_isolates = isolates.difference(A_isolates)
    
    G.remove_nodes_from(isolates)

    A = A.difference(A_isolates)
    A_leaves = []
    A_petioles = []
    for u in A:
        if G.degree[u] == 1:
            A_leaves.append(u)
            uniq_neigh = list(G.neighbors(u))[0]
            if G.degree[uniq_neigh] == 1:
                A_leaves.append(uniq_neigh)
            else:
                A_petioles.append(uniq_neigh)
                
            
    A = A.difference(set(A_leaves))
    A = A.union(set(A_petioles))

    return(G, A, G_isolates, A_isolates.union(A_leaves))

def single_expansion(G, A, A_removed, optimized=True):
    """
    Auxiliar method which give a single Expansion in the algorithm rigid
    expansion. Given a graph G and a subset of vertices A, returns the set 
    of vertices that can be uniquely determined using A
    Input: Graph G (dictionary), set A 
    Output: Set
    """
    N = set()

    for S in power_set(A):
        I = single_verification(G,S)
        if len(I)==1 and next(iter(I)) not in A:
            N = N.union(I)
 
    if len(N)>0:
        are_there_new_ones = True
        A = A.union(N)

    else:
        are_there_new_ones = False
        
    return (A,are_there_new_ones)

def single_expansion_complement(G, A, A_removed, optimized=True):
    """
    Auxiliar method which give a single Expansion in the algorithm rigid
    expansion. Given a graph G and a subset of vertices A, returns the set 
    of vertices that can be uniquely determined using A. This algorithm review
    which of the vertices in the complement can be u.d. using A
    Input: Graph G (dictionary), set A 
    Output: Set
    """
    A = A.union(A_removed)
    B = set(G.nodes) - A
    N = set()

    for b in B:
        AN = A.intersection(G.neighbors(b))
    
        for S in power_set(AN):
            if is_this_uniquely_det(G,S,b) == 1:
                N = N.union(set([b]))
                break
 
    if len(N)>0:
        are_there_new_ones = True
        A = A.union(N)

    else:
        are_there_new_ones = False
        
    return (A,are_there_new_ones)

def rigid_expansion(G,A,p,visual=True, optimized=True):
    """
    Given a graph G and a subset of vertices A, returns the set obtained by a
    sequence of rigid expansions
    Input: Graph G (dictionary), set A 
    Output: Set
    """
    #Optimization 1
    if optimized:
        G,A,G_r,A_r = sparse_graph_optimization(G,A)
    else:
        G_r, A_r = set(), set()
    
    #Call iterative method wich gives multiple single expansions
    return (iterative_expansion(G,A,G_r,A_r,p,visual,optimized))