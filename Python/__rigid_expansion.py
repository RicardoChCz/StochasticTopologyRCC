# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:55:08 2018
@author: Ricardo Chávez Cáliz - RCC

Computational experimentation to generate rigid expansions of a subgraph
"""
import imageio
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from __aux import random_set
from __aux import power_set
from __aux import power_set_efective
from __aux_graphic import layout_graph
from __aux_graphic import visual_rigid_exp

from numpy.random import randint
from numpy import log
from time import clock

def single_verification(G,S):
    """
    Given a graph G and a subset of vertices S, it says if the subset determine
    a vertex uniquely
    Input: Graph G (dictionary), S set
    Output: Set
    """
    if len(S)==0:
        return set()
    else:        
        n=len(S)
        neigh = [0]*n
        S = list(S)
        for i in range(0, n):
            neigh[i] = G.neighbors(S[i])

        #Optmizar empezando en el más pequeño y detener cuando sea vacio.
        intersection = G.neighbors(S[0])    
        for i in range(0,n):
            intersection = set(intersection).intersection(neigh[i])        
    
        return set(intersection)

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

def iterative_expansion(G,A,R, visual=True, method=1):
    """
    Auxiliar method which give a single Expansion in the algorithm rigid
    expansion. Given a graph G and a subset of vertices A, returns the set 
    of vertices that can be uniquely determined using A
    Input: Graph G (dictionary), set A 
    Output: Set
    """
    i=0
    if(visual):
        filenames = []
        print("We begin with A = "+ str(A))
        visual_rigid_exp(G,A.union(R),filenames,i)

    if method==1:    
        expand = single_expansion
    else:
        expand = single_expansion_complement

    A,are_there_new_ones = expand(G,A)

    while are_there_new_ones:
        if(visual):
            i+=1
            print("Expansion "+ str(i) + ", with A="+ str(A))
            visual_rigid_exp(G,A.union(R),filenames,i)

        A,are_there_new_ones = expand(G,A)

    #Generate gif
    if(visual):
        print("I'm done doing expansions!")
        images = []
        for filename in filenames:
            images.append(imageio.imread(filename))
            imageio.mimsave('Figures/rigid_expansion_gif/rigid_expansion.gif', images, duration=1)

    return A.union(R)

def sparse_graph_optimization(G,A):
    """
    ---------------------------------------------------------------------------
    Optimization.
    ---------------------------------------------------------------------------
    Hermits (no neighbors). This vertices don't have an effect to rigid 
    expansions, that's why it's posible to remove them.
    Leaves should be removed and petioles (neighbors of leaves) should be added
    to the set to expand.
    
    ---------------------------------------------------------------------------
    Given a graph G and a subset of vertices A, returns two parts:
    1. The usefull set A' without leves and isolated points but with petioles
    added.
    2. Leaves and isolated points.
    Input: Graph G (dictionary), set A 
    Output: (set, set)
    """
    
    A=list(A)
    hermits = []
    leaves = []
    for u in A:
        N = list(G.neighbors(u))
        if len(N) == 0:
            hermits.append(u)
        elif len(N) == 1:
            leaves.append(u)
            #Adding petioles
            A=A+N

    leaves= set(leaves)
    hermits= set(hermits)
    A=set(A)
    A = A.difference(leaves)
    A = A.difference(hermits)
    R=leaves.union(hermits)
    
    return(A, R)

def single_expansion(G,A):
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

def single_expansion_complement(G,A):
    """
    Auxiliar method which give a single Expansion in the algorithm rigid
    expansion. Given a graph G and a subset of vertices A, returns the set 
    of vertices that can be uniquely determined using A, this algorithm review
    which of the vertices in the complement can be u.d. using A
    Input: Graph G (dictionary), set A 
    Output: Set
    """
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

def rigid_expansion(G,A,p,visual=True):
    """
    Given a graph G and a subset of vertices A, returns the set obtained by a
    sequence of rigid expansions
    Input: Graph G (dictionary), set A 
    Output: Set
    """
    #Optimization 1
    A,R = sparse_graph_optimization(G,A)   
    
    #Optimization 2
    n,k=len(G.nodes),len(A)
    if k*log(2) < log(n-k) + (k*p)*log(2):
        method = 1
    else:
        method = 2

    #Call iterative method wich gives multiple single expansions
    return (iterative_expansion(G,A,R,visual,method))