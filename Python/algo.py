# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:55:08 2018

@author: Ricardo Chávez Cáliz
"""
import numpy as np

import networkx as nx
import matplotlib.pyplot as plt

def powerset(seq):
    """
    Returns all the subsets of this set.
    """
    if len(seq) <= 1:
        yield seq
        yield []
    else:
        for item in powerset(seq[1:]):
            yield [seq[0]]+item
            yield item
            
def revisa(G):
    """
    Dado v vertice de G revisa si hay un conjunto que lo determine de manera única
    """
    Unis=[]
    for v in G.nodes():
        A=G.neighbors(v)
        if len(uniDet(G,A))==1:
            Unis.append(v)
    return Unis
    
def expansion(G,A):
    dibuja(G,A)
    P = [x for x in powerset(A)]
    N = []
    
    for S in P:
        veci = uniDet(G,S)
        if len(veci)==1 and veci[0] not in A:
            N = list(set(N).union([veci[0]]))

    if len(N)>0:
        A = list(set(A).union(N))
        expansion(G,A)
    
def uniDet(G,A):
    """
    Given a graph G and a subset of vertices A, returns the 
    
    
    
    set of vertices 
    that can be uniquely determined using A.
    Input: graph G (dictionary netoworkx), list A 
    Output: list
    """
    n=len(A)
    Vec = [0]*n
    for i in np.xrange(0, n):
        Vec[i] = G.neighbors(A[i])
    
    Int = G.nodes()
    for i in np.xrange(0,n):
        Int = set(Int).intersection(Vec[i])
        
    return list (Int)
    
def dibuja(G,A):
    inte = uniDet(G,A)
    colores=['#75a3a3']*len(G.nodes())
    #colorea la intersección
    for j in inte:
        colores[j]='#ff6600'
    #colorea A
    for j in A:    
        colores[j]='#990033'   
    nx.draw_circular(G, node_color=colores)
    plt.show()

if __name__ == "__main__":
    n=15
    p=0.25
    G=nx.gnp_random_graph(n, p)
    U = revisa(G)
    
    for u in U:
        dibuja(G,G.neighbors(u))

    expansion(G,[0,1])    
    
    
    