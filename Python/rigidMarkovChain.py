# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:55:08 2018

@author: Ricardo Chávez Cáliz

Computational experimentation to describe de stochastic process of
rigid expansions.

"""
import matplotlib.pyplot as plt
import networkx as nx

from numpy import log
from rigidExpansions import singleExpansion
from rigidExpansions import singleExpansionC
from rigidExpansions import sparseGraphsOptimization
from auxiliar import randomSet

numExp = 50

def vectorSuccessJump(n,p,k):
    """
    Input: Graph G (dictionary), set A
    Output: Int
    """
    V = [0]*(n+1)
    
    if k*log(2) < log(n-k) + (k*p)*log(2):
        method = 1
        
    else:
        method = 2
        
        
    for i in range(numExp):
        G=nx.fast_gnp_random_graph(n, p)
        A=randomSet(k,n)
        V[jumpExperiment(G,A,method)] += 1
        
    return [x / numExp for x in V]
        
def jumpExperiment(G,A,method=1):
    """
    Given a graph G and a subset of vertices A, returns the size of the set 
    obtained after the first rigid expansion
    Input: Graph G (dictionary), set A
    Output: Int
    """
    #Optimization 1
    A,R = sparseGraphsOptimization(G,A)
    
    if method==1:    
        A = singleExpansion(G,A)[0]
    else:
        A = singleExpansionC(G,A)[0]
        
    return len(A.union(R))

if __name__ == "__main__":
    n=50
    p=0.1
    
    fig, ax = plt.subplots()
    M = [0]*(n+1)    
    #It's not possible to expand startig from empty set
    M[0] = [1] + ([0]*(n))
    
    for i in range(1,n):        
        M[i] = vectorSuccessJump(n,p,i)
        print (i)
        
    #It's already full
    M[n]=([0]*(n))+[1]
    
    ax.matshow(M, cmap=plt.cm.Blues)
    """
    for i in range(n+1):
        for j in range(n+1):
            c = M[j][i]
            if c>0:
                ax.text(i, j, str(c), va='center', ha='center')        
    """