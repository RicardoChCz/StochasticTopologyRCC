# -*- coding: utf-8 -*-
"""
Created on Tue May 15 16:58:06 2018
 
@author: rechavez
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
 
from auxiliar import randomSet
from auxiliar import powerset
from auxiliar import intervalo
from time import clock
from numpy.random import randint
from scipy.misc import comb
from math import exp

 
def layoutGraph(G,S):
    """
    Draw the given graph and a subset of it
    Input: Graph G (dictionary), S set
    Output: none
    """
    colores=['#339966']*len(G.nodes())
    #Paint S
    for j in S:    
        colores[j]='#99004d'
         
    plt.figure(num=None, figsize=(3, 3), dpi=80)
    nx.draw_circular(G, node_color=colores)
    plt.show()
 
def singleVerification(G,S):
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
        Neigh = [0]*n
        S = list(S)
        for i in range(0, n):
            Neigh[i] = G.neighbors(S[i])
 
        #Optmizar empezando en el más pequeño y detener cuando sea vacio.
        Int = G.neighbors(S[0])    
        for i in range(0,n):
            Int = set(Int).intersection(Neigh[i])        
     
        return set(Int)
 
def singleExpansion(G,A,R, visual=True):
    """
    Auxiliar method which give a single Expansion in the algorithm rigid
    expansion. Given a graph G and a subset of vertices A, returns the set 
    of vertices that can be uniquely determined using A
    Input: Graph G (dictionary), set A 
    Output: Set
    """   
    if (visual):
        layoutGraph(G,A.union(R))
    P = [x for x in powerset(A)]
    N = set()
     
    for S in P:
        I = singleVerification(G,S)
        if len(I)==1 and next(iter(I)) not in A:
            N = N.union(I)
 
    if len(N)>0:
        A = A.union(N)
        singleExpansion(G,A,R,visual)
 
         
def rigidExpansion(G,A, visual=True):
    """
    Given a graph G and a subset of vertices A, returns the set obtained by a
    sequence of rigid expansions
    Input: Graph G (dictionary), set A 
    Output: Set
    """
    if (visual):
        layoutGraph(G,A)
     
    """
    ---------------------------------------------------------------------------
    Optimization.
    ---------------------------------------------------------------------------
    Hermits (no neighbors). This vertices don't have an effect to rigid 
    expansions, that's why it's posible to remove them.
    Leaves should be removes and petioles (neighbors of leaves) should be added
    to the set to expand 
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
            A=A+N
 
    leaves= set(leaves)
    hermits= set(hermits)
    A=set(A)
    A = A.difference(leaves)
    A = A.difference(hermits)
     
    R=leaves.union(hermits)
    #Call iterative method wich gives multiple single expansions
    singleExpansion(G,A,R,visual)
 
def is_notempty(any_structure):
    if any_structure:
        return True
    else:
        return False
         
def isThisUniquelyDet(G,A,v):
    """
    Given a graph G and a subset of vertices A, it says if the subset A 
    determine the vertex v uniquely
    Input: Graph G (dictionary), A set, v int
    Output: int (0 or 1) 0 stands for failure and 1 for success 
    """
    k=len(A)
    Neigh = [0]*k
    A = list(A)
    for i in range(0, k):
        Neigh[i] = G.neighbors(A[i])
 
    Int = G.neighbors(A[0])
    i=0
     
    while (is_notempty(Int) and i<k):
        Int = set(Int).intersection(Neigh[i])
        i=i+1
     
    Int=list(Int)
     
    if (len(Int)!=1):
        return 0
    elif (v==Int[0]):
        return 1
    else:
        return 0
 
def doesThisUniquelyDetSomeone(G,A):
    """
    Given a graph G and a subset of vertices A, it says if the subset A 
    determine the vertex v uniquely
    Input: Graph G (dictionary), A set, v int
    Output: int (0 or 1) 0 stands for failure and 1 for success 
    """
    k=len(A)
    Neigh = [0]*k
    A = list(A)
    for i in range(0, k):
        Neigh[i] = G.neighbors(A[i])
 
    Int = G.neighbors(A[0])
    i=0
     
    while (is_notempty(Int) and i<k):
        Int = set(Int).intersection(Neigh[i])
        i=i+1
     
    Int=list(Int)
     
    if (len(Int)!=1):
        return 0
    else:
        return 1
         
def probabilityAnyone(n,p,k):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k uniqely determines some vertex.
    Input:(int, float, int)
    Output: float
    """
    exitos=0
    total=1000
     
    for i in range(0,total):
        G=nx.gnp_random_graph(n, p)
        A=randomSet(k,n-1)
        exitos = exitos + doesThisUniquelyDetSomeone(G,A)
         
    return exitos/float(total)
     
def probabilityLast(n,p,k):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k uniqely determines the nth-vertex.
    Input:(int, float, int)
    Output: float
    """
    exitos=0   
    total=1000
     
    for i in range(0,total):
        G=nx.gnp_random_graph(n, p)
        A=randomSet(k,n-1)
        exitos = exitos + isThisUniquelyDet(G,A,n-1)
         
    return exitos/float(total)
 
if __name__ == "__main__":
    """
    #Set of experiments
    N=[10,25,50]
    P=[0.1,0.5,0.9]
 
    r=3
    c=3
    f, axarr = plt.subplots(r, c,figsize=(6, 6), dpi=80)
    k=0
    l=0
    for n in N:
        I = intervalo(1,n,1)
        density= [0]*(len(I))
        cota1= [0]*(len(I))

        for p in P:
            for i in I:
                density[i-1] = probabilityAnyone(n,p,i)     
                cota1[i-1] = (p**(i))*(1-p**i)**(n-i-1)
            axarr[k, l].plot(I, density,linewidth=2,color='#339933')
            axarr[k, l].plot(I, cota1,linewidth=2,color='#ff6600')
            axarr[k, l].set_title('n='+ str(n) + ' p='+ str(p) )
            l=l+1
        l=0
        k=k+1
     
    f.subplots_adjust(hspace=0.35)
    plt.savefig('curvas.png')
    plt.tight_layout()
    
    #Experimentos individuales

    n=30
    p=0.5
    I = intervalo(1,n,1)
    density= [0]*(len(I))
    cota1= [0]*(len(I))

    for i in I:
        density[i-1] = probabilityLast(n,p,i)
        cota1[i-1] = (p**(i))*(1-p**i)**(n-i-1)
    
    plt.plot(I, density,linewidth=2,color='#339933')
    plt.plot(I, cota1,linewidth=2,color='#ff6600')
    plt.show()
    
    #"""
    
    I = intervalo(-0.05,1,0.5)
    f1= [0]*(len(I))
    f2= [0]*(len(I))
    f3= [0]*(len(I))

    for i in range(0, len(I)):
        x=I[i]
        d=3
        f1[i] = (d+1)*(x+1)*exp(-x) + x*(1-exp(-x))**(d+1)
        d=5
        f2[i] = (d+1)*(x+1)*exp(-x) + x*(1-exp(-x))**(d+1)
        d=10
        f3[i] = (d+1)*(x+1)*exp(-x) + x*(1-exp(-x))**(d+1)

    
    plt.plot(I, f1,linewidth=2,color='red')
    plt.plot(I, f2,linewidth=2,color='blue')
    plt.plot(I, f3,linewidth=2,color='green')
    plt.show()
