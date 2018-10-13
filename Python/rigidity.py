# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:55:08 2018

@author: Ricardo Chávez Cáliz

Computational experimentation to determinate if a subgraph is rigid
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from aux import randomSet
from aux import powerset
from time import clock
from numpy.random import randint


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
    
    
def f(Y,G):
    """
    Defines a random function f:Y ---> G where Y is a vertex-induced subcomplex
    of G.
    Input: Object from networkx library
    Output: List of pairs
    """
    y=Y.nodes()
    g=G.nodes()
    fun=[]
    for i in y:
        fun.append((i,randint(len(g))))
    return fun
   
def verifyIny1(f):
    """
    Verifies if the function f is inyective.
    Input:
    Output: Boolean
    """
    n=len(f)
    S=set()
    for i in f:
        S.add(i[1])
        
    if n!= len(S):
        return False
    else:
        return True
        

def verifyIny2(f):
    """
    Verifies if the function f is inyective.
    Input: List of pairs
    Output: Boolean
    """
    n=len(f)
    S=set()
    m=0
    r=0
    i=0
    
    while m==r:
        S.add(f[i][1])
        i = i + 1
        m = i + 1
        r = len(S)
          
    if n!= len(S):
        return False
    else:
        return True
  
def verifyLocIny(f,G):
    """
    Verifies if the function f is locally inyective.
    Input:
    Output: Boolean
    """
    y=Y.nodes()
    g=G.nodes()
    fun=[]
    for i in y:
        fun.append((i,randint(len(g))))
    return fun
    
    for i in range(0, len(f)):
        S =  G.neighbors(f[i][0])
        r = []
        for x in S:
            r.append((x,f))
        

if __name__ == "__main__":   
    n=15
    p=0.5
    G=nx.gnp_random_graph(n, p)
    #singleVerification(G, set([0,1,2]))
        
    rigidExpansion(G, randomSet(7,n))
    
 #   verifyLocIny(g,G)
            
    """
    #Compare the two implemented algorithms    
    
    R=200
    T = np.arange(5, R+5, 5)
    y1= np.zeros(len(T))
    y2= np.zeros(len(T))

    j=0
        
    for t in T:
        G=nx.complete_graph(t)
        Y=nx.complete_graph(t)
        g = f(Y,G)
            
        tiempo_inicial = clock()
        verifyIny1(g)
        y1[j] = clock() - tiempo_inicial
        tiempo_inicial = clock()
        verifyIny2(g)
        y2[j] = clock() - tiempo_inicial
        j=j+1
        
    plt.xlabel('Tamano de la grafica')
    plt.ylabel('Tiempo de ejecución')
    plt.plot(T, y1, linestyle='-', color='#ff3300', label="Algortimo 1")
    plt.plot(T, y2, linestyle='-', color='#003366', label="Algortimo 2")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("Cota.png")
    plt.show()
    """
    
    
    
    
    
    
    