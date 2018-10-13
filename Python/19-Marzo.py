# -*- coding: utf-8 -*-
"""
Created on Sat May 19 15:30:06 2018

@author: rechavez
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from auxiliar import randomSet
from auxiliar import powerset
from auxiliar import intervalo
from auxiliar import is_notempty
from time import clock
from numpy.random import randint
from scipy.misc import comb


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
            Neigh[i] = list(G.neighbors(S[i]))
        
        Neigh=sorted(Neigh, key=len)        
        #Optmizar empezando en el más pequeño y detener cuando sea vacio.
        Int = set(G.neighbors(S[0]))
        j=1
        while (is_notempty(Int) and j<len(Neigh)):
            Int = set(Int).intersection(Neigh[j])
            j=j+1
            
        return set(Int)

def doesItExpands(G,A):
    """
    Given a graph G and a subset A, returns succes or failure to the question
    of wheter A expands or doesn't.
    
    Input: Graph G (dictionary), set A 
    Output: int (0 or 1)
    """
    P = [x for x in powerset(A)]
    Y = 0
    i = 0
    
    #Hacer optimiación de acuerdo  donde es más posible encontrar expansiones
    #Buscar si hay hojas y ahpi preguntarte por expansiones
    while Y==0:
        S=P[i] 
        I = singleVerification(G,S)
        if len(I)==1 and next(iter(I)) not in A:
            Y=1
        else:
            i=i+1
            if i==len(P):
                break
    return Y

def probabiltyRigidExpansion(n,p,k):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k generates a rigid expantion
    Input: int n, float p, int k
    Output: float
    """
    exitos=0   
    total=100
     
    for i in range(0,total):
        G=nx.gnp_random_graph(n, p)
        A=randomSet(k,n)
        exitos = exitos + doesItExpands(G,A)
         
    return exitos/float(total)

def cotaExp(n,p,k):
    """
    Returns lower and upper bound for the probaility that in a ER graph G a set
    of size k generates a rigid expansion
    """
    M=range(0,k+1)
 
    print M   
    for m in M:
        M[m] = comb(k, m, exact=True) * (p**m) * (1-p)**(n-m-1)
    print M
    
    return (max(M), sum(M)) 
    
if __name__ == "__main__":
    """    
    #Set of experiments
    N=[6,10]
    P=[0.1,0.5,0.9]
 
    r=len(N)
    c=len(P)
    
    f, axarr = plt.subplots(r, c,figsize=(6, 6), dpi=80)
    k=0
    l=0
    for n in N:
        I = intervalo(1,n,1)
        density= [0]*(len(I))
        cota1= [0]*(len(I))

        for p in P:
            for i in I:
                density[i-1] = probabiltyRigidExpansion(n,p,i)   
            axarr[k, l].plot(I, density,linewidth=2,color='#339933')
            axarr[k, l].plot(I, cota1,linewidth=2,color='#ff6600')
            axarr[k, l].set_title('n='+ str(n) + ' p='+ str(p) )
            l=l+1
        l=0
        k=k+1
     
    f.subplots_adjust(hspace=0.35)
    plt.savefig('curvas.png')
    plt.tight_layout()
    """
    #Individual experiments
    n=10
    p=0.3
    I = intervalo(1,n,1)
    density= [0]*(len(I))
    cotaI= [0]*(len(I))
    cotaS= [0]*(len(I))
    

    for i in I:
        density[i-1] =  probabiltyRigidExpansion(n,p,i)
        A = cotaExp(n,p,i)
        cotaI[i-1] = A[0]
        cotaS[i-1] = A[1]

    plt.plot(I, density,linewidth=2,color='#660066', label="Experimentation")
    plt.plot(I, cotaI,linewidth=2,color='#339933', label="Lower bound")
    plt.plot(I, cotaS,linewidth=2,color='#ff6600', label="Upper bound")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)    
    plt.show()
    
    