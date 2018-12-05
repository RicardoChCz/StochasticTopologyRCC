# -*- coding: utf-8 -*-
"""
Created on Tue May 22 11:31:44 2018

@author: Ricardo Chávez Cáliz

Computational experimentation to verify probabilistic aproximations
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import auxiliar as au

from auxiliar import randomSet
from auxiliar import powerset
from auxiliar import intervalo
from auxiliar import is_notempty
from rigidExpansions import singleVerification
from rigidExpansions import isThisUniquelyDet
from time import clock
from numpy.random import randint
from scipy.special import comb

numExperiments = 500

"""
1. Functions to determine whether or not a vertex v is uniquely determinated by
a random subset of vertices of the graph in G(n,p).
"""

def probabilityLast(n,p,k):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k uniqely determines the nth-vertex.
    Input:(int, float, int)
    Output: float
    """
    exitos=0   
    total=numExperiments
    if k==n:
        return 0
    
    else:
        for i in range(0,total):
            G=nx.gnp_random_graph(n, p)
            A=randomSet(k,n-1)
            exitos = exitos + isThisUniquelyDet(G,A,n-1)
             
        return exitos/float(total)
    
def cotaLast(n,p,m):
    """
    Returns the probaility that in a ER graph G a set of size k uniquely deter-
    minates the n-th vertex.
    """
    if m==n:
        return 0
    else:
        return p**(m) * (1-p**m)**(n-m-1)
    
"""
2. Funciones para determinar si un conjunto de tamaño m puede determinar de 
manera única un vértice fuera del conjunto
"""
def doesThisUniquelyDetSomeone(G,A):
    """
    Given a graph G and a subset of vertices A, it says if the subset A 
    determine uniquely a vertex v outside of A
    Input: Graph G (dictionary), A set, v int
    Output: int (0 or 1) 0 stands for failure and 1 for success 
    """
    Int=list(singleVerification(G,A))
    
    if len(Int)==1 and next(iter(Int)) not in A:
        return 1
    else:
        return 0
    
def probabilityAnyone(n,p,k):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k uniqely determines some vertex.
    Input:(int, float, int)
    Output: float
    """
    exitos=0
    total=numExperiments
     
    for i in range(0,total):
        G=nx.gnp_random_graph(n, p)
        A=randomSet(k,n)
        exitos = exitos + doesThisUniquelyDetSomeone(G,A)
         
    return exitos/float(total)

def cotaAnyone(n,p,m):
    """
    Returns the probaility that in a ER graph G a set of size k uniquely deter-
    minates the last vertex.
    """
    return 1- (1 - cotaLast(n,p,m))**(n-m)

"""
3. Funciones para determinar si un conjunto de tamaño m puede expandirse
"""    

def doesItExpands(G,A):
    """
    Given a graph G and a subset A, returns succes or failure to the question
    of wheter A expands or doesn't.
    
    Input: Graph G (dictionary), set A 
    Output: int (0 or 1)
    """
    P = [x for x in powerset(A)]
    
    #Hacer optimiación de acuerdo donde es más posible encontrar expansiones
    
    #Buscar si hay hojas y ahí preguntarte por expansiones
    for i in range (0, len(P)):
        I = singleVerification(G,P[i])
        if len(I) == 1 and next(iter(I)) not in A:
            return 1
            break
        
    return 0

def probabilityRigidExpansion(n,p,k):
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

def cotaExpansion(n,p,k):
    """
    Returns lower and upper bound for the probaility that in a ER graph G a set
    of size k generates a rigid expansion
    Input: int n, float p, int k
    Output: float
    """
    if k==0:
        return 0
    else:
        #Probability that none of the vertex outside of A_k is u.d by B_m 
        # a subset of A_k
        rho = [(1-cotaLast(n,p,m))**(n-k) for m in range(1,k+1)]
        
        #Probability that none of the posible subsets of size m expand outside 
        #of A_k
        exp =[(rho[m-1])**comb(k,m) for m in range(1,k+1)]
        
        prod=1
        
        for e in exp:
            prod = prod*e
 
        return 1-prod

"""
Funciones para graficar
"""   

def setOfExperiments(N,P,evento):
    """
    Print curves with aproximated probabilities and its cuotes
    Input: list N, list P, int k
    Output: float
    """
    r=len(N)
    c=len(P)
        
    f, axarr = plt.subplots(r, c,figsize=(9, 6), dpi=100)
    k=0
    l=0
    for n in N:
        I = intervalo(0,n+1,1)
        
        density=[0]*(len(I))
        cota1= [0]*(len(I))

        for p in P:
            for i in I:
                if evento == 1:                    
                    density[i] = probabilityLast(n,p,i)
                    scatterColor = '#0074d9'
                    cota1[i] = cotaLast(n,p,i)
                    lineColor = '#001f3f'

                elif evento == 2:
                    density[i] = probabilityAnyone(n,p,i)
                    scatterColor = '#FF851B'
                    cota1[i] = cotaAnyone(n,p,i)
                    lineColor = '#FF4136'

                else:
                    density[i] = probabilityRigidExpansion(n,p,i)
                    scatterColor = '#eb7ab1'
                    cota1[i] = cotaExpansion(n,p,i)
                    lineColor = '#85144b'
                    
            axarr[k, l].plot(I, density,'o',markersize=5,alpha=0.8,
                 color=scatterColor, label="Experimentation")
            axarr[k, l].plot(I, cota1,linewidth=1.5,linestyle='-',
                 color=lineColor,label="Expected value")
            
            #Boxes
            if(k==0 and l==2):
                axarr[k, l].legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
                     borderaxespad=0.)
                                
            axarr[k, l].set_title('n='+ str(n) + ', p='+ str(p), 
                 backgroundcolor='silver' )
            l=l+1
        l=0
        k=k+1
     
    #Pretty
    f.subplots_adjust(hspace=0.4)    
    for i in range(0,r):
        for j in range(0,c):
            axarr[i,j].spines['top'].set_color('white')
            axarr[i,j].spines['right'].set_color('white')

    if evento == 1:
        plt.savefig('Figures/Uniquely-determinated-fixed-vertex.png')
    elif evento == 2:
        plt.savefig('Figures/Uniquely-determinated-any-vertex.png')
    else:
        plt.savefig('Figures/Expansion-probability.png')

    plt.tight_layout()
    
def individualExperiments(n,p,evento):
    """
    Print graphs with aproximated probabilities and its coutes
    Input: list N, list P, int k
    Output: float
    """
    I = intervalo(0,n+1,1)
    density= [0]*(len(I))
    cota1= [0]*(len(I))
    
    for i in I:
        if evento == 1:
            density[i] = probabilityLast(n,p,i)
            scatterColor = '#0074d9'
            cota1[i] = cotaLast(n,p,i)
            lineColor = '#001f3f'
                    
        elif evento == 2:
            density[i] = probabilityAnyone(n,p,i)
            scatterColor = '#FF851B'
            cota1[i] = cotaAnyone(n,p,i)
            lineColor = '#FF4136'
        else:
            density[i] = probabilityRigidExpansion(n,p,i)
            scatterColor = '#eb7ab1'
            cota1[i] = cotaExpansion(n,p,i)
            lineColor = '#85144b'
                
    plt.plot(I, density,'o',markersize=5,alpha=0.8,
                 color=scatterColor, label="Experimentation")
    plt.plot(I, cota1,linewidth=1.5,linestyle='-',
                color=lineColor,label="Expected value")

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
                 borderaxespad=0.)
            
    if evento == 1:
        plt.savefig('Figures/Single-uniquely-determinated-fixed-vertex.png')
    elif evento == 2:
        plt.savefig('Figures/Single-uniquely-determinated-any-vertex.png')
    else:
        plt.savefig('Figures/Single-expansion-probability.png')
        
    plt.show()
    
    
if __name__ == "__main__":
    """
    #Set of experiments
    # Try with the following functions:
    # 1. Fixed vertex
    # 2. Any vertex
    # 3. Expansion
    """
    N=[5,7,10]
    P=[0.2,0.5,0.8]
    individualExperiments(15,0.5,1)
    individualExperiments(15,0.5,2)
    individualExperiments(15,0.5,3)

    setOfExperiments(N,P,1)
    setOfExperiments(N,P,2) 
    setOfExperiments(N,P,3)
