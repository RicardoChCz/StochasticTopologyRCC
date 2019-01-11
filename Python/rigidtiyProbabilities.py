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
from bundles import bundleLast
from bundles import bundleAnyone
from bundles import bundleNewBySubset
from bundles import bundleExpansion
from math import floor
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
    
"""
2. Functions to determine if a set of size m can uniquely determine a vertex out-
side of the set
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

"""
3. Functions to determine if a subset B_m of A_k, of sizes m and k, can uniquely 
determine a vertex outside of A_k
"""
def doesSubsetUniquelyDetSomeone(G,A,m):
    """
    Given a graph G and a subset of vertices A, it says if B_m a subset of A_k 
    determine uniquely a vertex v outside of A_k
    Input: Graph G (dictionary), A set, m int
    Output: int (0 or 1) 0 stands for failure and 1 for success 
    """
    B = [list(A)[i] for i in range(0,m)]
    Int = singleVerification(G,set(B))
  
    if len(Int)==1 and next(iter(Int)) not in A:
        return 1
    else:
        return 0
    
def probabilitySubset(n,p,k,m):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k uniquely determines some vertex.
    Input:(int, float, int)
    Output: float
    """
    exitos=0
    if m==0:
        return 0
    
    else:
        for i in range(0,numExperiments):
            G=nx.fast_gnp_random_graph(n, p)
            A=randomSet(k,n)
            exitos = exitos + doesSubsetUniquelyDetSomeone(G,A,m)
            
        return exitos/float(numExperiments)

"""
4. Funciones para determinar si un conjunto de tamaño m puede expandirse
"""    

def doesItExpands(G,A):
    """
    Given a graph G and a subset A, returns succes or failure to the question
    of wheter A expands or doesn't.
    
    Input: Graph G (dictionary), set A 
    Output: int (0 or 1) 0 stands for failure and 1 for success 
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

"""
Funciones para graficar
"""   

def setOfExperiments(N,P,evento):
    """
    Print curves with aproximated probabilities and its cuotes
    Input: list N, list P, int k
    Output: -
    """
    r,c =len(N),len(P)
        
    f, axarr = plt.subplots(r, c,figsize=(9, 6), dpi=100)
    k=l=0
    for n in N:
        I = intervalo(0,n+1,1)
        density = [0]*(len(I))
        bundle1 = [0]*(len(I))

        for p in P:
            for i in I:
                if evento == 1:                    
                    density[i] = probabilityLast(n,p,i)
                    scatterColor = '#0074d9'
                    bundle1[i] = bundleLast(n,p,i)
                    lineColor = '#001f3f'

                elif evento == 2:
                    density[i] = probabilityAnyone(n,p,i)
                    scatterColor = '#FF851B'
                    bundle1[i] = bundleAnyone(n,p,i)
                    lineColor = '#FF4136'

                else:
                    density[i] = probabilityRigidExpansion(n,p,i)
                    scatterColor = '#eb7ab1'
                    bundle1[i] = bundleExpansion(n,p,i)
                    lineColor = '#85144b'
                    
            axarr[k, l].plot(I, density,'o',markersize=5,alpha=0.8,
                 color=scatterColor, label="Experimentation")
            axarr[k, l].plot(I, bundle1,linewidth=1.5,linestyle='-',
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
    Output: -
    """
    I = intervalo(0,n+1,1)
    density= [0]*(len(I))
    bundle1= [0]*(len(I))
    
    for i in I:
        if evento == 1:
            density[i] = probabilityLast(n,p,i)
            scatterColor = '#0074d9'
            bundle1[i] = bundleLast(n,p,i)
            lineColor = '#001f3f'
                    
        elif evento == 2:
            density[i] = probabilityAnyone(n,p,i)
            scatterColor = '#FF851B'
            bundle1[i] = bundleAnyone(n,p,i)
            lineColor = '#FF4136'
        else:
            density[i] = probabilityRigidExpansion(n,p,i)
            scatterColor = '#eb7ab1'
            bundle1[i] = bundleExpansion(n,p,i)
            lineColor = '#85144b'
                
    plt.plot(I, density,'o',markersize=5,alpha=0.8,
                 color=scatterColor, label="Experimentation")
    plt.plot(I, bundle1,linewidth=1.5,linestyle='-',
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
    
def subsetExperiments(n,P):
    """
    Print curves with aproximated probabilities and its bundles
    Input: int n , list P
    Output: -
    """
    r,c = 3,len(P)
    K=[i*floor(n/(r+1)) for i in range(1,r+1) ]
    
    f, axarr = plt.subplots(r, c,figsize=(9, 6), dpi=100)
    i=j=0
    for k in K:
        I = intervalo(1,k+1,1)
        density = [0]*(len(I))
        bundle1 = [0]*(len(I))

        for p in P:
            for m in I:
                density[i] = probabilitySubset(n,p,k,m)
                scatterColor = '#0074d9'
                bundle1[i] = bundleNewBySubset(n,p,k,m)
                lineColor = '#001f3f'
                
            axarr[i, j].plot(I, density,'o',markersize=5,alpha=0.8,
                 color=scatterColor, label="Experimentation")
            axarr[i, j].plot(I, bundle1,linewidth=1.5,linestyle='-',
                 color=lineColor,label="Expected value")
            #Boxes
            if(i==0 and j==2):
                axarr[i, j].legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
                     borderaxespad=0.)
            axarr[i, j].set_title('n='+ str(n) + ', k='+ str(k) + ', p='+ str(p), 
                 backgroundcolor='silver' )
            j+=1
        j=0
        i+=1
    #Pretty
    f.subplots_adjust(hspace=0.4)    
    for i in range(0,r):
        for j in range(0,c):
            axarr[i,j].spines['top'].set_color('white')
            axarr[i,j].spines['right'].set_color('white')
    plt.savefig('Figures/Expansion-by-subset-probability-'+str(n)+'-.png')
    plt.tight_layout()


if __name__ == "__main__":
    """
    #Set of experiments
    # Try with the following functions:
    # 1. Fixed vertex
    # 2. Any vertex
    # 3. By subset
    # 4. Expansion
    """

    """
    N=[5,7,10]
    P=[0.2,0.5,0.8]
    
    individualExperiments(15,0.5,1)
    individualExperiments(15,0.5,2)
    individualExperiments(15,0.5,3)

    setOfExperiments(N,P,1)
    setOfExperiments(N,P,2) 
    setOfExperiments(N,P,3)
    """
    
    n=20
    P=[0.2,0.5,0.8]
    
    subsetExperiments(n,P)