# -*- coding: utf-8 -*-
"""
Created on Tue May 22 11:31:44 2018

@author: Ricardo Chávez Cáliz

Computational experimentation to verify probabilistic aproximations
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import _aux as au

from _aux import random_set
from _aux import power_set
from _aux import intervalo
from _aux import is_notempty
from rigid_expansion import single_verification
from rigid_expansion import is_this_uniquely_det
from _bundles import bundle_last
from _bundles import bundle_anyone
from _bundles import bundle_new_by_subset
from _bundles import bundle_expansion
from math import floor
from time import clock
from numpy.random import randint
from scipy.special import comb

num_experiments = 500

"""
1. Functions to determine whether or not a vertex v is uniquely determinated by
a random subset of vertices of the graph in G(n,p).
"""
def probability_last(n,p,k):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k uniqely determines the nth-vertex.
    Input:(int, float, int)
    Output: float
    """
    exitos=0
    total=num_experiments
    if k==n:
        return 0
    
    else:
        for i in range(0,total):
            G=nx.gnp_random_graph(n, p)
            A=random_set(k,n-1)
            exitos = exitos + is_this_uniquely_det(G,A,n-1)
             
        return exitos/float(total)    
    
"""
2. Functions to determine if a set of size m can uniquely determine a vertex out-
side of the set
"""
def does_this_uniquely_det_someone(G,A):
    """
    Given a graph G and a subset of vertices A, it says if the subset A 
    determine uniquely a vertex v outside of A
    Input: Graph G (dictionary), A set, v int
    Output: int (0 or 1) 0 stands for failure and 1 for success 
    """
    Int=list(single_verification(G,A))
    
    if len(Int)==1 and next(iter(Int)) not in A:
        return 1
    else:
        return 0
    
def probability_anyone(n,p,k):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k uniqely determines some vertex.
    Input:(int, float, int)
    Output: float
    """
    exitos=0
    total=num_experiments
     
    for i in range(0,total):
        G=nx.gnp_random_graph(n, p)
        A=random_set(k,n)
        exitos = exitos + does_this_uniquely_det_someone(G,A)
         
    return exitos/float(total)

"""
3. Functions to determine if a subset B_m of A_k, of sizes m and k, can uniquely 
determine a vertex outside of A_k
"""
def does_subset_uniquely_det_someone(G,A,m):
    """
    Given a graph G and a subset of vertices A, it says if B_m a subset of A_k 
    determine uniquely a vertex v outside of A_k
    Input: Graph G (dictionary), A set, m int
    Output: int (0 or 1) 0 stands for failure and 1 for success 
    """
    B = [list(A)[i] for i in range(0,m)]
    Int = single_verification(G,set(B))
  
    if len(Int)==1 and next(iter(Int)) not in A:
        return 1
    else:
        return 0
    
def probability_subset(n,p,k,m):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k uniquely determines some vertex.
    Input:(int, float, int)
    Output: float
    """
    successes=0
    if m==0:
        return 0
    
    else:
        for i in range(0,num_experiments):
            G=nx.fast_gnp_random_graph(n, p)
            A=random_set(k,n)
            successes = successes + does_subset_uniquely_det_someone(G,A,m)
            
        return successes/float(num_experiments)

"""
4. Funciones para determinar si un conjunto de tamaño m puede expandirse
"""    

def does_it_expands(G,A):
    """
    Given a graph G and a subset A, returns success or failure to the question
    of wheter A expands or doesn't.
    
    Input: Graph G (dictionary), set A 
    Output: int (0 or 1) 0 stands for failure and 1 for success 
    """
    P = [x for x in power_set(A)]
    
    #Hacer optimiación de acuerdo donde es más posible encontrar expansiones
    
    #Buscar si hay hojas y ahí preguntarte por expansiones
    for i in range (0, len(P)):
        I = single_verification(G,P[i])
        if len(I) == 1 and next(iter(I)) not in A:
            return 1
            break
        
    return 0

def probability_rigid_expansion(n,p,k):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k generates a rigid expantion
    Input: int n, float p, int k
    Output: float
    """
    successes=0   
    total=100
     
    for i in range(0,total):
        G=nx.gnp_random_graph(n, p)
        A=random_set(k,n)
        successes = successes + does_it_expands(G,A)
         
    return successes/float(total)

"""
Funciones para graficar
"""   

def set_of_experiments(N,P,event):
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
                if event == 1:                    
                    density[i] = probability_last(n,p,i)
                    scatter_color = '#0074d9'
                    bundle1[i] = bundle_last(n,p,i)
                    line_color = '#001f3f'

                elif event == 2:
                    density[i] = probability_anyone(n,p,i)
                    scatter_color = '#FF851B'
                    bundle1[i] = bundle_anyone(n,p,i)
                    line_color = '#FF4136'

                else:
                    density[i] = probability_rigid_expansion(n,p,i)
                    scatter_color = '#eb7ab1'
                    bundle1[i] = bundle_expansion(n,p,i)
                    line_color = '#85144b'
                    
            axarr[k, l].plot(I, density,'o',markersize=5,alpha=0.8,
                 color=scatter_color, label="Experimentation")
            axarr[k, l].plot(I, bundle1,linewidth=1.5,linestyle='-',
                 color=line_color,label="Expected value")
            
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

    if event == 1:
        plt.savefig('Figures/Uniquely-determinated-fixed-vertex.png')
    elif event == 2:
        plt.savefig('Figures/Uniquely-determinated-any-vertex.png')
    else:
        plt.savefig('Figures/Expansion-probability.png')

    plt.tight_layout()
    
def individual_experiments(n,p,event):
    """
    Print graphs with aproximated probabilities and its coutes
    Input: list N, list P, int k
    Output: -
    """
    I = intervalo(0,n+1,1)
    density= [0]*(len(I))
    bundle1= [0]*(len(I))

    
    for i in I:
        if event == 1:
            density[i] = probability_last(n,p,i)
            bundle1[i] = bundle_last(n,p,i)
            line_color = '#001f3f'
                    
        elif event == 2:
            density[i] = probability_anyone(n,p,i)
            scatter_color = '#FF851B'
            bundle1[i] = bundle_anyone(n,p,i)
            line_color = '#FF4136'
        else:
            density[i] = probability_rigid_expansion(n,p,i)
            scatter_color = '#eb7ab1'
            bundle1[i] = bundle_expansion(n,p,i)
            line_color = '#85144b'
                
    plt.plot(I, density,'o',markersize=5,alpha=0.8,
                 color=scatter_color, label="Experimentation")
    plt.plot(I, bundle1,linewidth=1.5,linestyle='-',
                color=line_color,label="Expected value")

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
                 borderaxespad=0.)
            
    if event == 1:
        plt.savefig('Figures/Single-uniquely-determinated-fixed-vertex.png')
    elif event == 2:
        plt.savefig('Figures/Single-uniquely-determinated-any-vertex.png')
    else:
        plt.savefig('Figures/Single-expansion-probability.png')
        
    plt.show()
    
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
                density[i] = probability_subset(n,p,k,m)
                scatter_color = '#0074d9'
                bundle1[i] = bundle_new_by_subset(n,p,k,m)
                line_color = '#001f3f'
                
            axarr[i, j].plot(I, density,'o',markersize=5,alpha=0.8,
                 color=scatter_color, label="Experimentation")
            axarr[i, j].plot(I, bundle1,linewidth=1.5,linestyle='-',
                 color=line_color,label="Expected value")
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

    N=[5,7,10]
    P=[0.2,0.5,0.8]
    
    individual_experiments(15,0.5,1)
    individual_experiments(15,0.5,2)
    individual_experiments(15,0.5,3)

    set_of_experiments(N,P,1)
    set_of_experiments(N,P,2) 
    set_of_experiments(N,P,3)
    """
    
    n=20
    P=[0.2,0.5,0.8]
    
    subsetExperiments(n,P)
    """