# -*- coding: utf-8 -*-
"""
Created on Tue May 22 11:31:44 2018
@author: Ricardo Chávez Cáliz - RCC

Computational experiments related to rigid expansions
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from __aux import random_set
from __aux import power_set
from __aux import intervalo
from __aux import is_notempty
from __aux_text import write_table
from __rigid_expansion import single_verification
from __rigid_expansion import is_this_uniquely_det
from __bundles import bundle_last
from __bundles import bundle_anyone
from __bundles import bundle_new_by_subset
from __bundles import bundle_expansion
from math import floor
from time import clock
from numpy.random import randint
from scipy.special import comb
from decimal import Decimal

num_experiments = 500
caption_goodness_of_fit = "Supremum of absolute differences between hypothesized and empirical probabilities"

"""
1. Functions to determine whether or not a vertex v is uniquely determinated 
by a random subset of vertices of the graph in G(n,p).
"""
def probability_last(n,p,k):
    """
    Returns the empirical probability that in a G(n,p) i.e. ER-model a 
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
4. Functions to determine if a set of size m can expand
"""    

def does_it_expands(G,A):
    """
    Given a graph G and a subset A, returns success or failure to the question
    of wheter A expands or doesn't.
    
    Input: Graph G (dictionary), set A 
    Output: int (0 or 1) 0 stands for failure and 1 for success 
    """
    P = [x for x in power_set(A)]
    
    #Optmize acording with tresholds with more likely subsets
    
    #Look for leaves and check for expansions
    for i in range (0, len(P)):
        I = single_verification(G,P[i])
        if len(I) == 1 and next(iter(I)) not in A:
            return 1
            break
        
    return 0

def probability_rigid_expansion(n,p,k):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k generates a rigid expansion
    Input: int n, float p, int k
    Output: float
    """
    successes=0   
    total=50
     
    for i in range(0,total):
        G=nx.gnp_random_graph(n, p)
        A=random_set(k,n)
        successes = successes + does_it_expands(G,A)
         
    return successes/float(total)


def set_of_experiments(N,P,event):
    """
    Print curves with aproximated probabilities and its cuotes
    Input: list N, list P, int k
    Output: -
    """
    r,c =len(N),len(P)

    goodness_of_fit = [[0]*c]*r
    
    f, axarr = plt.subplots(r, c,figsize=(9, 6), dpi=100)
    k=l=0
    for n in N:
        I = intervalo(0,n+1,1)
        density = [0]*(len(I))
        bundle = [0]*(len(I))
        for p in P:
            for i in I:                
                if event == 1:                    
                    density[i] = probability_last(n,p,i)
                    bundle[i] = bundle_last(n,p,i)
                    scatter_color = '#0074d9'
                    line_color = '#001f3f'

                elif event == 2:
                    density[i] = probability_anyone(n,p,i)
                    bundle[i] = bundle_anyone(n,p,i)
                    line_color = '#FF4136'
                    scatter_color = '#FF851B'

                else:
                    density[i] = probability_rigid_expansion(n,p,i)
                    bundle[i] = bundle_expansion(n,p,i)
                    line_color = '#85144b'
                    scatter_color = '#eb7ab1'
                    
            goodness_of_fit[N.index(n)][P.index(p)] = '%.2E' % Decimal(max([abs(x - y) for x, y in zip(density, bundle)]))
            
            axarr[k, l].plot(I, density,'o',markersize=5,alpha=0.8,
                 color=scatter_color, label="Experimentation")
            axarr[k, l].plot(I, bundle,linewidth=1.5,linestyle='-',
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

    print("GOF", goodness_of_fit)
    
    if event == 1:
        plt.savefig('Figures/Uniquely-determinated-fixed-vertex.png', bbox_inches="tight")
        write_table("Uniquely-determinated-fixed-vertex-table-errors",caption_goodness_of_fit,"gofExp1",N,P,goodness_of_fit)
        plt.show()
    elif event == 2:
        plt.savefig('Figures/Uniquely-determinated-any-vertex.png', bbox_inches="tight")
        write_table("Uniquely-det-any-table-errors",caption_goodness_of_fit,"gofExp2",N,P,goodness_of_fit)
        plt.show()

    else:
        plt.savefig('Figures/Expansion-probability.png', bbox_inches="tight")
        write_table("Expansion-probability-table-errors",caption_goodness_of_fit,"gofExp3",N,P,goodness_of_fit)
        plt.show()

    plt.tight_layout()

    return goodness_of_fit
    
def individual_experiments(n,p,event):
    """
    Print graphs with aproximated probabilities and its coutes
    Input: list N, list P, int k
    Output: -
    """
    I = intervalo(0,n+1,1)
    density= [0]*(len(I))
    bundle= [0]*(len(I))
    
    for i in I:
        if event == 1:
            density[i] = probability_last(n,p,i)
            bundle[i] = bundle_last(n,p,i)
            scatter_color = '#FF851B'
            line_color = '#001f3f'
        elif event == 2:
            density[i] = probability_anyone(n,p,i)
            scatter_color = '#FF851B'
            bundle[i] = bundle_anyone(n,p,i)
            line_color = '#FF4136'
        else:
            density[i] = probability_rigid_expansion(n,p,i)
            scatter_color = '#eb7ab1'
            bundle[i] = bundle_expansion(n,p,i)
            line_color = '#85144b'

    plt.plot(I, density,'o',markersize=5,alpha=0.8,
                 color=scatter_color, label="Experimentation")
    plt.plot(I, bundle,linewidth=1.5,linestyle='-',
                color=line_color,label="Expected value")

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
                 borderaxespad=0.)

    if event == 1:
        plt.savefig('Figures/Single-uniquely-determinated-fixed-vertex.png', bbox_inches="tight")
    elif event == 2:
        plt.savefig('Figures/Single-uniquely-determinated-any-vertex.png', bbox_inches="tight")
    else:
        plt.savefig('Figures/Single-expansion-probability.png', bbox_inches="tight")
        
    plt.show()

if __name__ == "__main__":
    """
    -------------------------------------------------------------------------
      Computational experiments related to rigid expansions.
    -------------------------------------------------------------------------    
      Event 1: Fixed vertex
      Event 2: Any vertex
      Event 3: Expansion
    -------------------------------------------------------------------------
    """
    n,p=15,0.5
    print("EMPIRICAL VS ESTIMATED PROBABILITIES. n="+ str(n)+ " , p=" + str(p))
    print("Event 1. Prob. of a fixed vertex being uniquely det.")
    individual_experiments(15,0.5,event=1)

    print("Event 2. Prob. of any vertex being uniquely det.")
    individual_experiments(15,0.5,event=2)
    
    print("Event 3. Prob. that a subset generates a rigid expansion.")
    individual_experiments(15,0.5,event=3)

    #------------------------------------------------------------------------
    N,P=[8,15,20],[0.1,0.5,0.75]
    print("EMPIRICAL VS ESTIMATED PROBABILITIES. N="+ str(n)+ " , P=" + str(p))
    print("Event 1. Prob. of a fixed vertex being uniquely det.")
    print(set_of_experiments(N,P,1))

    print("Event 2. Prob. of any vertex being uniquely det.")
    set_of_experiments(N,P,2) 

    print("Event 3. Prob. that a subset generates a rigid expansion.")
    set_of_experiments(N,P,3)