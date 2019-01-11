# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:55:08 2018

@author: Ricardo Chávez Cáliz

Computational experimentation to describe de stochastic process of
rigid expansions.

"""
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from numpy import log
from numpy import mean
from numpy import var
from rigidExpansions import singleExpansion
from rigidExpansions import singleExpansionC
from rigidExpansions import sparseGraphsOptimization
from auxiliar import randomSet
from scipy.special import comb
from auxiliar import intervalo

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter



numExp = 50

def sampleJump(n,p,k,size):
    """
    Input: Graph G (dictionary), set A
    Output: Array
    """
    S = [0]*(size)
    
    if k*log(2) < log(n-k) + (k*p)*log(2):
        method = 1
    else:
        method = 2
        
        
    for i in range(size):
        G=nx.fast_gnp_random_graph(n, p)
        A=randomSet(k,n)
        S[i]= jumpExperiment(G,A,method)
        
    return S


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

def binomial(n,p,k):
    return comb(n,k)*(p**k)*(1-p)**(n-k)

if __name__ == "__main__":
    n=30
    p=0.3
    
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
    plt.show()
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Make data.
    X = np.arange(0, n+1, 1)
    Y = np.arange(0, n+1, 1)
    X, Y = np.meshgrid(X, Y)
    Z = np.asarray(M)
    
    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.set_zlim(0, 1)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()


    """
    k = 50
    sample = sampleJump(n,p,k,1000)
    
    m = mean(sample)
    v = var(sample)
    print(m/v)
    
    n_est = m**2 / (m - v)
    p_est = (m - v)/m
    print(n_est)
    print(p_est)
    
    I= intervalo(min(sample),max(sample)+1,1)
    density = [binomial(n_est,p_est,k) for k in I]
    plt.hist(sample, alpha=0.5, facecolor='#cc0000',edgecolor='#800000', linewidth=1, normed=1)
    plt.plot(I, density, color="#800000",linewidth=1.5, label="densidad")
    plt.title("Histogram")
    plt.show()
    """