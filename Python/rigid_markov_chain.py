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
from rigid_expansion import single_expansion
from rigid_expansion import single_expansion_complement
from rigid_expansion import sparse_graph_optimization
from scipy.special import comb
from _aux import intervalo
from _aux import random_set

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.pyplot import grid
from matplotlib.ticker import LinearLocator, FormatStrFormatter

numExp = 50

def sample_jump(n,p,k,size):
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
        A=random_set(k,n)
        S[i]= jump_experiment(G,A,method)
        
    return S


def vector_success_jump(n,p,k):
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
        A=random_set(k,n)
        V[jump_experiment(G,A,method)] += 1
        
    return [x / numExp for x in V]
        
def jump_experiment(G,A,method=1):
    """
    Given a graph G and a subset of vertices A, returns the size of the set 
    obtained after the first rigid expansion
    Input: Graph G (dictionary), set A
    Output: Int
    """
    #Optimization 1
    A,R = sparse_graph_optimization(G,A)
    
    if method==1:    
        A = single_expansion(G,A)[0]
    else:
        A = single_expansion_complement(G,A)[0]
        
    return len(A.union(R))

def binomial(n,p,k):
    return comb(n,k)*(p**k)*(1-p)**(n-k)


def heatmap(data, ax=None, cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)

    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")


    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)
    
    return im, cbar

if __name__ == "__main__":
    n=30
    p=0.2
    
    M = [0]*(n+1)    
    #It's not possible to expand startig from empty set
    M[0] = [1] + ([0]*(n))
    
    for i in range(1,n):        
        M[i] = vector_success_jump(n,p,i)
        print (i)
        
    #It's already full
    M[n]=([0]*(n))+[1]

    #Plot heat map
    fig, ax = plt.subplots(figsize=(6, 4), dpi=100)
    im, cbar = heatmap(M, ax=ax, cmap="Blues", cbarlabel="Empirical probability")
    fig.tight_layout()
    plt.savefig('Figures/Transition-matrix-secuence-of-rigid-expansions.png')
    plt.show()

    
    
    #3D plotting
    fig = plt.figure(figsize=(6, 4), dpi=100)
    ax = fig.gca(projection='3d')
    plt.rcParams['grid.color'] = "#e6e6e6"
    
    # Make data.
    X = np.arange(0, n+1, 1)
    Y = np.arange(0, n+1, 1)
    X, Y = np.meshgrid(X, Y)
    Z = np.asarray(M)
    
    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap="Blues",
                       linewidth=0, antialiased=False)
    
    # Customize the z axis.
    ax.set_zlim(0, 1)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    plt.savefig('Figures/3D-Transition-matrix-secuence-of-rigid-expansions.png')
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