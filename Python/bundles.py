# -*- coding: utf-8 -*-
"""
Created on Tue Dec 05 11:31:44 2018

@author: Ricardo Chávez Cáliz

Bundles of probability aproximations
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import auxiliar as au

from auxiliar import randomSet
from auxiliar import powerset
from auxiliar import powersetEfective
from auxiliar import intervalo
from auxiliar import is_notempty
from rigidExpansions import singleVerification
from rigidExpansions import isThisUniquelyDet
from time import clock
from numpy.random import randint
from scipy.special import comb

numExperiments = 500

def cotaLast(n,p,m):
    """
    Returns the probaility that in a ER graph G a set of size k uniquely deter-
    minates the n-th vertex.
    """
    if m==n:
        return 0
    else:
        return p**(m) * (1-p**m)**(n-m-1)

def cotaAnyone(n,p,m):
    """
    Returns the probaility that in a ER graph G a set of size k uniquely deter-
    minates the last vertex.
    """
    return 1- (1 - cotaLast(n,p,m))**(n-m)

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
        #a subset of A_k
        rho = [1-cotaExpansionBySubset(n,p,k,m) for m in range(1,k+1)]
        
        #Probability that none of the posible subsets of size m expand outside 
        #of A_k
        exp =[(rho[m-1])**comb(k,m) for m in range(1,k+1)]
        
        prod=1
        
        for e in exp:
            prod = prod*e
 
        return 1-prod

def cotaExpansionBySubset(n,p,k,m):
    """
    Returns the probaility that in a ER graph G a A_m \subset A_k uniquely de-
    terminates the a vertex outside of A_k.
    """
    return 1 - (1- cotaLast(n,p,m))**(n-k) 

def efectiveBundles(n,p,k):
    """
    Returns the lower and upper bundles which determine efective interval to 
    search for rigid expansions
    """
    tooSmall = True
    bigEnough = True
    i=1

    while (tooSmall):
        actual = cotaExpansionBySubset(n,p,k,i)
        if(actual>0.05):
            kMin = i
            i += 1
            tooSmall = False

    while (bigEnough):
        actual = cotaExpansionBySubset(n,p,k,i)
        if(actual<0.05):
            kMax = i
            i += 1
            bigEnough = False

    return (kMin,kMax)