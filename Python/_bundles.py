# -*- coding: utf-8 -*-
"""
Created on Tue Dec 05 11:31:44 2018

@author: Ricardo Chávez Cáliz

Bundles of probability aproximations
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import _aux as au

from _aux import random_set
from _aux import power_set
from _aux import power_set_efective
from _aux import intervalo
from _aux import is_notempty
from rigid_expansion import single_verification
from rigid_expansion import is_this_uniquely_det
from time import clock
from numpy.random import randint
from scipy.special import comb

num_experiments = 500

def bundle_last(n,p,m):
    """
    Returns the probaility that in a ER graph G a set of size k uniquely deter-
    minates the n-th vertex.
    """
    if m==n:
        return 0
    else:
        return p**(m) * (1-p**m)**(n-m-1)

def bundle_anyone(n,p,m):
    """
    Returns the probaility that in a ER graph G a set of size k uniquely deter-
    minates the last vertex.
    """
    return 1- (1 - bundle_last(n,p,m))**(n-m)

def bundle_new_by_subset(n,p,k,m):
    """
    Returns the probaility that in a ER graph G a A_m \subset A_k uniquely de-
    terminates the a vertex outside of A_k.
    """
    if m==0:
        return 0
    else:
        return 1 - (1- ((n-k)/n)*bundle_last(n,p,m))**(n-k) 

def bundle_expansion(n,p,k):
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
        rho = [1-bundle_new_by_subset(n,p,k,m) for m in range(1,k+1)]
        
        #Probability that none of the posible subsets of size m expand outside 
        #of A_k
        exp =[(rho[m-1])**comb(k,m) for m in range(1,k+1)]
        
        prod=1
        
        for e in exp:
            prod = prod*e
 
        return 1-prod

def efective_bundles(n,p,k):
    """
    Returns the lower and upper bundles which determine efective interval to 
    search for rigid expansions
    """
    too_small = True
    big_enough = True
    i=1

    while (too_small):
        actual = bundle_new_by_subset(n,p,k,i)
        if(actual>0.05):
            k_min = i
            i += 1
            too_small = False

    while (big_enough):
        actual = bundle_new_by_subset(n,p,k,i)
        if(actual<0.05):
            k_max = i
            i += 1
            big_enough = False

    return (k_min,k_max)