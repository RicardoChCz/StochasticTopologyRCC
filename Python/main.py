#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 23:30:28 2019

@author: ricardochcz
"""
import networkx as nx
from _aux import random_set
from rigid_expansion import rigid_expansion

if __name__ == "__main__":   
    n=15
    p=0.5
    G=nx.gnp_random_graph(n, p)
    A=random_set(7,n)
    print("A = ", A)
    print("expanded A = ", rigid_expansion(G, A))
    