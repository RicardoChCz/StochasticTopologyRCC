#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 17:32:12 2019

@author: ricardochcz
"""

import matplotlib.pyplot as plt
from auxiliar import intervalo
from math import log
import numpy as np

"""
1. Functions to determine whether or not a vertex v is uniquely determinated by
a random subset of vertices of the graph in G(n,p).
"""

def probabilityDiameter(n,c):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k uniqely determines the nth-vertex.
    Input:(int, float, int)
    Output: float
    """

    return ((n*log(n**2/c))**(1/2.0)) / n


def graphCurves():
    """
    Print graphs with aproximated probabilities and its coutes
    Input: list N, list P, int k
    Output: -
    """
    I = intervalo(2,1000,10)
    curve = [0]*(len(I))
    lineColor="red"
    for c in [x/100.0 for x in range(1,10)]:
        print(c)
        for i in range(0, len(I)):
            curve[i] = probabilityDiameter(I[i],c)
            plt.plot(I, curve,linewidth=1,linestyle='-',color=lineColor,label="c= 0.5")

    plt.show()
    

if __name__ == "__main__":
    graphCurves()