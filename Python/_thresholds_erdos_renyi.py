#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 17:32:12 2019
@author: Ricardo Chávez Cáliz - RCC

Proability thresholds for Erdös-Renyi model
"""

import matplotlib.pyplot as plt
from __aux import intervalo
from math import log

"""
1. Lower bundles for properties in ER model
"""

def bundle_connectivity(n):
    """
    Returns the lower bundle for diameter treshold
    Input: n
    Output: p(n)
    """
    return log(n) / n

def bundle_diameter(n):
    """
    Returns the lower bundle for diameter treshold
    Input: n
    Output: p(n)
    """
    return 1.0 / n**(0.5)

def bundle_locally_infinite(n):
    """
    Returns the lower bundle for locally infinite treshold
    Input: n
    Output: p(n)
    """
    return 1.0 / n

def graph_curves():
    """
    Print bundles
    Input: -
    Output: -
    """
    I = intervalo(5,100,1)
    curve = [0]*(len(I))
    properties = [
      {
        'bundle': bundle_locally_infinite,
        'color': '#ebc844',
        'label': 'Locally infinite'
      },
      {
        'bundle': bundle_connectivity,
        'color': '#f16c20',
        'label': 'Connected'
      },
      {
         'bundle': bundle_diameter,
         'color': '#0d3c55',
         'label': 'Diameter 2'
      }
    ]
    for prop in properties:
        for i in range(0,I.size):
            curve[i]=prop['bundle'](I[i])
            
        plt.plot(I, curve,linewidth=1.5,linestyle='-',
                color=prop['color'],label=prop['label'])
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
                 borderaxespad=0.)
        plt.fill_between(I, curve, 0.5, facecolor=prop['color'], alpha=0.4)
    
    plt.savefig('Figures/ER-thresholds.png')
    plt.show()
if __name__ == "__main__":
    graph_curves()