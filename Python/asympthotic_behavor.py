#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 14 12:55:39 2018

@author: Ricardo Esteban Chávez Cáliz

Graphs to understand asympthotic behavor

"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import log 

#------------------------------------------------------------------------------
#Posible values for delta. Functions that tends to infinity.
#------------------------------------------------------------------------------

def deltaLog2(x):
    return log(x/2)


#------------------------------------------------------------------------------
#Definitions for functions
#------------------------------------------------------------------------------    
def conexity(d,I):
    """
    Corresponding values for the conectivty threshold of random Ërdos-Rényi 
    graphs. Whith barely lower and upper cases.
    input: (float, float, array) = (delta, interval)
    output: (array, array, array)
    """
    f = np.zeros(len(I))
    u = np.zeros(len(I))
    l = np.zeros(len(I))
    i = 0
    
    for t in I:
        f[i] = log(t)/t 
        u[i] = f[i] +  d(t)/t
        l[i] = f[i] -  d(t)/t
        i = i+1
    
    return (u,f,l)

def binomialRV(I,a,k):
    """

    input: (float, float, array) = (argument, delta, interval)
    output: (array, array, array)
    """
    f = np.zeros(len(I))
    i = 0
    
    for n in I:
        f[i] = (n**(k*(1-a))) *((1- 1/(n**a))**n)
                
        i = i+1
    
    return f

def T1(I,a,k):
    """
    input: (float, float, array) = (argument, delta, interval)
    output: (array, array, array)
    """
    f = np.zeros(len(I))
    i = 0
    
    for n in I:
        f[i] = (n**(k*(1-a)))
                
        i = i+1
    
    return f

def T2(I,a):
    """
    input: (float, float, array) = (argument, delta, interval)
    output: (array, array, array)
    """
    f = np.zeros(len(I))
    i = 0
    
    for n in I:
        f[i] = ((1- 1/(n**a))**n)
                
        i = i+1
    
    return f
#------------------------------------------------------------------------------
#Methods
#------------------------------------------------------------------------------   
    

#------------------------------------------------------------------------------
if __name__ == "__main__":    
    #Initialize values
    """
    I = np.arange(2, 5000, 10)
    functions = conexity(deltaLog2,I)
    
    #Plot the graphs
    plt.xlabel('Number of vertices')
    plt.ylabel('Value of p')
    
    plt.plot(I, functions[0], linestyle='-.', color='#339933', label="Over threshold")
    plt.plot(I, functions[1], linestyle='-',  color='#0066cc', label="Threshold")
    plt.plot(I, functions[2], linestyle='--', color='#cc0000', label="Below threshold")
             
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.savefig("Cota.png")
    plt.show()
    """
    a = 0.99
    k=10
    #Initialize values
    I = np.arange(2, 500, 10)
    functions = (T2(I,0.2), T2(I,0.5), T2(I,0.8))
    
    #Plot the graphs
    plt.xlabel('Number of vertices')
    plt.ylabel('Asymphtotic behavor of prob of k degree')
    
    plt.plot(I, functions[0], linestyle='-.', color='#339933', label="termino 1")
    plt.plot(I, functions[1], linestyle='-',  color='#0066cc', label="termino 2")
    plt.plot(I, functions[2], linestyle='--', color='#cc0000', label="multiplicación")
             
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.savefig("Cota.png")
    plt.show()