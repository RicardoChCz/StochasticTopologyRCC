
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 13:58:03 2018

@author: ricardochcz
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import log
from numpy import exp

def f(n):
    return log(n)

def g(n):
    return n/2

def h(n):
    return 2*log(n) 

def cota(x,vx):
    k=8
    return x * (vx)**k * exp(-vx)

if __name__ == "__main__":
    """
    Obtener gráfica de comportamiento de la cota y los umbrales
    """
    k=5
    I = np.arange(2,5100000, 1000)
    y1= np.zeros(len(I))
    y2= np.zeros(len(I))
    y3= np.zeros(len(I))
    for t in range(0,len(I)):
        y1[t] = cota(I[t], f(I[t]))
        y2[t] = cota(I[t], g(I[t]))
        y3[t] = cota(I[t], h(I[t]))
        
    plt.xlabel('Tamano de la grafica')
    plt.ylabel('Valor de n')
    #plt.plot(I, y1, linestyle='-.', color='#339933', label="f(n)=log(n)/n")
    plt.plot(I, y2, linestyle='-', color='#0066cc', label="g(n)=cte")
    plt.plot(I, y3, linestyle='--', color='#cc0000', label="h(n)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("Cota.png")
    plt.show()
    
    #"""
