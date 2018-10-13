#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 13:58:03 2018

@author: ricardochcz
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import log


if __name__ == "__main__":
    """
    Obtener gráfica de comportamiento de la cota y los umbrales
    """
    I = np.arange(2,51, 1)
    y1= np.zeros(len(I))
    #y2= np.zeros(len(I))
    #y3= np.zeros(len(I))
    for t in range(0,len(I)):
        d = 1/(I[t])
        y1[t] = (I[t]*log((I[t]**2)/0.02))**(1/d) /(log(I[t]))**3  
        #y2[t] = log(I[t])/I[t] 
        #y3[t] = (log(I[t]) - log(I[t]/2.0))/I[t] 
        
    plt.xlabel('Tamano de la grafica')
    plt.ylabel('Valor de n')
    plt.plot(I, y1, linestyle='-.', color='#339933', label="Sobre Umbral")
    #plt.plot(I, y2, linestyle='-', color='#0066cc', label="Cota")
    #plt.plot(I, y3, linestyle='--', color='#cc0000', label="Bajo Umbral")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("Cota.png")
    plt.show()
    
    #"""
