#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 23:49:09 2018

@author: ricardochcz
"""

import numpy as np

import networkx as nx
import matplotlib.pyplot as plt
from time import clock
from numpy import log
from numpy.random import rand
from numpy.random import randint

def Bernoulli(p):
    """
    Función que simula una v.a distribución Bernoulli con parámetro p
    Input: float
    Output: int (éxito ó fracaso)
    """
    M = rand(1)
    if p > 1 or p<0:
        raise ValueError("El parametro no es apropiado")
    if M < p:
        M=1
    else:
        M=0        
    return M

def erdosRenyi(n,p):
    """
    Devuelve una gráfica aleatoría del modelo Erdos Renyi con parámetros n y p
    Input: int n, float p en [0,1]
    Output: Objeto G de networkX
    """
    #Crea gráfica vacía
    G=nx.Graph()    
    #Añade n nodos
    G.add_nodes_from(list(range(0, n)))
    #Añade aristas según e reultado de v.a. ~ Bernoulli(p)
    for i in range(0,n):
        for j in range(i+1,n):
            if Bernoulli(p)==1:
                G.add_edge(i,j)                            
    return G

def graficaGrafos():
    """
    Grafica un conjunto de gráficas aleatorias Erdos-Renyi variando p, con n=10
    Input: ninguno
    Output: genera un archivo.png
    """
    r=3
    c=3
    k=1.0
    f, axarr = plt.subplots(r, c,figsize=(6, 6.3), dpi=80)
    N = [10, 20, 50]
    A = [0.8, 1, 1.1]
    
    for i in range(0,3):
        for j in range(0,3):
            G = erdosRenyi(N[i], 1/(N[i])**A[j] )
            nx.draw_shell(G,node_size=50, node_color='#003366',ax=axarr[i, j])
            axarr[i, j].set_title('a = ' + str(A[j]) + ", n = " + str(N[i]))
            k=k+1
    plt.savefig('Variables.png')
    f.subplots_adjust(hspace=0.4)
    plt.show()
    

if __name__ == "__main__":
    """
    Muestra de gráficas obtenidas con erdosRenyi
    """
    graficaGrafos()
    