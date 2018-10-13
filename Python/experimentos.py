# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 11:14:30 2018

@author: rechavez
"""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import collections

from numpy import log
from numpy.random import rand
from numpy.random import randint
from networkx.generators.degree_seq import expected_degree_graph

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
    for i in xrange(0,n):
        for j in xrange(i+1,n):
            if Bernoulli(p)==1:
                G.add_edge(i,j)                            
    return G

def conexidad(n,p,M):
    """
    Estima la probabilidad de obtener una gráfica aleatoria conexa en G(n,p)
    Input: int n, float p entre 0 y 1, int M - número de repeticiones
    Output: float (probabilidad estimada)
    """
    r=0.0
    for i in xrange(0,M):
        G=erdosRenyi(n,p)
        if nx.is_connected(G)==True:
            r=r+1
    return r/M

def mapeaHistogramas(n,r):
    K = np.arange(1,n*r,r)        

    for k in K:
        s = np.random.poisson(k, 10000)
        count, bins, ignored = plt.hist(s, 14, normed=True, color = (0.2,float(k)/K[len(K)-1],0.8), edgecolor=(0.2,float(k)/K[len(K)-1],0.8), alpha=0.8)
    plt.show()


if __name__ == "__main__":
    mapeaHistogramas(10,2)

    s = np.random.poisson(10000, 10000)
    count, bins, ignored = plt.hist(s, 14, normed=True, color = '#339933', alpha=0.8)
    plt.show()
    
    print "Gráfica de comportamiento de las cota para conexidad y grado en el modelo Erdös-Rényi"
    I = np.arange(2000, 2100, 1)
    y1= np.zeros(len(I))
    y2= np.zeros(len(I))
    for t in xrange(0,len(I)):
        y1[t] = (2.0) /I[t]  
        y2[t] = log(I[t])/I[t] 
        
    plt.xlabel('Tamano de la grafica')
    plt.ylabel('Valor de p')
    plt.plot(I, y1, linestyle='-', color='#ff2d2d', label="Grado")
    plt.plot(I, y2, linestyle='-', color='#0066cc', label="Conexidad")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()

    print "Grafica conexa"
    n=50
    p=0.9
    
    G = nx.gnp_random_graph(n, p)
    plt.figure(num=None, figsize=(5, 5), dpi=80)
    nx.draw_shell(G,node_size=10, node_color='#003366', width=0.1 )
    plt.show()
    print nx.is_connected(G)
    
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    # print "Degree sequence", degree_sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())

    fig, ax = plt.subplots()
    plt.bar(deg, cnt, width=0.80, color='#003366', edgecolor='#003366', alpha=0.8)

    plt.title("Degree Histogram")
    plt.ylabel("Count")
    plt.xlabel("Degree")
    ax.set_xticks([d + 0.4 for d in deg])
    ax.set_xticklabels(deg)
    
    # draw graph in inset
    plt.axes([0.4, 0.4, 0.5, 0.5])
    Gcc = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)[0]
    pos = nx.spring_layout(G)
    plt.axis('off')
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos, alpha=0.4)
    
    plt.show()
    