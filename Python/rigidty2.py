# -*- coding: utf-8 -*-
"""
Created on Tue May 22 11:31:44 2018

@author: Ricardo Chávez Cáliz

Computational experimentation to determinate if a subgraph is rigid
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import auxiliar as au

from auxiliar import randomSet
from auxiliar import powerset
from auxiliar import intervalo
from auxiliar import is_notempty
from auxiliar import singleVerification
from time import clock
from numpy.random import randint
from scipy.misc import comb


"""
1. Funciones para determinar si un vertice v está determinado de manera única por
un subconjunto en una gráfica aleatoria Erdös-Renyi
"""

def isThisUniquelyDet(G,A,v):
    """
    Given a graph G and a subset of vertices A, it says if the subset A 
    determine the vertex v uniquely
    Input: Graph G (dictionary), A set, v int
    Output: int (0 or 1) 0 stands for failure and 1 for success 
    """
    Int = list(singleVerification(G,A))
     
    if (len(Int)!=1):
        return 0
    elif (v==Int[0]):
        return 1
    else:
        return 0

def probabilityLast(n,p,k):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k uniqely determines the nth-vertex.
    Input:(int, float, int)
    Output: float
    """
    exitos=0   
    total=1000
     
    for i in range(0,total):
        G=nx.gnp_random_graph(n, p)
        A=randomSet(k,n-1)
        exitos = exitos + isThisUniquelyDet(G,A,n-1)
         
    return exitos/float(total)
    
def cotaLast(n,p,m):
    """
    Returns the probaility that in a ER graph G a set of size k uniquely deter-
    minates the last vertex.
    """
    return p**(m) * (1-p**m)**(n-m-1)
    
"""
2. Funciones para determinar si un conjunto de tamaño m puede determinar de 
manera única un vértice fuera del conjunto
"""
def doesThisUniquelyDetSomeone(G,A):
    """
    Given a graph G and a subset of vertices A, it says if the subset A 
    determine uniquely a vertex v outside of A
    Input: Graph G (dictionary), A set, v int
    Output: int (0 or 1) 0 stands for failure and 1 for success 
    """
    Int=list(singleVerification(G,A))
    
    if len(Int)==1 and next(iter(Int)) not in A:
        return 1
    else:
        return 0
    
def probabilityAnyone(n,p,k):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k uniqely determines some vertex.
    Input:(int, float, int)
    Output: float
    """
    exitos=0
    total=1000
     
    for i in range(0,total):
        G=nx.gnp_random_graph(n, p)
        A=randomSet(k,n)
        exitos = exitos + doesThisUniquelyDetSomeone(G,A)
         
    return exitos/float(total)

def cotaAnyone(n,p,m):
    """
    Returns the probaility that in a ER graph G a set of size k uniquely deter-
    minates the last vertex.
    """
    return (n-m) * p**(m) * (1-p**m)**(n-m-1)

"""
3. Funciones para determinar si un conjunto de tamaño m puede expandirse
"""    

def doesItExpands(G,A):
    """
    Given a graph G and a subset A, returns succes or failure to the question
    of wheter A expands or doesn't.
    
    Input: Graph G (dictionary), set A 
    Output: int (0 or 1)
    """
    P = [x for x in powerset(A)]
    Y = 0
    i = 0
    
    #Hacer optimiación de acuerdo  donde es más posible encontrar expansiones
    
    #Buscar si hay hojas y ahí preguntarte por expansiones
    while Y==0:
        S=P[i] 
        I = singleVerification(G,S)
        if len(I)==1 and next(iter(I)) not in A:
            Y=1
        else:
            i=i+1
            if i==len(P):
                break
    return Y

def probabilityRigidExpansion(n,p,k):
    """
    Returns the probability (aproximated) that in a G(n,p) i.e. ER-model a 
    random set of size k generates a rigid expantion
    Input: int n, float p, int k
    Output: float
    """
    exitos=0   
    total=100
     
    for i in range(0,total):
        G=nx.gnp_random_graph(n, p)
        A=randomSet(k,n)
        exitos = exitos + doesItExpands(G,A)
         
    return exitos/float(total)

def cotaExpansion(n,p,k):
    """
    Returns lower and upper bound for the probaility that in a ER graph G a set
    of size k generates a rigid expansion
    Input: int n, float p, int k
    Output: float
    """
    
    E = [comb(k, m, exact=True) *(n-k)* p**m * (1-p**m)**(n-m-1) for m in range(1,k+1)]
    
    rho = [(n-k)*p**m * (1-p**m)**(n-m-1) for m in range(1,k+1)]
    
    ant = [(1-rho[m-1])**comb(k, m, exact=True) for m in range (1,k)]
    a=1
    for t in ant:
        a=a*t

    return (sum(E), max(rho), 1 - a)
    
def setOfExperiments(N,P,evento):
    """
    Print graphs with aproximated probabilities and its coutes
    Input: list N, list P, int k
    Output: float
    """
    r=len(N)
    c=len(P)
    
    f, axarr = plt.subplots(r, c,figsize=(6, 6), dpi=80)
    k=0
    l=0
    for n in N:
        I = intervalo(1,n,1)
        density=[0]*(len(I))
        cota1= [0]*(len(I))
        if evento==3:
            cota2=[0]*(len(I))

        for p in P:
            for i in I:
                if evento == 1:
                    density[i-1] = probabilityLast(n,p,i)
                    cota1[i-1] = cotaLast(n,p,i)
                    
                elif evento == 2:
                    density[i-1] = probabilityAnyone(n,p,i)
                    cota1[i] = cotaAnyone(n,p,i)
                else:
                    density[i-1] = probabilityRigidExpansion(n,p,i)
                    #C = cotaExpansion(n,p,i)
                    cota1[i-1] = cotaExpansion(n,p,i)
                    #cota2[i-1] = C[1]
                    
            axarr[k, l].plot(I, density,linewidth=2,color='#660066', label="Experimentation")
            axarr[k, l].plot(I, cota1,linewidth=2,color='#339933',label="Cota")
            if evento==3:
                #axarr[k, l].plot(I, cota2,linewidth=2,color='#ff6600', label="Cota sup")
                print("Holi")
                
            axarr[k, l].set_title('n='+ str(n) + ' p='+ str(p) )
            l=l+1
        l=0
        k=k+1
     
    f.subplots_adjust(hspace=0.35)
    plt.tight_layout()    
    
def individualExperiments(n,p,evento):
    """
    Print graphs with aproximated probabilities and its coutes
    Input: list N, list P, int k
    Output: float
    """
    I = intervalo(1,n+1,1)
    density= [0]*(len(I))
    cota1= [0]*(len(I))
    if evento==3:
            cota2=[0]*(len(I))
    
    for i in I:
        if evento == 1:
            density[i-1] = probabilityLast(n,p,i)
            cota1[i-1] = cotaLast(n,p,i)
                    
        elif evento == 2:
            density[i-1] = probabilityAnyone(n,p,i)
            cota1[i-1] = cotaAnyone(n,p,i-1)
        else:
            density[i-1] = probabilityRigidExpansion(n,p,i)
            C = cotaExpansion(n,p,i)
            cota1[i-1] = C[0]
            cota2[i-1] = C[1]
            plt.plot(I, cota2,linewidth=2,color='#ff6600', label="Upper bound")

    if evento==2:
        print(np.argmax(cota1))
                
    plt.plot(I, density,linewidth=2,color='#660066', label="Experimentation")
    plt.plot(I, cota1,linewidth=2,color='#339933', label="Calculada")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)    
    plt.show()
    
    
if __name__ == "__main__":
    """
    #Set of experiments
    # Try with the following functions:
    # 1.Ultimo
    # 2.Cualquiera
    # 3.Expansion
    """
    N=[10,20,30]
    P=[0.1,0.5,0.9]    
    #setOfExperiments(N,P,2)
    """
    r=len(N)
    c=len(P)
    
    f, axarr = plt.subplots(r, c,figsize=(6, 6), dpi=80)
    k=0
    l=0
    for n in N:
        I = intervalo(1,n,1)
        funcion=[0]*(len(I))
        
        for p in P:
            for i in I:
                funcion=au.g(n,p,i)
                
            axarr[k, l].plot(I, funcion,linewidth=2,color='#660066', label="Experimentation")                
            axarr[k, l].set_title('n='+ str(n) + ' p='+ str(p) )
            au.newton(n,p,0.001)
            l=l+1
        l=0
        k=k+1
     
    f.subplots_adjust(hspace=0.35)
    plt.tight_layout()
    """
    n=15
    p=0.7
    
    I = intervalo(1,n+1,1)    
    funcion=[0]*(len(I))
       
    for k in I:
        funcion[k-1] = au.g(n,p,k)
    
    plt.plot(I, funcion,linewidth=2,color='#660066', label="funcion a max")
    plt.show()
    au.newton(n,p,0.001)
    
    individualExperiments(n,p,2)
    
    """
    #individualExperiments(20,0.8,2)
    n=8
    p=0.3
    
    I = intervalo(1,n+1,1)
    cota1= [0]*(len(I))
    cota2=[0]*(len(I))
    cota3=[0]*(len(I))
    
    density=[0]*(len(I))
       
    for k in I:
        density[k-1] = probabilityRigidExpansion(n,p,k)
        C = cotaExpansion(n,p,k)
        cota1[k-1] = C[0]
        cota2[k-1] = C[1]
        cota3[k-1] = C[2]

    
    plt.plot(I, density,linewidth=2,color='#660066', label="Experimentation")
    #plt.plot(I, cota1,linewidth=2,color='#339933',label="Cota sup")
    plt.plot(I, cota2,linewidth=2,color='blue',label="Cota inf") 
    plt.plot(I, cota3,linewidth=2,color='red',label="Verga con venas")     
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)        
    plt.show()
    
    """    