# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 10:48:49 2017

@author: Ricardo Chávez Cáliz
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt
from numpy.linalg import qr
from numpy import log
from time import clock
from numpy import sort

limite=100000
t=1e-8

def normaSup(A):
    """
    Halla la norma en supremo de A (matrices)
    Input: Matriz A
    Output:  float
    """
    m = len (A)
    n = len (A.T)
    
    dist = 0
    for i in xrange(0,m):
        for j in xrange(0,n):
            if dist < abs(A[i,j]):
                dist= abs(A[i,j])
    return dist

def diagonal(A):
    """
    Obtiene los elementos en la diagonal de una matriz
    Input: Matriz a obtener elementos en diagonal.
    Output: Array con elementos de la diagonal.
    """
    n=len(A)
    D=np.zeros(n)
    for i in xrange (0,n):
        D[i]=A[i,i]
    return D

def mdiagonal(A):
    """
    Obtiene una matriz diagonal cuya diagonal es como la de A.
    Input: Matriz a obtener elementos en diagonal.
    Output: Matriz con elementos de la diagonal.
    """
    n=len(A)
    D=np.zeros((n,n))
    for i in xrange (0,n):
        D[i,i]=A[i,i]
    return D

def iteracionQR(A):
    """
    Obtiene eigenvalores y eigenvectores de una matriz,
    Input: array
    Output: Eigenvalores(array), eigenvectores(array n*n), num iteracion hasta
            rebasar tolerancia (int), norma(A-diag(A)) (float)
    """
    EV = np.identity(len(A))
    k=0
    
    while normaSup(A-mdiagonal(A)) > t:
        S=qr(A)
        Q=S[0]
        R=S[1]
        A=np.dot(R,Q)
        EV = np.dot(EV, Q)
        k=k+1
        if k == limite:
            print("El número de iteraciones ha exedido el límite establecido")
            print("norma de A-diag(A) es ") + str(normaSup(A-mdiagonal(A)))
            break
    return diagonal(A), EV, k , normaSup(A-mdiagonal(A)) 
    
def signo(d):
    """
    Obtiene el signo de un número (+1 0 -1), si es 0 entonces le asigna 1
    Input: float, int, etc
    Output: -1.0 o 1.0 
    """
    if d<0:
        return -1.0
    else:
        return 1.0
    
def wilkinson(A):
    """
    Calcula el desplazamiento de Wilkinson de una matriz A. 
    Input: array
    Output: float 
    """
    n = len(A)
    d = (A[n-2,n-2] - A[n-1,n-1])/2.0
    
    return A[n-1,n-1] - (signo(d)*A[n-1, n-2]**2)/(abs(d) + sqrt(d**2 + A[n-1, n-2]**2))


def iteracionDespQR(A, wilki=True):
    """
    Obtiene eigenvalores y eigenvectores de una matriz, decide el shift a realizar
    Input: arange, boolean
    Output: Eigenvalores, eigenvectores, num iteracion hasta rebasar tolerancia
            norma(A-diag(A))
    """
    n = len(A)
    EV = np.identity(n)
    k=0
    if (wilki): 
        s=wilkinson(A)
    else:
        s= A[n-1, n-1]
    
    while normaSup(A-mdiagonal(A)) > t:
        S=qr(A-s*np.identity(n))
        Q=S[0]
        R=S[1]
        A=np.dot(R,Q) + s*np.identity(n)
        EV = np.dot(EV, Q)
        k=k+1
        if k == limite:
            print("El número de iteraciones ha exedido el límite establecido")
            print("norma de A-diag(A) es ") + str(normaSup(A-mdiagonal(A)))
            break
    return diagonal(A), EV, k, normaSup(A-mdiagonal(A))

if __name__ == "__main__":
    
    """
    2. Implementa la iteración QR con shift. Aplícala a la matriz A del Ejer-
    cicio 1 con  e = 10**N para N = 1, .., 5.
    """
    #Comparativa 3 métodos.
    x = np.arange(1, 6, 1);
    R = len(x)
    y1= np.zeros(R)
    y2= np.zeros(R)
    y3= np.zeros(R)
    
    N=[1,2,3,4,5]
    for n in N:
        e = 10**n
        A = np.array([[8,1,0],
                      [1,4,e],
                      [0,e,1]])
        
        T = iteracionQR(A)
        y1[n-1] = log(T[2])
        T = iteracionDespQR(A, wilki=False)
        y2[n-1] = log(T[2])
        T = iteracionDespQR(A, wilki=True)
        y3[n-1] = log(T[2])

    plt.xlabel('valor de n')
    plt.ylabel('numero de iteraciones')
    plt.plot(x, y1, 'go', label="Sin desplazamiento")
    plt.plot(x, y2, 'r^', label="Desplazamiento simple")
    plt.plot(x, y3, 'bs', label="Desplazamiento Wilkinson")

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("eigen")
    plt.show()

    #Wilkinson solo.
    y1= np.zeros(R)
    
    for n in N:
        e = 10**n
        A = np.array([[8,1,0],
                      [1,4,e],
                      [0,e,1]])
        T = iteracionDespQR(A, wilki=True)
        y1[n-1] = T[2]

    plt.xlabel('valor de n')
    plt.ylabel('numero de iteraciones')
    plt.plot(x, y1, 'bs', label="Desplazamiento Wilkinson")

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("eigen2")
    plt.show()
        
    x = np.arange(1, 6, 1);
    R = len(x)
    
    #Contra numpy.
    y1= np.zeros(R)
    y2= np.zeros(R)
    
    N=[1,2,3,4,5]
    for n in N:
        e = 10**n
        A = np.array([[8,1,0],
                      [1,4,e],
                      [0,e,1]])
        
        tiempo_inicial = clock()
        E1 = iteracionDespQR(A, wilki=True)[0]
        y1[n-1] = (clock() - tiempo_inicial)
        
        tiempo_inicial = clock()
        E2 = np.linalg.eig(A)[0]
        y2[n-1] = (clock() - tiempo_inicial)
        
        print np.allclose(sort(E1), sort(E2))


    plt.xlabel('valor de n')
    plt.ylabel('tiempo de ejecucion')
    plt.plot(x, y1, 'bs', label="Desplazamiento Wilkinson")
    plt.plot(x, y1, '-', label="Algoritmo de numpy")


    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("WilkinsonVsNumpy")
    plt.show()

    """
    3. Determina todos los eigenvalores y eigenvectores de una matriz de
    Householder.
    """
    print "Muestra de eigenvalores para matrices de Housholder construidas a partir de vectores aleatorios" 
    n = 5
    for k in xrange (0,10):     
        v = np.random.rand(n).reshape((n,1))
        H = np.identity(n) - 2* np.dot(v,v.T)/np.dot(v.T,v)
        print np.sort(iteracionDespQR(H)[0])
    
    
    """
    5. ¿Qué pasa si aplicas la iteración QR sin shift a una matriz ortogonal?
    """
    a = 5.0
    O = np.array([[ a,a],
                  [-a,a]])
    print "Matriz ortogonal de ejemplo"
    print O
    iteracionQR(O)