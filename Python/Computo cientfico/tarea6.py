# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 12:08:19 2017

@author: Ricardo Chávez Cáliz
"""
from numpy.random import rand
from numpy.random import beta
from numpy.random import normal

from numpy import cos
from numpy import pi
from numpy import zeros
from numpy import arange
from numpy import mean

import matplotlib.pyplot as plt

def Bernoulli(p,n):
    """
    Función que simula una muestra de tamaño n con distribución Bernoulli con parámetro p
    Input: int
    Output: array [muestra, núm de éxitos en la muestra]
    """
    M = rand(n)
    if p > 1 or p<0:
        raise ValueError("El parametro no es apropiado") 
    r = 0
    for i in xrange(0, n):
        if M[i] < p:
            M[i]=1
            r=r+1
        else:
            M[i]=0
    return M,r
        
def priori(p,n,r):
    if p<0.5 and p>0.0:
        return (p**r)*((1-p)**(n-r))*cos(pi*p)
    else:
        return 0.0

def rhoBeta(p,q):
    if q>0 and q<0.5:
        return min(1, cos(pi*q)/cos(pi*p))
    elif q<0 or q>0.5:
        return 0
    else:
        return 0

def rhoUnif(p,q,n,r):
    if q>0 and q<0.5:
        return min(1, ((q**r)*((1-q)**(n-r))*(cos(pi*q)))/((p**r)*((1-p)**(n-r))*(cos(pi*p))))
    elif q<0 or q>0.5:
        return 0
    else:
        return 0

def MHMC(n,r,tam,nombre):
    """
    Simula usando el algoritmo de Metropolis Hasting para la priori definida en el problema.
    Input: (int, int, int, funcion, funcion)
    Output: Imprime función e histograma
    """
    f, axarr = plt.subplots(2, sharex=True)
    
    #Grafica priori con r1
    x = arange(0, 0.51, 0.01)
    axarr[0].plot(x, [priori(i,n,r) for i in x])
    axarr[0].set_title('Para n=' + str(n) + " y r=" + str(r) )

    x = zeros(tam)
    x0 = (rand(1)[0])/2.0
    x[0] = x0
    
    for i in xrange(1,tam):
        if nombre == "beta":     
            y= beta(r+1,n-r+1,1)[0]
            ro = rhoBeta(x[i-1],y)
        elif nombre == "uniforme": 
            y= rand(1)[0]/2
            ro = rhoUnif(x[i-1],y,n,r)
        else:
            print "No reconozco esa instrumental, usaré beta"
            y= beta(r+1,n-r+1,1)[0]
            ro = rhoBeta(x[i-1],y)
        
        if Bernoulli(ro,1)[0][0] == 1.0:
            x[i] = y
        else:
            x[i] = x[i-1]
            

    axarr[1].hist(x)
    plt.show()    
    
if __name__ == "__main__":
    """
    1. Simular n = 5 y n = 30 v.a Bernoulli Be(1/3); sea r el 
    número de exitos en cada caso.
    """
    p=1/3.0

    #Prueba de Bernoulli
    x = arange(1, 5000, 1)
    y = zeros(len(x))
    r = zeros(len(x))
    
    for i in x:
        y[i-1] = mean(Bernoulli(p,i)[0])
        r[i-1] = p
    plt.plot(x, y, 'g', x, r, 'r')
    plt.xlabel('Tamano de la muestra')
    plt.ylabel('Promedio de la muestra')
    plt.title('Promedio v.a. Bernoulli(1/3)')
    plt.show()
    
    tam=5
    M = Bernoulli(p, tam)
    r1= M[1]
    print "para n=5 obtuvimos r=" + str(r1)
        
    tam=30
    M = Bernoulli(p, tam)
    r2= M[1]
    print "para n=30 obtuvimos r=" + str(r2)

    """
    2. Implementar el algoritmo Metropolis-Hastings para simular
    de la posterior
    
    f(p|x¯) ∝ p^r(1 − p)^(n-r)cos(πp)I[0,1/2](p),

    con los dos casos de n y r de arriba. Para ello poner la 
    propuesta (p'|p) = p' ∼ Beta(r + 1, n−r + 1) y la distribución
    inicial de la cadena µ ∼ U(0,1/2)
    """
    tam = 10000
    MHMC(5,r1,tam,"beta")
    MHMC(30,r2,tam,"beta")
        
    """
    4. Implementar el algoritmo Metropolis-Hastings con la posterior de arriba
    tomando una propuesta diferente.
    """
    MHMC(5,r1,tam,"uniforme")
    MHMC(30,r2,tam,"uniforme")
