# -*- coding: utf-8 -*-
"""
Created on Thu Oct 05 16:28:03 2017

@author: Ricardo Chávez Cáliz
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import log
from numpy import exp
from time import clock
from scipy.special import gamma
from scipy import stats

def digitosMedios(n):
    M = len(str(n))/2
    return float(str(n)[M/2:M/2 + M])

def uniformeMCL(n):
    """
    Simula una muestra  de tamaño n de variables aleatorias con distribución 
    uniforme con el método de congruencia lineal 
    Input: int(n) = tamaño de la muestra
    Output: array(tam n) = muestra
    """
    U=np.zeros(n)
    
    #Modulo temporal inicializado con números aleatorios por método del cuadrado medio
    T=[0]*5
    T[0] = digitosMedios(clock()**2)
    for i in xrange (1,5):
        T[i] = digitosMedios(T[i-1]**2)
    
    a = 10774182.0
    b = 104420.0
    m = 2**31 - 1
    
    for i in xrange(0,n):
        s = (a*T[4] + b*T[0]) % m
        U[i] = s/m
        #recorrer el estado
        for j in xrange(0,4):
            T[j]=T[j+1]
        T[4] = s
        
    return U
    
    
def exponencial(lam,n):
    """
    Simula una v.a. exponencial con parametro lam y tamaño n.
    Input: float, int
    Output: array
    """
    y=np.random.rand(n)
    return (-1.0/float(lam))*log(1.0-y)

#------------------------------------------------------------------------------

def recta(x,a,b):
    """
    Devuelve el valor de la recta formada por los puntos (a,logf(a)) y (b,logf(b)).
    Input: tritupla (3-floats)
    Output: float
    """
    fa=log(f(a))
    return ((fa-log(f(b)))/(a-b))*(x-a) + fa

def alphasBetas(x):
    """
    Calcula los alpha y betas asociados a g_n(x) para cada intervalo [xi,x(i+1)]
    Input: array (tam n+2) (puntos en el soporte)
    Output (array tam n+1,array tam n+1) 
    """
    r = len(x) - 1 #r=n+1
    a = np.zeros(r)
    b = np.zeros(r)
    
    for i in xrange(r):
        a[i] = (log(f(x[i])) - log(f(x[i+1])))/(x[i]-x[i+1])
        b[i] = log(f(x[i])) - a[i]*x[i]

    return a,b

def f(y):
    """
    Densidad de gama con parámetro de forma  a y de razón b
    """
    a=2
    b=1
    return (b**a) * (y**(a-1)) * exp(-b*y) / (gamma(a))
    
def f_(z,x):
    """
    Cota inferior y superior para densidad
    """
    k=len(x)-1 #k=n+1

    if z<x[0] or z>x[k]:
        return 0.0,exp(min(recta(z,x[0],x[1]),recta(z,x[k-1],x[k])))

    elif z<x[1]:
        return exp(recta(z,x[0],x[1])),exp(recta(z,x[1],x[2]))
    elif z<x[k] and z>x[k-1]:
        return exp(recta(z,x[k-1],x[k])), exp(recta(z,x[k-2],x[k-1]))
    else:
        i=1
        while z>x[i]:
            i=i+1
        i=i-1
        return exp(recta(z,x[i],x[i+1])), exp(min((recta(z,x[i-1],x[i])),(recta(z,x[i+1],x[i+2]))))
    
def omega(x,a,b):
    """
    Calcula la constante de normalización para la g_k
    Input: 3 array (n) 
    Output: float
    """    
    suma = 0.0
    for i in range(len(a)):
        suma = suma + exp(b[i]) * (exp(a[i]*x[i+1]) - exp(a[i]*x[i]))/a[i]
    
    return suma
    
def probabilidades(x,w,a,b):
    """
    Deueleve un vector con la probabilidades acumuladas para escojer el intervalo
    [xi,x(i+1)] con probabilidad pi, es decir array [p0,p0+p1, p0+p1+p2,....,1]
    """
    k=len(x)-1 #n+1
    p=np.zeros(k)
    for i in xrange(0,k):
        p[i] = exp(b[i]) * (exp(a[i]*x[i+1]) - exp(a[i]*x[i])) / (w*a[i])
    return np.cumsum(p)

def escojeI(p):
    """
    Recibe un vector de probabilidades acumuladas y devuelve el valor i del intervalo
    a seleccionar [xi, x(i+1)]
    """
    o = np.random.rand(1)[0]
    i=0
    while o>p[i]:
        i=i+1        
    return i

def simulaG(x,w,p,a,b):
    """
    Función auxiliar para paso 1 del algoritmo
    """
    i=escojeI(p)
    u = np.random.rand(1)[0]
    return (a[i])**(-1) * log(exp(a[i]*x[i]) + u*(exp(a[i]*x[i+1]) - exp(a[i]*x[i])))

def ARS(N):
    """
    Algoritmo para un Adpative Rejection Sampler
    Input: int (tamaño de la muestra)
    Output: array (tamaño N)
    """
    M= np.zeros(N)
        
    #Inicializar varaibles auxiliares
    ef = 0 
    i = 0
    
    #1 - Inizializar n y x
    
    n = 3
    x=uniformeMCL(n)*15
    x = np.array(sorted(np.concatenate((x,np.array([0.01,15.0])),axis=0)))
    AlfaBetas= alphasBetas(x)
    a=AlfaBetas[0]
    b=AlfaBetas[1]
    w=omega(x,a,b)
    p=probabilidades(x,w,a,b)
    
    while M[N-1]==0:
        #2 - Generar X y U
        z = simulaG(x,w,p,a,b)
        u = np.random.rand(1)[0]
    
        #3 - Decide si acepta, adapta, rechaza.
        cotas = f_(z,x)
        fI=cotas[0]
        fS=cotas[1]
        if u<= (fI/fS):
            M[i] = z
            i=i+1
        elif u<= (f(z)/fS):
            M[i] = z
            x = np.array(sorted(np.concatenate((x,np.array([z])),axis=0)))
            AlfaBetas= alphasBetas(x)
            a=AlfaBetas[0]
            b=AlfaBetas[1]
            w=omega(x,a,b)
            p=probabilidades(x,w,a,b)
            i=i+1
        else:
            ef=ef+1
        
    return M,x

if __name__ == "__main__":
    
    """
    2.Implementar el siguiente algoritmo para simular variables aleatorias uniformes:
    $$x_i = 107374182x_{i-1} + 104420x_{i-5} \text{ (mod } 2^{31} - 1)$$
    regresa $x_i$ y recorrer el estado, esto es $x_{j-1} = x_j; j = 1, 2, 3, 4, 5;$ ¿parecen U(0, 1)?
    """
    n = 10000
    tiempo_inicial = clock()
    muestra=uniformeMCL(n)
    print "Tiempo de ejecución con MCL ", (clock() - tiempo_inicial) , " para N=", n
    print stats.kstest(muestra, 'uniform')
    
    tiempo_inicial = clock()
    control=np.random.rand(n)
    print "Tiempo de ejecución con numpy ", (clock() - tiempo_inicial) , " para N=", n
    print stats.kstest(control, 'uniform')
    I = np.arange(0, 1.2, 0.2)

    fig, axarr = plt.subplots(1, 2, figsize=(7,3.5),dpi=80)
    axarr[0].hist(muestra, bins=25, alpha=0.5, facecolor='#cc0000',edgecolor='#800000', linewidth=1, normed=1)
    axarr[1].hist(control, bins=25, alpha=0.5, facecolor='#006699',edgecolor='#003366', linewidth=1, normed=1)
    axarr[0].set_title('MCL')
    axarr[0].set_ylabel('Frecuencia relativa')
    axarr[0].set_xlabel('Valor')    
    axarr[1].set_title('Numpy')
    axarr[1].set_xlabel('Valor')
    axarr[0].plot(I, np.ones(len(I)), color="#000000",linewidth=1.5, label="densidad")
    axarr[1].plot(I, np.ones(len(I)), color="#000000",linewidth=1.5, label="densidad")
    plt.savefig("HistoUni")
    plt.show()
    
    n = 1000
    muestra=uniformeMCL(n)
    m2=np.zeros(n)
    for i in range(n-1):
        m2[i]=muestra[i+1]
    m2[n-1]=muestra[0]    
    #Scatter
    plt.figure(num=None, figsize=(3, 3), dpi=80)
    colors = np.arange(0, 1, 1.0/n)
    area = 20*np.ones(n)

    control=np.random.rand(n)
    c2=np.zeros(n)
    for i in range(n-1):
        c2[i]=control[i+1]
    c2[n-1]=control[0]    

    #Scatter
    colors = np.arange(0, 1, 1.0/n)
    area = 20*np.ones(n)
    fig, axarr = plt.subplots(1, 2, figsize=(7,3.5),dpi=80)
    axarr[0].scatter(muestra, m2, s=area, c=colors)
    axarr[1].scatter(control, c2, s=area, c=colors) 
    axarr[0].set_title('MCL') 
    axarr[1].set_title('Numpy')
    plt.savefig("pruebaUni.png")
    plt.show()

    """
    5.Implementar el algoritmo Adaptive Rejection Sampling y simular de una Gama(2, 1) 10,000 muestras. ¿cuándo es conveniente dejar de adaptar la envolvente?
    """
    N = 10000
    
    tiempo_inicial = clock()
    expo=exponencial(1,N) + exponencial(1,N)
    print "Tiempo de ejecución con exponenciales ", (clock() - tiempo_inicial) , " para N=", N
    print stats.kstest(expo, 'gamma', args=(2, 0, 1))
    
    tiempo_inicial = clock()
    cont=np.random.gamma(2,1,N)
    print "Tiempo de ejecución con Numpy ", (clock() - tiempo_inicial) , " para N=", N
    print stats.kstest(cont, 'gamma', args=(2, 0, 1))
    
    tiempo_inicial = clock()
    resultado=ARS(N)
    print "Tiempo de ejecución con ARS ", (clock() - tiempo_inicial) , " para N=", N
    arsM=resultado[0]
    print stats.kstest(arsM, 'gamma', args=(2, 0, 1))
    
    I = np.arange(0, max(max(arsM),max(cont),max(expo)), 0.2)    

    fig, axarr = plt.subplots(1, 3, figsize=(9,3),dpi=80)
    axarr[0].hist(expo, bins=25, alpha=0.5, facecolor='#cc0000',edgecolor='#800000', linewidth=0.7, normed=1)
    axarr[1].hist(cont, bins=25, alpha=0.5, facecolor='#006699',edgecolor='#003366', linewidth=0.7, normed=1)
    axarr[2].hist(arsM, bins=25, alpha=0.5, facecolor='#339933',edgecolor='#006600', linewidth=0.7, normed=1) 
    axarr[0].plot(I, f(I), color="#000000",linewidth=1.5, label="densidad")
    axarr[1].plot(I, f(I), color="#000000",linewidth=1.5, label="densidad")
    axarr[2].plot(I, f(I), color="#000000",linewidth=1.5, label="densidad")
    axarr[0].set_title('Exponenciales')    
    axarr[1].set_title('Numpy')
    axarr[2].set_title('ARS')
    axarr[0].set_ylabel('Frecuencia relativa')
    axarr[0].set_xlabel('Valor')
    axarr[1].set_xlabel('Valor')
    axarr[2].set_xlabel('Valor')
    plt.savefig("Gammas.png")
    plt.show()
    
    print "Haciendo pruebas de comparación de tiempos entre Numpy y ARS..."
    
    R=10500
    #Contra numpy.
    
    T = np.arange(1, R, 500)    
    y1= np.zeros(len(T))
    y2= np.zeros(len(T))

    for t in xrange(0,len(T)):
        tiempo_inicial = clock()
        cont=np.random.gamma(2.5,1.5,T[t])
        y1[t] = log(clock() - tiempo_inicial)
        
        tiempo_inicial = clock()
        ARS(T[t])
        y2[t] = log(clock() - tiempo_inicial)
        
    plt.xlabel('tamano de la muestra')
    plt.ylabel('log(tiempo de ejecucion)')
    plt.plot(T, y1, linestyle='--', color='#003366', label="Numpy")
    plt.plot(T, y2, linestyle='-', color='#006600', label="ARS")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("NumpyVsARS")
    plt.show()

