# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 16:12:30 2017

@author: Ricardo Chávez Cáliz
"""
import numpy as np
from numpy import zeros
from numpy import log
from numpy import exp
from numpy.random import rand
from numpy.random import beta
from numpy.random import normal
from numpy.random import gamma
from numpy.random import exponential
from numpy.random import weibull
from numpy.random import poisson
from numpy.random import hypergeometric
from scipy.special import  gammaln
from math import factorial

import matplotlib.pyplot as plt

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

def Reidemacher(p):
    """
    Función que simula una v.a distribución Reidemacher con parámetro p 
    (prob de subir 1)
    Input: float
    Output: int: -1 o +1 
    """
    
    M = rand(1)
    if p > 1 or p<0:
        raise ValueError("El parametro no es apropiado") 
    if M < p:
        M=1
    else:
        M=-1
    return M

def Eco2(N,p,a,b,datos):
    m =len(datos)
    s=sum(datos)
    A=0
    for d in datos:
        A=A+gammaln(N-d+1)
    return m*gammaln(N+1) + s*log(p) + N*m*log(1-p) - s*log(1-p) + 2*(a-1)*log(p) + 2*(b-1)*log(1-p) - A 


def Eco5(N,p,a,b,datos):
    m =len(datos)
    s=sum(datos)
    A=0
    for d in datos:
        A=A+gammaln(N-d+1)
    return m*gammaln(N+1) + s*log(p) + N*m*log(1-p) - s*log(1-p) + (a-1)*log(p) + (b-1)*log(1-p) - A

def propuesta(numP,p,N,datos):
    """
    Devuelve distintas propuestas de acuerdo a la opción selecccionada.
    Input: int (opción deseada)
    Output: float, float, float (pp,Np,ro)
    """    
    if numP>3 or numP<1:
        raise ValueError("No conozco esa propuesta") 
    
    suma= sum(datos)
    #mg=max(datos)
    m=len(datos)
    
    #Modificable
    Nmax=1000
    a=1
    b=m
    #lam=1

    #Condicional total de Gibs(mueve a p)
    if numP==1:
        pp=beta(a + suma,b + m*N*suma)
        ro=1
        return pp, N,ro
    
    #Priori (mueve a las dos)
    elif numP==2:        
        pp= beta(a,b)
        Np= int((Nmax+1)*rand(1)[0]) 
        c = min(0, Eco2(Np,pp,a,b,datos) - Eco2(N,p,a,b,datos))
        
        return pp, Np,exp(c)
    
    #Caminata aleatoria (mueve a N)
    else:
        Np=N+Reidemacher(0.5)
        pp=p
        c = min(0, Eco5(Np,pp,a,b,datos) - Eco5(N,p,a,b,datos))        
        return pp, Np,exp(c) 
    
    """
    #hipergeométrica (mueve a N)
    elif numP==3:
        Np= hypergeometric(N*p,N*(1-p),N)
        return p, Np,ro
    
    #Poisson (mueve a N)
    elif numP==4:
        Np=mg+poisson(lam)
        return p, Np,ro
    """

def EcologiaMCMC(datos,tam):
    """
    Simula valores de la distribución posterior f(N,P|x), usando
    un kernel híbrido que considera las propuestas:
    Con parámetros a priori dados
    
    Input: a (float), c(float), d (float), tam(int), w1(float)
    Output: array(tam x 2) (muestra)
    """
    #Modificable
    Nm=1000
    
    M = zeros((tam, 2))#p,N
    ef=0.0
    mg = max(datos)
    
    #Inical.
    M[0,0] = rand(1)#p inicial
    M[0,1] = int((Nm-mg)*rand(1)[0])+mg #N inicial    

    for i in xrange(1,tam):
        p=M[i-1][0]
        N=M[i-1][1]
        numP = int(3*rand(1)[0])+1
        
        R = propuesta(numP,p,N,datos)
        pp = R[0]
        Np = R[1]
        ro = R[2]
        if Bernoulli(ro) == 1.0:
            M[i,0] = pp
            M[i,1] = Np

        else:
            M[i,0] = M[i-1,0]
            M[i,1] = M[i-1,1]
            ef=ef+1
    
    print "Se rechazaron el "+ str(ef*100.0/tam)+ "% de las propuestas"
    return M

#-----------------------------------------------------------------------------
def funcionLiga(x,a,b,c):
    return c* exp(-(x-a)**2/(2*b**2))

def posteriori2(a,b,c,X,Y):
    return log(c)*sum(Y) + (1/(2*(b**2)))*sum(Y*(-(X-a)**2)) - (c/(2*b**2))*sum(exp(-((X-a)**2))) - ((a-35.0)**2)/(2.0*5**2) - (5.0/2.0)*b + log(b) - (950.0/3)*c + 2*log(c)

def propuesta2(numP,a,b,c):
    """
    Devuelve distintas propuestas de acuerdo a la opción selecccionada.
    Input: int (opción deseada)
    Output: float, float, float (pp,Np,ro)
    """    
    if numP>3 or numP<1:
        raise ValueError("No conozco esa propuesta") 
    
    if numP==1:
        ap = a + normal(0,1)
        return ap, b,c
    
    elif numP==2:        
        bp = b + normal(0,0.1)
        return a, bp, c
    
    else:
        cp = c + normal(0,2)
        return a, b,cp 


def MercadoMCMC(X,Y, tam):
    """
    Simula valores de la distribución posterior f(a,b,c|x), usando
    Con parámetros a priori dados
    
    Input: datos, tam(int)
    Output: array(tam x 2) (muestra)
    """
    M = zeros((tam, 3))#a,b,c
    ef=0.0
    
    #Iniciales
    M[0,0] = normal(35,5)#a inicial
    M[0,1] = 4 #gamma(2,5.0/2)#b inicial
    M[0,2] = max(Y) #gamma(3,950.0/3)#c inicial

    for i in xrange(1,tam):
        a=M[i-1][0]
        b=M[i-1][1]
        c=M[i-1][2]

        numP = int(3*rand(1)[0])+1
        
        R = propuesta2(numP,a,b,c)

        ap = R[0]
        bp = R[1]
        cp = R[2]
    
        ro = exp(min(0, posteriori2(ap,bp,cp,X,Y) - posteriori2(a,b,c,X,Y)))
        
        if Bernoulli(ro) == 1.0:
            M[i,0] = ap
            M[i,1] = bp
            M[i,2] = cp

        else:
            M[i,0] = M[i-1,0]
            M[i,1] = M[i-1,1]
            M[i,2] = M[i-1,2]
            ef=ef+1
    
    print "Se rechazaron el "+ str(ef*100.0/tam)+ "% de las propuestas"
    return M


if __name__ == "__main__":
    """    
    1. (Problema en ecología) Sean X1, ...Xm variables aleatorias donde
    Xi denota el número de individuos de una especie en cierta region.
    Suponga que Xi|N, p ∼ Binomial(N, p), entonces
    f(x|N, p) = prod [N!/(xi!(N − xi)!)] p**xi (1 − p)**N−xi

    A partir del algoritmo MH, simule valores de la distribucion posterior
    usando un kernel hıbrido. 
    """
    datos= np.array([4, 9, 6, 7, 8, 2, 8, 7, 5, 5, 3, 9, 4, 5, 9, 8, 7, 5, 3, 2])
    M = EcologiaMCMC(datos,5000)

    X=M[1:,0]
    Y=M[1:,1]
    print np.mean(X)
    print np.mean(Y)
    
    n, bins, patches = plt.hist(Y, 40, normed=1, facecolor='#006666', alpha=0.75,edgecolor='#003366', linewidth=1.2)
    plt.title('Histograma para N')
    plt.show()

    n, bins, patches = plt.hist(X, 40, normed=1, facecolor='#006666', alpha=0.75,edgecolor='#003366', linewidth=1.2)
    plt.title('Histograma para p')
    plt.show()

    n=len(X)
    colors = np.arange(0, 1, 1.0/n)
    area = 30*np.ones(n)    
    plt.xlabel('valores de N')
    plt.ylabel('valores de p')
    plt.title('Muestra de la postterior de tamano ' + str(n+1))
    plt.scatter(Y, X, s=area, c=colors)
    plt.show()
    
    """
    2. (Estudio de mercado) Se tiene un producto y se realiza una encuesta
    con el fin de estudiar cuanto se consume dependiendo de la edad. Sea
    Yi el monto de compra y Xi la covariable la cual representa la edad.
    Suponga que Yi ∼ P o(λi) (distribuci´on Poisson con intensidad λi)
    λi = cgb(xi − a)
    para gb la siguiente funcion de liga
    gb(x) = exp −x22b2
    
    .
    O sea, se trata de regresion Poisson con una función liga no usual. Si
    λi = 0 entonces P(Yi = 0) = 1. a = a˜nos medio del segmento (a˜nos),
    c = gasto promedio (pesos), b = “amplitud” del segmento (a˜nos).
    Considere las distribuciones a priori
    a ∼ N(35, 5), c ∼ Gama(3, 3/950), b ∼ Gama(2, 2/5).
    El segundo par´ametro de la normal es desviaci´on estandard y el segundo
    par´ametro de las gammas es taza (rate).
    Usando MH simule de la distribuci´on posterior de a, c y b.
    Los datos son estos, n = 100:
    """
    X = np.array([ 17, 14, 28, 51, 16, 59, 16, 54, 52, 16, 31, 31, 54, 26, 19, 13, 59, 48, 54, 23, 50, 59, 55, 37, 61, 53, 56,
                31, 34, 15, 41, 14, 13, 13, 32, 46, 17, 52, 54, 25, 61, 15, 53, 39, 33, 52, 65, 35, 65, 26, 54, 16, 47, 14, 42, 47, 48,
                25, 15, 46, 31, 50, 42, 23, 17, 47, 32, 65, 45, 28, 12, 22, 30, 36, 33, 16, 39, 50, 13, 23, 50, 34, 19, 46, 43, 56, 52,
                42, 48, 55, 37, 21, 45, 64, 53, 16, 62, 16, 25, 62])
    
    Y = np.array([ 165, 9, 493, 0, 72, 0, 89, 0, 0, 70, 79, 96, 0, 1127, 548, 4, 0, 0, 0, 1522, 0, 0, 0, 0, 0, 0, 0, 80, 5, 38,
               0, 11, 8, 4, 31, 0, 174, 0, 0, 1305, 0, 39, 0, 0, 18, 0, 0, 4, 0, 1102, 0, 94, 0, 13, 0, 0, 0, 1308, 33, 0, 90, 0, 0, 1466,
               156, 0, 39, 0, 0, 496, 2, 1368, 190, 0, 12, 76, 0, 0, 5, 1497, 0, 6, 533, 0, 0, 0, 0, 0, 0, 0, 0, 1090, 0, 0, 0, 93, 0, 88,
               1275, 0])
    
    tam=1000
    
    #Qué pasa?
    delta=0.25
    s= np.arange(min(X), max(X)+delta, delta)
    a=normal(35,5)
    b=gamma(2,5.0/2)
    c=gamma(3,950.0/3)
    print("Parámetros aleatorios ~ priori")
    print a
    print b
    print c
    f=funcionLiga(s,a,b,c)  
    plt.scatter(X,Y,c="#003366")
    plt.plot(s, f)
    plt.show()

    #Ejecutado
    M = MercadoMCMC(X,Y,tam)
    A=M[0:,0]
    B=M[0:,1]
    C=M[0:,2]
    
    delta=0.25
    s= np.arange(min(X), max(X)+delta, delta)
    a=np.mean(A)
    b=np.sqrt(np.mean(B))
    c=np.mean(C)

    print "promedios"
    print a
    print b
    print c
    f=funcionLiga(s,a,b,c)  
    plt.scatter(X,Y,c="#003366")
    plt.plot(s, f)
    plt.show()