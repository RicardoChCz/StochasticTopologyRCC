# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 21:36:35 2017

@author: Ricardo Chávez Cáliz
"""
import numpy as np
import scipy
from scipy import stats
from scipy.stats import poisson
from scipy.special import gamma
from numpy import exp
from numpy import log
import matplotlib.pyplot as plt


#Parámetros
lam = 1
r=10
p=0.6
shape= 2.0
scale= 2.0
mu= 2.0
sigma= 2.0

#Configuración
tama = [50, 250, 1000, 2000]
k = len(tama)
rep = 100

def FDE (muestra, I):
    """
    Provee la función de distribución empírica (FDE) a partir de una muestra
    Input:  Un array de tamaño n (muestra) y un array de tamaño N que
            representa el dominio de la función (intervalo discretizado).
    Output: Un arreglo con los valores de la FDE asociados al intervalo 
            discretizado.            
    """
    
    n = len(muestra)
    N = len(I)
    F = np.zeros(N)
    
    for t in xrange(0, N):
        cuenta = 0.
        for i in xrange(0,n):
            if muestra[i] <= I[t]:
                cuenta = cuenta + 1.
                
        F[t] = (1.0/n)*cuenta
    
    return F  

def Intervalo (x0, x1, d):
    """
    Devuelve un intervalo discretizado que empieza en x0 y termina en x1
    con separación d entre cada punto de I.
    Input:  3 números: extremo izquierdo x0, extremo derecho x1,
            separación d.
    Output: Un arreglo con valores que empiezan en x0, terminan en x1 y
            difieren en d
    """
    I = np.arange(x0, x1, d)
    return I

def factorial(n):
    if n == 0:
        return 1
    
    else:
        return n * factorial(n-1)
    
def simula(nombre, j):
    if nombre == "Binomial":
        return np.random.binomial(r , p , tama[j])
    elif nombre == "BinomialNeg":
        return np.random.negative_binomial(r , p , tama[j])
    elif nombre == "Poisson":
        return np.random.poisson(lam, tama[j])
    elif nombre == "Gama":
        return np.random.gamma(shape, scale, tama[j])
    elif nombre == "Gumbel":
        return np.random.gumbel(mu, sigma, tama[j])
    else:
        raise ValueError("No conozco esa distribución") 
        

def graficas(nombre):
    """
    Devuelve las gráficas  
    Input:  Nombre de la distribución a graficar
    Output: Gráficas de FD, FGM, FGP, y log(FGP)
    """
    
    # "rep" empiricas para cada n
    FDEs = [0]*rep        
    muestras = [0]*k
    for i in xrange (0,k):
        muestras[i] =  simula(nombre, i)
    
    if nombre == "Binomial":
        #Intervalo FD de Binomial
        I1 = Intervalo (0,r+1,1)
        N1= len(I1)
        #FD Teórica 
        FD = np.zeros(N1)
        for t in xrange(0,N1):
            FD[t] =  stats.binom.cdf(I1[t], r, p)
                
    elif nombre == "BinomialNeg":
        #Intervalo FD de Binomial negativa
        I1 = Intervalo (0,20,0.5)
        N1= len(I1)
        #FD Teórica 
        FD = np.zeros(N1)
        for t in xrange(0,N1):
            FD[t] = stats.nbinom.cdf(I1[t], r,p)
                
    elif nombre == "Poisson":
        #Intervalo FD de Poisson 
        I1 = Intervalo (0,4,0.2)
        N1= len(I1)
        #FD Teórica 
        FD = np.zeros(N1)
        for t in xrange(0,N1):
            FD[t] = poisson.cdf(I1[t], lam)
                
    elif nombre == "Gama":
        #Intervalo FD de Gama
        I1 = Intervalo (0,11,0.3)
        N1= len(I1)
        #FD Teórica 
        FD = np.zeros(N1)
        for t in xrange(0,N1):
            FD[t] = scipy.stats.gamma.cdf(I1[t], shape, 0, scale)
                
    elif nombre == "Gumbel":
        #Intervalo FD de Gumbel
        I1 = Intervalo (0,8,0.3)
        N1= len(I1)
        #FD Teórica 
        FD = np.zeros(N1)
        for t in xrange(0,N1):
            FD[t] = exp(-exp(-(I1[t]-mu)/sigma))
                
    else:
        raise ValueError("No conozco esa distribución")
        
    """
    # FUNCIONES DE DISTRIBUCIÓN
    """
    #Distancia Kolmogorov-Smirnov
    
    DKS = np.zeros(rep)
    print "Distancias DKS para " + nombre
    f, axarr = plt.subplots(2, 2,figsize=(8,6),dpi=80)

    #Para cada tamaño
    for j in xrange(0,len(tama)):
        t=j%2
        s=j/2
        axarr[s, t].set_title(str(rep) + ' FDE para n=' + str(tama[j]) )
        #Para cada reptición
        for i in xrange(0,rep):
            S= simula(nombre, j)
            FDEs[i] = FDE(S, I1)
            #Distancia KS
            DKS[i] = max(FDE(S,I1)-FD)
            #Empieza a graficar
            axarr[s, t].plot(I1, FDEs[i], color='#006622',alpha=0.3)
        #Grafica la teoríca para cada tamaño            
        axarr[s, t].plot(I1, FD, color='#ff6600',linewidth=2)
        print "para n="+ str(tama[j])+ " el promedio de D_n fue " + str(np.mean(DKS))
        print "para n="+ str(tama[j])+ " el promedio de sqrt(n)*D_n fue " + str( (tama[j]**(1.0/2.0))*np.mean(DKS))
        
    plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
    plt.savefig("FD"+ nombre)
    plt.show()
            
if __name__ == "__main__":
    graficas("Binomial")
    graficas("BinomialNeg")
    graficas("Poisson")
    graficas("Gama")
    graficas("Gumbel")
    