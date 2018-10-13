# -*- coding: utf-8 -*-
"""
Created on Thu Sep 07 20:14:02 2017

@author: Equipo JJR
"""
import numpy as np
from numpy import exp
from numpy import log
import matplotlib.pyplot as plt

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


def FGME (muestra, I):
    """
    Provee la función generadora de momentos empirica (FGME) a partir de una muestra
    Input:  Un array de tamaño n (muestra) y un array de tamaño N que
            representa el dominio de la función (intervalo discretizado).
    Output: Un arreglo con los valores de la FGME asociados al intervalo 
            discretizado.            
    """
    
    n = len(muestra)
    N = len(I)
    M = np.zeros(N)
    muestra = muestra.astype(float)
        
    for t in xrange(0, N):
        suma = 0.0
        for i in xrange(0,n):
            suma = suma + exp (I[t]*muestra[i])
        M[t] = (1.0/n)*suma
    
    return M    

def FGPE (muestra, I):
    """
    Provee la función generadora de probabilidades empírica (FGPE) a partir 
    de una muestra
    Input:  Un array de tamaño n (muestra) y un array de tamaño N que
            representa el dominio de la función (intervalo discretizado).
    Output: Un arreglo con los valores de la FGPE asociados al intervalo 
            discretizado.
    """
    
    n = len(muestra)
    N = len(I)
    P = np.zeros(N)
    muestra = muestra.astype(float)
        
    for t in xrange(0, N):
        suma = 0.0
        for i in xrange(0,n):
            suma = suma + I[t]**(muestra[i])
        P[t] = (1.0/n)*suma
    
    return P


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
        
def testMV(M):
    return np.mean(M)/np.var(M)     

def grafica(y,b,bn,p,I):
    fig0 = plt.figure(figsize = (16,16))
    fig0.subplots_adjust(hspace=0.5)
    f, axarr = plt.subplots(2, 2)
    axarr[0, 0].plot(I, y, 'g')
    axarr[0, 1].plot(I, b, 'b')
    axarr[1, 0].plot(I, bn, 'y')
    axarr[1, 1].plot(I, p, 'm')
    axarr[0, 0].set_title('Muestra')
    axarr[0, 1].set_title('Binomial')
    axarr[1, 0].set_title('Binomial negativa')
    axarr[1, 1].set_title('Poisson')
    plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
    plt.show()

def compara(fatales, totales, artif=True):
    tam = len(fatales)
    if(artif):
        controlB = np.random.binomial(10, 0.5 , tam)
        controlBN = np.random.negative_binomial(10 , 0.5 , tam)
        controlP = np.random.poisson(1, tam)

    else: 
        T = np.mean(totales)
        F = np.mean(fatales)
        controlB = np.random.binomial(T, F/T , tam)
        controlBN = np.random.negative_binomial(T-F , 1-F/T , tam)
        controlP = np.random.poisson(F, tam)
    
    """
    PRUEBA DE FGC
    """
    I = Intervalo (-1,1,0.1)
    y = log(FGPE(fatales, I))
    b = log(FGPE(controlB, I))
    bn = log(FGPE(controlBN, I))
    p = log(FGPE(controlP, I))

    grafica(y,b,bn,p,I)
    """"
    Histogramas
    """
    fig0 = plt.figure(figsize = (16,16))
    fig0.subplots_adjust(hspace=0.5)
    f, axarr = plt.subplots(2, 2)
    axarr[0, 0].hist(fatales, color ='green')
    axarr[0, 1].hist(controlB, color = 'blue')
    axarr[1, 0].hist(controlBN, color = 'yellow')
    axarr[1, 1].hist(controlP, color= 'magenta')
    axarr[0, 0].set_title('Muestra')
    axarr[0, 1].set_title('Binomial')
    axarr[1, 0].set_title('Binomial negativa')
    axarr[1, 1].set_title('Poisson')
    plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)

    plt.show()

    """
    Media-Varianza
    """
    print "TEST COCIENTE MEDIA-VARIANZA"
    print "mayor que 1 - Binomial"
    print "igual a 1 - Poisson"
    print "menor que 1 - Binomial Negativa"

    t=[0]
    plt.plot(t, 1, 'rs', label="referencia")
    plt.plot(t, testMV(fatales), 'g^', label="muestra")
    plt.plot(t, testMV(controlB), 'bo', label="binomial")
    plt.plot(t, testMV(controlBN), 'yo', label="binomial-negativa")
    plt.plot(t, testMV(controlP), 'mo', label="poisson")
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()



def graficas(muestra):
    m = np.mean(muestra)
    v = np.var(muestra)
    
    p= (m/v)
    r= m**2/(v-m)
    print p
    print r
    
    I1 = Intervalo (0,max(muestra)+10,0.2)
    I4 = Intervalo (1.0/(1-p) - 1,1.0/(1-p),0.01)

    fde = FDE(muestra, I1)
    loge = log(FGPE(muestra, I4))
        
    """"
    Empírcas
    """
    f, (ax1,ax2) = plt.subplots(1, 2, figsize = (7,3.5))
    f.subplots_adjust(hspace=0.25)

    ax1.set_title('Distribucion')
    ax2.set_title('Log-probabilidad')    
    ax1.plot(I1,fde,linewidth=2,color='#cc0000')
    ax2.plot(I4,loge,linewidth=2,color='#cc0000')
    
    plt.show()

    """
    Media-Varianza
    """
    print "TEST COCIENTE MEDIA-VARIANZA"
    print "mayor que 1 - Binomial"
    print "igual a 1 - Poisson"
    print "menor que 1 - Binomial Negativa"

    t=[0]
    plt.figure(num=None, figsize=(1, 3.5), dpi=80)
    plt.plot([-0.5,0.5], [1,1], 'r', label="referencia")    
    plt.plot(t, testMV(muestra), 'g^', label="muestra")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()


        
def analiza(muestra):
    m = np.mean(muestra)
    v = np.var(muestra)
    
    p= (m/v)
    r= m**2/(v-m)
    print p
    print r
    
    I1 = Intervalo (0,max(muestra)+10,0.2)
    #I2 = Intervalo (-log(p)-0.2,-log(p),0.01)
    #I3 = Intervalo (-1.0/p,1/p,0.02)
    I4 = Intervalo (1.0/(1-p) - 1,1.0/(1-p),0.01)

    fde = FDE(muestra, I1)
    #fgm = FGME(muestra, I2)
    #fgp = FGPE(muestra, I3)
    loge = log(FGPE(muestra, I4))
        
    """"
    Empírcas
    """
    f, (ax1,ax2) = plt.subplots(1, 2, figsize = (7,3.5))
    f.subplots_adjust(hspace=0.25)

    ax1.set_title('Distribucion')
    ax2.set_title('Log-probabilidad')
    
    rep = 30
    for i in xrange(0,rep):
        S = np.random.negative_binomial(r , p , len(muestra))
        e = FDE(S, I1)
        ax1.plot(I1, e, color='#003366',linewidth=1,alpha=0.3)
        l = log(FGPE(S, I4))
        ax2.plot(I4, l, color='#003366',linewidth=1,alpha=0.3)

    ax1.plot(I1,fde,linewidth=2,color='#cc0000')
    ax2.plot(I4,loge,linewidth=2,color='#cc0000')
    
    plt.show()

    """
    Media-Varianza
    """
    print "TEST COCIENTE MEDIA-VARIANZA"
    print "mayor que 1 - Binomial"
    print "igual a 1 - Poisson"
    print "menor que 1 - Binomial Negativa"

    t=[0]
    plt.figure(num=None, figsize=(1, 3.5), dpi=80)
    plt.plot([-0.5,0.5], [1,1], 'r', label="referencia")
    
    rep = 30
    for i in xrange(0,rep):
        S = np.random.negative_binomial(r , p , len(muestra))
        plt.plot(t, testMV(S), 'o',markersize=5,alpha=0.3,color='#ffcc00')

    S = np.random.negative_binomial(r , p , len(muestra))
    plt.plot(t, testMV(S), 'o',markersize=5,color='#ffcc00',label="control")
    
    plt.plot(t, testMV(muestra), 'g^', label="muestra")

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()








if __name__ == "__main__":
    grupo1 = np.array([19,	8,	10,	22,	12,	20,	26,	16,	19,	23,	22, 15,	22,	16,	17,	15,	19,	20,	17,	16,	8.,	21, 12,	16,	14,	14,	9.,	9.,	19,	20,	21,	20,	17, 8.,	20,	14,	15,	14,	15,	12,	19,	19,	19,	21, 16,	16,	5.,	11,	13,	15,	9.,	18,	23,	18,	20, 20,	11,	11,	11,	13,	16,	11,	10,	25,	21,	16, 13,	12,	10,	10,	8.,	17,	15,	11,	9.,	13,	12, 14,	31,	8.,	10,	8.,	11,	14,	13,	12,	16,	10, 9.,	12,	8.,	15,	16,	15,	16,	18,	10,	16,	10, 10,	11,	17,	9.,	10,	9.,	18,	14,	12,	22,	20, 13,	20,	11,	12,	15,	14,	18,	28,	13,	30,	19, 9.,10,	18,	8.,	24,	34,	22,	17,	30,	22,	26])
    
    n, bins, patches = plt.hist(grupo1, 20, normed=1, facecolor='#006666',edgecolor='#003366', linewidth=1)
    plt.title('Histograma de datos del grupo 1')
    plt.savefig("Histograma1")
    plt.show()

    graficas(grupo1)
    analiza(grupo1)

    grupo2 = np.array([40,	54,	26,	49,	30,	26,	34,	20,	27, 19, 31,	20,	46,	23,	25,	39,	25,	26, 45,	28,	22,	67,	26,	31,	24,	23,	42, 23,	21,	41,	67,	32,	33,	33,	16,	30, 35,	28,	35,	42,	26,	34,	31,	20,	24, 24,	23,	28,	59,	58,	28,	22,	19,	21, 23,	26,	17,	33,	44,	25,	15,	17,	15, 65,	22,	43,	21,	15,	36,	12,	18,	16, 41,	21,	43,	23,	24,	32,	19,	11,	20, 36,	27,	22,	32,	18,	30,	17,	17,	21, 57,	19,	19,	26,	25,	34,	26,	22,	17, 35,	64,	33,	40,	33,	38,	26,	44,	27])    
    analiza(grupo2)
    graficas(grupo2)
    
    n, bins, patches = plt.hist(grupo2, 20, normed=1, facecolor='#006666',edgecolor='#003366', linewidth=1)
    plt.title('Histograma de datos del grupo 2')
    plt.savefig("Histograma2")
    plt.show()
