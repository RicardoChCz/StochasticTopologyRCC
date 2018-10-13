# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 19:46:48 2017

@author: Ricardo Chávez Cáliz
"""
from numpy.random import rand
from numpy.random import normal
from numpy.random import exponential
from numpy.linalg import cholesky

from numpy.random import gamma as Rgamma
from numpy import exp
from numpy import log
from numpy import zeros
from numpy import arange
from numpy import array
from numpy import mean
from numpy import dot
from numpy import transpose
from numpy import meshgrid

import numpy as np

from scipy.special import gamma
from matplotlib.mlab import bivariate_normal

import matplotlib.pyplot as plt

import matplotlib.collections as mcoll
import matplotlib.path as mpath


def densObj1(a,b,r1,r2,n):
    c = n*a*log(b) + (a-1)*log(r1) - b*(r2+1) - n*log(gamma(a))
    return c

def vectorGauss(k):
    m = np.array([[3],
                  [5]])
    sigma = np.array([[1, 0.8],
                      [0.8, 1]])
    #Aplicar Cholesky
    U = np.linalg.cholesky(sigma)
    
    #Generar vector aleatorio de dimensión k con distribución N(0,1)
    Z = np.zeros([k, 1], dtype=float)
    for i in range(0, k):
        Z[i]=np.random.normal(0,1)

    #Generar vector aleatorio de dimensión 1k con ditribución N(m, sigma)
    X = m + np.dot(np.transpose(U),Z)
    return X

def muestraNMV(n,k):    
    A=np.zeros((n,k))
    for i in range(0, n):
        X=vectorGauss(k)
        for j in range(0,k):
            A[i,j] = X[j,0]
    return A

def grafBivar(cov_test):
    """
    Función que grafica bivariada
    """
    # Devueleve eigenvalores
    Eig = np.linalg.eig(cov_test)
    
    # Construir eigenvectores normalizados
    V1 = np.array([Eig[1][0][0], Eig[1][1][0]])
    V2 = np.array([Eig[1][0][1], Eig[1][1][1]])
    
    # Modificar eigen valores N_{i} = alpha*raíz de lambda_{1} * V_{i}
    alpha = 3
    N1 = np.array([V1[0] * alpha * np.sqrt(abs(Eig[0][0])) , V1[1] * alpha * np.sqrt(abs(Eig[0][0])) ])
    N2 = np.array([V2[0] * alpha * np.sqrt(abs(Eig[0][1])) , V2[1] * alpha * np.sqrt(abs(Eig[0][1])) ])
    
    # Define los límites considerando la mayor proyección de los vectores modificados
    # N1 y N2 en x y en y respectivamente.
    lx = max(abs(N1[0]), abs(N2[0]))
    ly = max(abs(N1[1]), abs(N2[1]))
    
    #Da el dominio de graficación X Y usando los limites y trasladando a la media.
    x = np.arange(-lx + 3, lx + 3, 0.1)
    y = np.arange(-ly + 5, ly + 5, 0.1)
    X, Y = np.meshgrid(x, y)
    
    #Te da una bivariada con dominio  X Y, sigma1, sigma2, mu1, mu2, pho*s1*s2 y la grafíca.
    Z = bivariate_normal(X, Y, 1, 1, 3, 5, 0.8)
    plt.contour(X,Y,Z)
    plt.savefig('normalBivariada.png')

def colorline(x, y, z, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0),linewidth=3, alpha=0.6):
    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha)

    ax = plt.gca()
    ax.add_collection(lc)

    return lc


def make_segments(x, y):
    """
    Crea una lista de segmentos de cordenaas x y y, en el formato apropiado 
    para LineCollection: un array de la forma num lineas * (num de puntos 
    por línea) * 2 (x y y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


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
        

def normalMV(m,sigma):
    k= len(m)
    #Aplicar Cholesky
    U = cholesky(sigma)
    
    #Generar vector aleatorio de dimensión 10 con distribución N(0,1)
    Z = zeros([k, 1])
    for i in range(0, k):
        Z[i]=normal(0,1)

    #Generar vector aleatorio de dimensión 10 con ditribución N(m, sigma)
    X = m + dot(transpose(U),Z)

    return X
    

def MHMC(tam,s1,s2,r1,r2,n):
    """
    Simula usando el algoritmo de Metropolis Hasting para la priori definida en el problema.
    Input: (int, int, int, funcion, funcion)
    Output: array
    """
    dim = 2
    M = zeros((tam, dim))
    #Distribuciones a priori
    a0 = 1.0 + 3.0*rand(1)[0]
    b0 = exponential(1)

    M[0,0] = a0
    M[0,1] = b0
    
    ef=0.0    
    
    for i in xrange(1,tam):
        a=M[i-1,0]
        b=M[i-1,1]
    
        ap= a+normal(0,s1)
        bp= b+normal(0,s2)
        
        if ap<1 or ap>4 or bp<1:
            ro = 0
        else:
            w = n*ap*log(bp) + (ap-a)*log(r1) +(b-bp) + n*log(gamma(a)) -n*a*log(b) - n*log(gamma(ap))
            if w<0:
                ro=exp(w)
            else:
                ro=1
        
        if Bernoulli(ro,1)[0][0] == 1.0:
            M[i,0] = ap
            M[i,1] = bp
        else:
            M[i,0] = M[i-1,0]
            M[i,1] = M[i-1,1]
            ef=ef+1
    
    print "se rechazaron el " + str(ef/tam)+'%de los puntos propuestos'
    return M

def gamaEntera(z,tam):
    return Rgamma(z,1,tam)
    

def MHMC2(tam,a,a0):
    """
    Simula usando el algoritmo de Metropolis Hasting para la priori definida en el problema.
    Input: (int, int, int, funcion, funcion)
    Output: array
    """
    M = zeros(tam)
    M[0] = a0
    
    ef=0.0
    for i in xrange(1,tam):
        y= gamaEntera(a,1)[0]
        d=a-int(round(a))
        ro = min(1, (y**(d))/(M[i-1]**(d)) )
        
        if Bernoulli(ro,1)[0][0] == 1.0:
            M[i] = y
        else:
            M[i] = M[i-1]
            ef=ef+1
            
    print "Porcentaje de rechazos " + str(ef/tam)
    return M

def ev4(y1,y2):
    return (-1/(2*0.36))*((y1-3.0)**2 + (y2-5.)**2 - 2.*0.8*(y1-3.)*(y2-5.))

def MHMC3(tam,a0,b0,s):
    """
    Simula usando el algoritmo de Metropolis Hasting para la priori definida en el problema.
    Input: (int, int, int, funcion, funcion)
    Output: array
    """
    dim = 2
    M = zeros((tam, dim))
    M[0,0] = a0
    M[0,1] = b0
    
    for i in xrange(1,tam):        
        y= M[i-1] + normal(0,s,2)
        y0= y[0]
        y1= y[1]
        x0=M[i-1][0]
        x1=M[i-1][1]
        
        Ey=ev4(y0, y1)
        Ex=ev4(x0, x1)
        
        if (Ey-Ex<0):
            ro = exp(Ey-Ex)
        else:
            ro = 1
        
        if Bernoulli(ro,1)[0][0] == 1.0:
            M[i,0] = y0
            M[i,1] = y1

        else:
            M[i,0] = M[i-1,0]
            M[i,1] = M[i-1,1]
    
    return M

def curvasGamma(n,r1,r2):
    delta = 0.025
    x = np.arange(1, 4+delta, delta)
    y = np.arange(80, 120, delta)
    X, Y = np.meshgrid(x, y)
    Z = densObj1(X,Y,r1,r2,n)    
    ma= np.amax(Z)
    mi= np.amin(Z)
    
    plt.figure()
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    ax1.contour(X, Y, Z, levels=arange(mi,ma,1))
    ax1.set_title('log(f) para n = ' + str(n))
    
    #Graficar original
    Z = exp(Z)
    ax2.contourf(X, Y, Z)
    ax2.set_title('f para n = ' + str(n))    
    plt.savefig('CurvasNivel'+str(n)+'.png')
    plt.show()

def erres(x):
    r1=1.0
    r2=0.0
    for i in xrange(0, len(x)):
        r1=r1*x[i]
        r2=r2+x[i]
        
    return r1,r2

def RWMH(M,N,a0,b0):
    cov_test = np.array([[1, 0.8],
                         [0.8, 1]])
    A= M[:,0]
    B= M[:,1]
    x = (A).tolist()
    y = (B).tolist()

    #Scatter    
    colors = np.arange(0, 1, 1.0/N)
    area = 50*np.ones(N)
    plt.scatter(x, y, s=area, c=colors, alpha=0.5)
    grafBivar(cov_test)
    plt.savefig('ScatterLevels'+str(N)+'-'+str(a0)+','+str(b0)+'.png')
    plt.show()
    
    x = (A).tolist()
    y = (B).tolist()

    path = mpath.Path(np.column_stack([x, y]))
    verts = path.interpolated(steps=3).vertices
    x, y = verts[:, 0], verts[:, 1]    
    colorline(x, y, colors, cmap=plt.get_cmap('jet'), linewidth=1)
    grafBivar(cov_test)
    plt.savefig('RandomWalk'+str(N)+'-'+str(a0)+','+str(b0)+'.png')

    plt.show()
    
def evalua(M):
    n=len(M)
    x = np.arange(0., n, 1)
    y = log(M)
    colors = np.arange(0, 1, 1.0/n)
    area = 30*np.ones(n)
    plt.scatter(x, y, s=area, c=colors, alpha=0.8)
    plt.show()
    

if __name__ == "__main__":
    """
    Con el algoritmo Metropolis-Hastings (MH), simular lo siguiente:
    1. Sean xi ∼ Ga(α, β); i = 1, 2, . . . , n. Simular datos xi con α = 3 y
    β = 100 considerando los casos n = 5 y 30.
    Con α ∼ U(1,4), β ∼ exp(1) distribuciones a priori, se tiene la posterior
    f(α, β|x¯) ∝ (β**nα / Γ(α)**n) * r**α−1 * e**(−β(r2+1)) * 1(1 ≤ α ≤ 4)1(β > 1),
    
    con r2 = sum x_i
    y r1 = prod x_i
    .
    En ambos casos, grafica los contornos para visualizar d´onde est´a concentrada
    la posterior.
    Utilizar la propuesta
    
    q(αp,βp|α,β) = (α,β) + (ε1,ε2)

    donde (ε1,ε2) ∼ N2 ((0,0),(σ1^2,0
                               0,σ2^2))
    """
    s1=0.2
    s2=1
    tam=10000
    N=tam
    
    n=5
    x = Rgamma(3,(1.0/100),n)
    r1=erres(x)[0]
    r2=erres(x)[1]
    
    print "r1 vale " + str(r1)
    print "r2 vale " + str(r2)
    
    curvasGamma(n,r1,r2)
    M1 = MHMC(tam,s1,s2,r1,r2,n)
    
    A= M1[:,0]
    B= M1[:,1]
    x = (A).tolist()
    y = (B).tolist()

    #Scatter    
    colors = np.arange(0, 1, 1.0/N)
    area = 50*np.ones(N)
    plt.scatter(x, y, s=area, c=colors, alpha=0.5)
    plt.savefig('ScatterGamma5.png')
    plt.show()
    
    n=30
    x = Rgamma(3,(1.0/100),n)
    r1=erres(x)[0]
    r2=erres(x)[1]
    curvasGamma(n,r1,r2)
    M1 = MHMC(tam,s1,s2,r1,r2,n)
    
    A= M1[:,0]
    B= M1[:,1]
    x = (A).tolist()
    y = (B).tolist()

    #Scatter    
    colors = np.arange(0, 1, 1.0/N)
    area = 50*np.ones(N)
    plt.scatter(x, y, s=area, c=colors, alpha=0.5)
    plt.savefig('ScatterGamma30.png')
    plt.show()
    
    evalua(M1[:,0])
    evalua(M1[:,1])
    print "estimada de alfa" + str(mean(M1[:,0]))
    print "estimada de beta" + str(mean(M1[:,1]))
    
    """
    2. Simular de la distribuci´on Gamma(α,1) con la propuesta Gamma([α],1),
    donde [α] denota la parte entera de α.
    Adem´as, realizar el siguiente experimento: poner como punto inicial
    x0 = 1, 000 y graficar la evoluci´on de la cadena, es decir,
    f(Xt) vs t.
    """
    n=1000
    
    a1= MHMC2(n,3.4,1)  
    a2= MHMC2(n,3.4,10)  
    a3= MHMC2(n,3.4,100)  
    a4= MHMC2(n,3.4,1000)  
    
    f, axarr = plt.subplots(2, 2)
    axarr[0, 0].hist(a1, bins='auto')
    axarr[0, 0].set_title('Estado inicial en 1')
    axarr[0, 1].hist(a2, bins='auto')
    axarr[0, 1].set_title('Estado inicial en 10')
    axarr[1, 0].hist(a3, bins='auto')
    axarr[1, 0].set_title('Estado inicial en 100') 
    axarr[1, 1].hist(a4, bins='auto')
    axarr[1, 1].set_title('Estado inicial en 1000') 
    plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
    plt.savefig('histo2.png')
    plt.show()

    x = np.arange(0., n, 1)
    y = log(a1)
    colors = np.arange(0, 1, 1.0/n)
    area = 30*np.ones(n)
    
    f, axarr = plt.subplots(2, 2)
    axarr[0, 0].scatter(x, log(a1), s=area, c=colors, alpha=0.8)
    axarr[0, 0].set_title('Estado inicial en 1')
    axarr[0, 1].scatter(x, log(a2), s=area, c=colors, alpha=0.8)
    axarr[0, 1].set_title('Estado inicial en 10')
    axarr[1, 0].scatter(x, log(a3), s=area, c=colors, alpha=0.8)
    axarr[1, 0].set_title('Estado inicial en 100') 
    axarr[1, 1].scatter(x, log(a4), s=area, c=colors, alpha=0.8)
    axarr[1, 1].set_title('Estado inicial en 1000') 
    plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
    plt.savefig('scatter2.png')
    plt.show()

    """
    3. Implementar Random Walk Metropolis Hasting (RWMH) donde la distribuci´on
    objetivo es N2(µ, Σ), con
    µ = (3,5) Σ = (1, 0.8
                   0.8, 1)
    
    Utilizar como propuesta εt ∼ N2 (0, σI). ¿C´omo elegir σ para que la
    cadena sea eficiente? ¿Qu´e consecuencias tiene la elecci´on de σ?
    Como experimento, elige como punto inicial xo = (1000,1)
    y comenta los resultados.
    
    Para todos los incisos del ejercicio anterior:
    • Establece cual es tu distribuci´on inicial.
    • Grafica la evoluci´on de la cadena.
    • Indica cu´al es el Burn-in.
    • Comenta qu´e tan eficiente es la cadena.
    • Implementa el algoritmo MH considerando una propuesta diferente.
    """
    
    N=1000
    s=1
    M = MHMC3(N,3,5,s)
    RWMH(M,N,3,5)
    
    M = MHMC3(N,1000,1,20)
    RWMH(M,N,1000,1)

    M = MHMC3(N,100,1,1)
    RWMH(M,N,100,1)
    
    evalua(M[:,0])
    evalua(M[:,1])
