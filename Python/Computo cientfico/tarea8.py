# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 17:55:59 2017

@author: Ricardo Chávez Cáliz
"""
import numpy as np

from numpy import zeros
from numpy import log
from numpy import exp

from numpy.random import rand
from numpy.random import normal
from numpy.random import gamma
from numpy.random import exponential
from numpy.random import weibull
from numpy.random import poisson
from scipy.special import gammaln as loggamma

from matplotlib.mlab import bivariate_normal

import matplotlib.pyplot as plt

def Bernoulli(p):
    """
    Función que simula una v.a distribución Bernoulli con parámetro p
    Input: int
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

def vectorGauss(m,sigma):
    k=len(m)
    
    #Aplicar Cholesky
    U = np.linalg.cholesky(sigma)
    
    #Generar vector aleatorio de dimensión k con distribución N(0,1)
    Z = np.zeros([k, 1], dtype=float)
    for i in range(0, k):
        Z[i]=np.random.normal(0,1)

    #Generar vector aleatorio de dimensión 1k con ditribución N(m, sigma)
    X = m + np.dot(np.transpose(U),Z)
    return X

def muestraNMV(m,sigma,n):
    k=len(m)
    A=np.zeros((n,k))
    for i in range(0, n):
        X=vectorGauss(m,sigma)
        for j in range(0,k):
            A[i,j] = X[j,0]
    return A

def grafBivariada(m,sigma):
    """
    Función que grafica contornos de nivel  de bivariada
    """
    m1=m[0,0]
    m2=m[1,0]
    s1=np.sqrt(sigma[0,0])
    s2=np.sqrt(sigma[1,1])
    ro=(sigma[0,1])/(s1*s2)
    
    # Devueleve eigenvalores
    Eig = np.linalg.eig(sigma)
    
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
    x = np.arange(-lx + m1, lx + m1, 0.1)
    y = np.arange(-ly + m2, ly + m2, 0.1)
    X, Y = np.meshgrid(x, y)
    
    #Te da una bivariada con dominio  X Y, sigma1, sigma2, mu1, mu2, pho*s1*s2 y la grafíca.
    Z = bivariate_normal(X, Y, s1, s2, m1, m2, ro*s1*s2)
    plt.contour(X,Y,Z, alpha=0.5)
    
    
def graficaMuestraBi(m,sigma,w1,n):
    """
    Grafica una muestra normal bivriada de tamaño n con parámetros m y sigma
    usando el metodo de MH con kerneles híbridos de parametro w1, junto con 
    los contornos de nivel de la densidad correspondiente.
    Input: array, array, float, int (media, matriz de covarianza, probabiidad de
           tomar el primer Kernel, tamaño)
    Output: Muestra gráfica
    """    
    M1 = NormalMHMC(m,sigma,w1,n)
    
    A= M1[:,0]
    B= M1[:,1]
    x = (A).tolist()
    y = (B).tolist()

    #Scatter
    colors = np.arange(0, 1, 1.0/n)
    area = 50*np.ones(n)
    plt.scatter(x, y, s=area, c=colors, alpha=0.8)
    grafBivariada(m,sigma)
    plt.savefig('bivariadaMH'+str(ro)+'-'+str(n)+'.png')
    plt.title('Muestra de tamano '+str(n)+ ' con rho = '+ str(ro))
    plt.show()

def densidad2(a,l,T,b,c):
    n=len(T)
    r1=1.0
    for t in T:
        r1 = r1 * t
        
    sa=0
    for t in T:
        sa = sa + t**a
        
    suma=0
    for t in T:
        suma = suma + log(t)

    return ((n+a-1)*log(l)) - (l*(b+sa)) + (log(b)*a ) + (log(a)*n) + (a-1)*suma - (c*a) - loggamma(a)
    
    
def grafTiempos(M,T,b,c):
    """
    Función que grafica bivariada
    """
    A=min(0.1, min(M[:,0]))
    B=max(M[:,0])
    C=min(0.1, min(M[:,1]))
    D=max(M[:,1])

    #Contornos
    delta = 0.25
    x = np.arange(A, B+ 2*delta, delta)
    y = np.arange(C, D+ 2*delta, delta)
    X, Y = np.meshgrid(x, y)
    Z = densidad2(X,Y,T,b,c)
    ma= np.amax(Z)
    mi= np.amin(Z)
    
    plt.figure()
    plt.contour(X, Y, Z, levels=np.arange(mi,ma,5))

def NormalMHMC(m,sigma,w1,tam):
    """
    Aplica el algoritmo de Metropolis-Hastings considerando como funcion 
    objetivo la distribución normal bivariada con Kerneles híbridos dados
    por las siguientes propuestas:      
        q1 ((x01, x02)|(x1, x2)) = f_X1|X2(x01|x2)1(x02 = x2)
        q2 ((x01, x02)|(x1, x2)) = fX2|X1(x02|x1)1(x01 = x1)
    Input: media, matriz de covarianza, tamaño (array, array, int)
    Output: muestra (array)
    """
    dim = 2
    M = zeros((tam, dim))

    m1=m[0,0]
    m2=m[1,0]
    s1=np.sqrt(sigma[0,0])
    s2=np.sqrt(sigma[1,1])
    ro=(sigma[0,1])/(s1*s2)
    
    #Inical.
    M[0,0] = m1
    M[0,1] = m2
    
    for i in xrange(1,tam):
        if Bernoulli(w1)==1:
            x2=M[i-1,1]
            M[i,0]=normal(m1+ro*s1*(x2-m2)/s2, (s1**2)*(1-ro**2))
            M[i,1]=M[i-1,1]
            
        else:
            x1=M[i-1,0]
            M[i,1]=normal(m2+ro*s2*(x1-m1)/s1, (s2**2)*(1-ro**2))
            M[i,0]=M[i-1,0]
            
    return M

def propuesta(numP,a,l,b,c,sigma,T):
    """
    Devuelve distintas propuestas de acuerdo a la opción selecccionada.
    Input: int (opción deseada)
    Output: float, float (propuesta, ro)
    """
    n=len(T)
    
    if numP>4 or numP<1:
        raise ValueError("No conozco esa propuesta") 
    
    if numP==1:
        sa=0
        for t in T:
            sa = sa + t**a
        #propuestas
        lp= gamma(a + n , 1.0/(b + sa))
        ap= a
        #rho
        ro = 1
        return ap,lp,ro

    elif numP==2:
        r1=1.0
        for t in T:
            r1 = r1 * t
        #propuestas
        ap = gamma(n + 1 , 1.0/(-log(b)-log(r1)+c)) 
        lp = l
        
        #rho
        sap=0
        for t in T:
            sap = sap + t**ap
        sa=0
        for t in T:
            sa = sa + t**a

        aux = float(loggamma(a)) + (ap-a)*l - l*sap-float(loggamma(ap)) + l*sa
        c=min(0,aux)
        return ap,lp,exp(c)

    elif numP==3:
        ap = exponential(c)
        lp = gamma(ap , 1.0/b)

        #rho
        sap=0.0
        for t in T:
            sap = sap + t**ap
        sa=0.0
        for t in T:
            sa = sa + t**a
            
        r1=1.0
        for t in T:
            r1= r1*t

        aux = n*log((ap*lp)/(a*l)) + (ap-a)*log(r1)- lp*sap + l*sa
        c=min(0,aux)
        return ap,lp,exp(c)

    else:
        ap=a + normal(0,sigma)
        lp= l
        
        suma=0
        for t in T:
            suma = suma + log(t)
        sap=0
        for t in T:
            sap = sap + t**ap
        sa=0
        for t in T:
            sa = sa + t**a

        aux =  ap*log(l)-l*(b+sap)+ap*log(b)+n*log(ap)+(ap-1)*suma-c*ap-float(loggamma(ap))-a*log(l)+l*(b+sa)-a*log(b)-n*log(a)-(a-1)*suma+c*a+float(loggamma(a))
        
        c=min(0,aux)
        return ap,lp,exp(c)


def TiemposMHMC(c,b,sigma,tam,T):
    """
    Aplica el algoritmo MH usando Kerneles híbridos para simular valores de
    la distribución posterior f(α, λ|t ̄) ∝ f(t ̄|α, λ)f(α, λ), considerando 
    las siguientes propuestas:
    Propuesta 1:
        λp|α,t ̄∼ Gamaα + n , b +Xni=1tαi!
    Propuesta 2:
        αp|λ,t ̄ ∼ Gama (n + 1 , −log(b) − log(r1) + c)
    Propuesta 3:
        αp ∼ exp(c) y λp|αp ∼ Gama(αp, b)
    Propuesta 4 (RWMH):
        αp = α + , con  ∼ N(0, σ) y dejando λ fijo.
    
    con distribuciones a priori datos simulado usando α = 1 y λ = 1 con n = 20
    c = 1 y b = 1.
    Input:
    Output: muestra (array)
    """
    dim = 2
    M = zeros((tam, dim))
    ef=0.0
    
    #Inical.
    M[0,0] = exponential(1) 
    M[0,1] = gamma(M[0,0], 1)
    
    for i in xrange(1,tam):
        a=M[i-1][0]
        l=M[i-1][1]
        numP = int(4*rand(1)[0])+1
        
        R = propuesta(numP,a,l,b,c,sigma,T)
        ap = R[0]
        lp = R[1]
        ro = R[2]
        
        if Bernoulli(ro) == 1.0:
            M[i,0] = ap
            M[i,1] = lp

        else:
            M[i,0] = M[i-1,0]
            M[i,1] = M[i-1,1]
            ef=ef+1
    
    print "Se rechazaron el "+ str(ef*100.0/tam)+ "% de las propuestas"
    return M
    
def bombasAguaMHMC(a,c,d,tam,w1):
    """
    Simula valores de la distribución posterior f(λ1, . . . , λn, β|p ̄), usando
    un kernel híbrido que considera las propuestas:
        λi ̄∼ Gama(t_ip_i + α , β + 1)
        β ∼ Gama(nα + γ , δ + Sum(λ)
    Con parámetros a priori dados
    
    Input: a (float), c(float), d (float), tam(int), w1(float)
    Output: array(tam x 11) (muestra)
    """
    D =  np.array([[94.32, 5],[15.72, 1],[62.88, 5],[125.76, 14],[5.24, 3],
                   [31.44, 19],[1.05, 1],[1.05, 1],[2.1, 4],[10.48, 22]])
    n=len(D)
    dim = n+1
    M = zeros((tam, dim))
    
    #Inical.
    M[0,n] = gamma(c, 1.0/(d)) #Beta a priori
    M[0,0:n] = gamma(a,1.0/(M[0,n]),n) #Lambdas a priori
    
    for i in xrange(1,tam):
        b= M[i-1,n]
        L=M[i-1,0:n]
        if Bernoulli(w1) == 1.0:
            #deja betas
            M[i,n]=b
            #Mueve lambdas
            for j in xrange(0,n):
                #M[i,j]= gamma(D[j,0]*D[j,1] + a, 1.0/(b+1) )
                M[i,j]= gamma(D[j,1] + a, 1.0/(b+D[j,0]) )
        else:
            #Mueve beta
            M[i,n]=gamma(n*a + c,1.0/(d+sum(L)))
            #Deja lambdas
            M[i,0:n] = L
    return M
    
def graficaColumnas(M):
    """
    Simula valores de la distribución posterior f(λ1, . . . , λn, β|p ̄), usando
    un kernel híbrido que considera las propuestas:
        λi ̄∼ Gama(t_ip_i + α , β + 1)
        β ∼ Gama(nα + γ , δ + Sum(λ)
    Con parámetros a priori dados
    
    Input: a (float), c(float), d (float), tam(int), w1(float)
    Output: array(tam x 11) (muestra)
    """
    r=4
    c=3
    f, axarr = plt.subplots(r, c,figsize=(8, 6), dpi=80)
 
    C=np.ones(len(M))
    for j in xrange(0,9):
        p= np.mean(M[10:,j])
        axarr[j/c, j%c].plot(M[10:,j],'o',markersize=2.5,alpha=0.3,color='#009999')
        axarr[j/c, j%c].plot(C*p,linewidth=3,color='#990033')
        axarr[j/c, j%c].set_title('lambda'+str(j+1))
        axarr[j/c, j%c].set_xticklabels([])
        print "El promedio de lambda"+str(j+1)+" es " + str(p)

    p= np.mean(M[10:,9])
    axarr[3, 1].plot(M[10:,9],'o',markersize=2.5,alpha=0.3,color='#009999')
    axarr[3, 1].plot(C*np.mean(M[10:,9]),linewidth=3,color='#990033')
    axarr[3, 1].set_title('lambda10')
    print "El promedio de lamda10 es " + str(p)
    
    #Bonito
    f.subplots_adjust(hspace=0.4)
    
    for t in [0,2]:
        axarr[3,t].spines['bottom'].set_color('white')
        axarr[3,t].spines['left'].set_color('white')
        axarr[3,t].spines['top'].set_color('white')
        axarr[3,t].spines['right'].set_color('white')
        for s in axarr[3,t].xaxis.get_ticklines(): s.set_color('white')
        for s in axarr[3,t].yaxis.get_ticklines(): s.set_color('white')

        plt.setp([axarr[3,t].get_yticklabels()], visible=False)
    
    for j in xrange(0,r):
        plt.setp([a.get_xticklabels() for a in axarr[j, :]], visible=False)

    plt.savefig('lambdas.png')
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
    1. Aplique el algoritmo de Metropolis-Hastings considerando como funcion 
    objetivo la distribución normal bivariada:
  
    Considere las siguientes propuestas:
        q1 ((x01, x02)|(x1, x2)) = f_X1|X2(x01|x2)1(x02 = x2)
        q2 ((x01, x02)|(x1, x2)) = fX2|X1(x02|x1)1(x01 = x1)

    A partir del algoritmo MH usando Kerneles híbridos simule valores de
    la distribución normal bivariada, fijando σ1 = σ2 = 1, considere los
    casos ρ = 0.8 y ρ = 0.99.
    """
    for n in [1000,10000]:
        m = np.array([[0],
                      [0]]) 
        w1=0.5
        #CASO 0.8
        ro = 0.8
        sigma = np.array([[1, ro],
                          [ro, 1]])
        
        graficaMuestraBi(m,sigma,w1,n)
        
        #CASO 0.99 
        ro = 0.99
        sigma = np.array([[1, ro],
                          [ro, 1]])

        graficaMuestraBi(m,sigma,w1,n)
            
    """
    2. Considere los tiempos de falla t1, . . . tn con distribución Weibull(α, λ):

        f(t|α, λ) = αλt^(α−1)e^(−tαiλ)

    Se asumen como a priori α ∼ exp(c) y λ|α ∼ Gama(α, b), por lo tanto,
    f(α, λ) = f(λ|α)f(α). Así, para la disitribución posterior se tiene:
   
        f(α, λ|t ̄) ∝ f(t ̄|α, λ)f(α, λ)

    A partir del algoritmo MH usando Kerneles híbridos simule valores de
    la distribución posterior f(α, λ|t ̄), considerando las siguientes propues-
    tas:

    Propuesta 1: 
        λp|α,t ̄∼ Gamaα + n , b +Xni=1tαi! y dejando α fijo.
    Propuesta 2:
        αp|λ,t ̄ ∼ Gama (n + 1 , −log(b) − log(r1) + c) , y dejando λ fijo.
    Propuesta 3:
        αp ∼ exp(c) y λp|αp ∼ Gama(αp, b).
    Propuesta 4 (RWMH):
        αp = α + , con  ∼ N(0, σ) y dejando λ fijo.

    Simular datos usando α = 1 y λ = 1 con n = 20. Para la a priori usar    
    c = 1 y b = 1.
    """
    
    n=1000
    T = weibull(1,20)
    M = TiemposMHMC(1,1,0.1,n,T)
    grafTiempos(M,T,1,1)
 
    A= M[:,0]
    B= M[:,1]
    x = (A).tolist()
    y = (B).tolist()

    #Scatter    
    colors = np.arange(0, 1, 1.0/n)
    area = 50*np.ones(n)
    plt.scatter(x, y, s=area, c=colors, alpha=0.5)
    plt.savefig('tiemposMH'+str(n)+'.png')
    plt.title('Muestra de tamano '+str(n))
    plt.show()
    
    
    """
    3. Considere el ejemplo referente al n ́umero de fallas de bombas de agua
    en una central nuclear, donde pi representa el n ́umero de fallas en el
    tiempo de operaci ́on ti, con i = 1, . . . n.

    Se considera el modelo pi ∼ P oisson(λiti), (las λi son independien-
    tes entre si), con distribuciones a priori λi|β ∼ Gama(α, β) y β ∼Gama(γ, δ),
    por lo tanto:

        f(λ1, . . . , λn, β) = f(λ1|β)f(λ2|β). . . f(λn|β)f(β)

    Para la distribución posterior se tiene:
        
        f(λ1, . . . , λn, β|p ̄) ∝ L( ̄p, λ, β  ̄ )f(λ1, . . . , λn, β)

    Simule valores de la distribuci ́on posterior f(λ1, . . . , λn, β|p ̄), usando
    un kernel h ́ıbrido, considerando las propuestas:

        λi|λ ̄−i, β,t ̄ ∼ Gama(tipi + α , β + 1)β|λ,  ̄ t ̄ ∼ Gama

        nα + γ , δ +Xni=1λi!

    Verifique que estas son propuestas Gibbs.

    Use los datos del Cuadro 1 con los parámetros a priori α = 1.8, γ = 0.01
    y δ = 1.
    """
    tam=1000
    a=1.8
    c=0.01
    d=1
    w1=0.5
    M = bombasAguaMHMC(a,c,d,tam,w1)    
    graficaColumnas(M)
    
    C=np.ones(len(M))
    plt.figure(num=None, figsize=(3, 3), dpi=80)
    plt.plot(M[10:,10],'o',markersize=5,alpha=0.3,color='#ffcc00')
    p = np.mean(M[10:,10])
    plt.plot(C*p,linewidth=3,color='#990033')
    plt.title('Betas simuladas con Gibs sampler')
    plt.savefig('betas.png')
    plt.show()
    print "El promedio de beta es " + str(p)

    lambdas= np.zeros(11)
    for j in xrange(0,11):
        lambdas[j]=np.mean(M[10:,j])