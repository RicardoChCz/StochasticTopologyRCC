# -*- coding: utf-8 -*-
"""
Created on Thu Sep 07 22:28:09 2017

@author: Ricardo Chávez Cáliz
"""
import numpy as np
from time import clock
from numpy import sqrt
from numpy import pi
from numpy import sin
from numpy.linalg import qr
import matplotlib.pyplot as plt

#Tolerencia
t = 0.001
def VerificaTS(M):
    """
    Verifica si la matriz M es Triangular superior.
    Input:  Una matriz de n x n.
    Output: -
    """
    n = len(M)
    for k in xrange(0,n):
        for j in xrange(0,k):
            if M[k,j] > t:
                raise ValueError("La matriz no es triangular superior")
                return False
            else:
                return True
                
def BackwardSubst (tS, b):
    """
    Backward Substitution.
    Input:  Una matriz tS de n x n triangular superior.
            y b un vector de n x 1.
    Output: Un vector X que solucione tS*X=b.
    """
    n = len(tS)
    if b.size != n:
        raise ValueError("Las dimensiones de tS y b no son compatibles,", n, b.size)    
    if VerificaTS(tS) == True:   
        X = np.zeros((n,1))
        #Paso inicial
        X[n-1,0] = float(b[n-1,0])/tS[n-1, n-1]

        #Pal real
        for i in range (1, n):
            k = n-i-1
            suma = 0.0
            for i in range (k+1, n):
                suma = suma + float(tS[k,i])*X[i,0]
            X[k,0] = float(b[k,0]- suma)/tS[k,k]
    else:
        raise ValueError("La matriz no es triangular superior")   

    return X


def norma(V):
    """
    Calcula la norma de V
    Input: Un array de n*1
    Output: Un real.
    """
    N = 0.
    n=len(V)
    for i in xrange(0,n):
        N = N + V[i]**2

    return sqrt(N)

def gs(A):
    """
    Algoritmo modificado de Gram-Schmidt para obtener vectores
    q1,...qn a partir de una matriz A con q_j pertenezca a lo
    generado por a1.... a_j y sea ortonormal a los anteriores
    Input: Una matriz A de rango completo de tamaño mxn
    Output: Matrices Q y R tales que A=QR donde las columnas
            de Q son como se dice en la descripcion (Q*Q=I) y R es 
            triangular sup    
    """
    m = len(A)
    n = len(A.T)
    if m >= n:
        Q =  np.zeros((m, n))
        R =  np.zeros((n, n)) 
        for j in xrange(n):
            v = A[:,j]
            for i in xrange(j):
                R[i,j] =  np.dot(Q[:,i].T , A[:,j])
                v = v.squeeze() - (R[i,j] * Q[:,i])
            R[j,j] =  np.linalg.norm(v)
            Q[:,j] = (v / R[j,j]).squeeze()
        return Q, R
    else:
        raise ValueError("Las dimensiones no satisfacen m>=n")   

def construyeX(n):
    """
    Construye un conjunto de observaciones de tamaño n
    Input: int n: Tamaño de la muestra
    Output: lista tamaño n
    """
    X=np.zeros(n)
    for i in xrange(0,n):
        X[i] = 4*pi*i/n
    return X
    
def construyeY(X, sigma):
    """
    Construye un vector aleatorio Y con forma sen(x) más un ruido
    gaussiano de media 0 y varianza sigma, partiendo de X.
    Input: Intervalo discretizado (array, tam n)
    Output: array tamaño n*1 con forma sen(x) + N(0,sigma)
    """
    n=len(X)
    e=np.zeros(n)
    for i in xrange(0,n):
        e[i] = np.random.normal(0, sigma)
    
    Y = np.zeros((n,1))
    for i in xrange(0, n):
        Y[i,0] = sin(X[i]) + e[i]
        
    return Y

def matrizDeDiseno(X, p):
    """
    Construye la matriz de diseño con observaciones X mediante 
    un polinomio de orden p-1.
    Input: X array de tamaño n y p orden aumentado en uno del polinomio
    Output:  Matriz de tamaño n*p
    """
    n=len(X)
    M = np.zeros((n, p))
    
    M[:,0] = 1
    
    #En cada renglon
    for i in xrange(0,n):
        for j in xrange(1,p):
            M[i,j] = X[i]**(j)

    return M

def estimadorMC(X, Y, p, artesanal=True):
    """
    Halla b tal que ||Y-Mb||_2 sea mín. Donde M es matriz de diseño obtenida 
    del vector de observaciones X y Y es el vector aleatorio a estimar. 
    Para esto llama al método que construye la matriz de diseño de X,
    y al método que calcula descomposición QR de esta, resuelve b en Rb=Q*Y 
    usando Backward sustitution y devuelve Mb.
    
    Input: X array de tamaño n= número de observaciones y Y aleatorio a estimar, 
            parametro de estimación p (grado + 1 de la matriz de diseño) y boolean
            que determina que método de descomposición QR usar. True para el descrito
            en este código y False para el propio de Numpy
    Output:  b array estimador que minimiza ||Y-Xb||
    """    
    M = matrizDeDiseno(X,p)
    
    if (artesanal):
        S = gs(M)

    else:
        S = np.linalg.qr(M)
        
    Q = S[0]
    R = S[1]        
    
    B = BackwardSubst (R, (np.dot(Q.T,Y)))

    return np.dot(M,B)

def mediana(A):
    O=sorted(A)

    n=len(O)
    m=n/2
    med=0
    
    # codigo para calcular la mediana
    if n%2==0:
        med = (O[m+1] + O[m+2]) / 2.0
    else:
        med = O[m+1]
        
    return med
    
if __name__ == "__main__":    

    """
    #1. Implementar el algoritmo de Gram-Schmidt modificado 8.1 del Trefethen
    #(p. 58) para generar la descomposición QR.
    """
    m=5
    n=4
    for i in xrange(10):
        A = np.random.rand(m,n)
        S = gs(A)   
        Q = S[0]
        R = S[1]
        
        if np.allclose(A, np.dot(Q, R)) == True and np.allclose(np.dot(Q.T, Q) , np.identity(n)) == True and VerificaTS(R) == True:
            print "Descomposición exitosa"
    
        
    """
    #2. Implementar el algoritmo que calcula el estimador de mínimos cuadrados
    #en una regresión usando la descomposición QR.
    
    Algortimo se encuentra en la función estimadorMC
    """
    
    """
    #3. Generar Y compuesto de yi = sen(xi) + e_i donde e_i ~ N(0, σ) con
    #σ = 0.1, para xi = 4πin para i = 1, . . . , n.
    #Hacer un ajuste de mínimos cuadrados a Y, con descomposición QR,
    #ajustando un polinomio de grado p − 1.
    #• Considerar los 12 casos: p = 3, 4, 5, 100 y n = 100, 1000, 10000.
    #• Graficar el ajuste en cada caso.
    #• Medir tiempo de ejecución de su algoritmo, comparar con descomposición
    #QR de scipy y graficar los resultados.
    """
    tama = [100,1000,10000]
    p = [3,4,5,100]
    g= len(p)

    for n in tama:
        X = construyeX(n)
        Y = construyeY(X,0.1)
        fig0 = plt.figure(figsize = (16,16))
        fig0.subplots_adjust(hspace=0.5)
        f, axarr = plt.subplots(2, 2)
        axarr[0, 0].plot(X, Y, 'bo', markersize=1)
        axarr[0, 1].plot(X, Y, 'bo', markersize=1)
        axarr[1, 0].plot(X, Y, 'bo', markersize=1)
        axarr[1, 1].plot(X, Y, 'bo', markersize=1)

        for i in xrange(0,g):
            E = estimadorMC(X, Y, p[i], True)
            t=i%2
            s=i/2
            axarr[s, t].plot(X, E, 'r')
            axarr[s, t].set_title('Para n=' + str(n) + " y p=" + str(p[i]))
            
        plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
        plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
        plt.savefig("MinCuad" + str(n))
        plt.show()

    #Uno con n = 1000 y con p=10

    fig = plt.figure()
    fig.subplots_adjust(bottom=0.025, left=0.025, top = 0.975, right=0.975)

    X = construyeX(100)
    Y = construyeY(X,0.1)
    E = estimadorMC(X, Y, 10, True)
    
    plt.subplot(2, 3, 4)
    plt.plot(X, Y, 'bo', markersize=1)
    plt.plot(X, E, 'r')
    
    X = construyeX(1000)
    Y = construyeY(X,0.1)
    E = estimadorMC(X, Y, 10, True)
    
    plt.subplot(2, 3, 5)
    plt.plot(X, Y, 'bo', markersize=1)
    plt.plot(X, E, 'r')

    X = construyeX(10000)
    Y = construyeY(X,0.1)
    E = estimadorMC(X, Y, 10, True)
    plt.subplot(2, 3, 6)
    plt.plot(X, Y, 'bo', markersize=1)
    plt.plot(X, E, 'r')

    plt.show()
        
    c = 100
    x = np.arange(0, c, 1);
    y1= np.zeros(c)
    y2= np.zeros(c)
    
    for k in xrange(2, c):
        X = construyeX(k)
        Y = construyeY(X,0.1)
        
        tiempo_inicial = clock()  
        E1 = estimadorMC(X, Y, 2, artesanal=True)
        y1[k] = clock() - tiempo_inicial
 
        tiempo_inicial = clock()  
        E2 = estimadorMC(X, Y, 2, artesanal=False)    
    
        y2[k] = clock() - tiempo_inicial
        
    plt.plot(x, y1, 'r', label="Algortimo propio")
    plt.plot(x, y2, 'b--', label="Algoritmo scipy")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("Compara")
    plt.show()
    
    #Tiempos medios
    c = 100
    r = 10
    x = np.arange(0, c, 1);
    y1= np.zeros(c)
    y2= np.zeros(c)
    
    for i in xrange(2, c):
        Art = [0]
        Num = [0]
        for k in xrange(1, r):
            X = construyeX(c)
            Y = construyeY(X,0.1)
        
            tiempo_inicial = clock()  
            E1 = estimadorMC(X, Y, 2, artesanal=True)
            a= clock() - tiempo_inicial
            Art.append(a)
 
            tiempo_inicial = clock()  
            E2 = estimadorMC(X, Y, 2, artesanal=False)    
            b= clock() - tiempo_inicial
            Num.append(b)
            
        y1[i]=mediana(Art)
        y2[i]=mediana(Num)
        
    plt.plot(x, y1, 'r', label="Algortimo propio")
    plt.plot(x, y2, 'b--', label="Algoritmo scipy")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()

    
    """#
    #4. Hacer p = 0.1n, o sea, diez veces más observaciones que coeficientes en
    #la regresión, ¿Cual es la n máxima que puede manejar su computadora?
    #"""
    """
    for n in xrange(10128, 1000000):
        p=int(0.1*n)
        tiempo_inicial = clock()  
        X = construyeX(n)
        Y = construyeY(X,0.1)
        E1 = estimadorMC(X, Y, p, artesanal=True)
        b= clock() - tiempo_inicial
        print b
        print n
    """