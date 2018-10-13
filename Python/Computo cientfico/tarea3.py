# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 22:25:06 2017

@author: Ricardo Chávez Cáliz
"""
import numpy as np
import scipy as sc
from time import clock
from numpy import sqrt
from numpy.linalg import qr
import matplotlib.pyplot as plt
from numpy import log
from scipy.linalg import cholesky

#Tolerencia
t = 0.001

def construyeDiag(S):
    m = len(S)
    D =  np.identity(m)
    for i in xrange(0,m):
        D[i,i] = S[i]
    
    return D
                
def Cholesky(A):
    """
    Factorización Cholesky.
    Input:  Una matriz A de n x n hermitiana definida postiva.
    Output: Una matriz R tal que A=R*R^{*} 
    """
    n =  len(A)
    R = np.copy(A)
 
    for k in xrange(0,n):
        for j in xrange(k+1, n):
            R[j,j:n] = R[j,j:n] - R[k,j:n]*R[k,j]/R[k,k]
        R[k,k:n] = R[k,k:n]/np.sqrt(abs(R[k,k]))
    
    for k in xrange(1,n):
        for j in xrange(0, k):
            R[k,j]= 0.0
    
    return R

def comparaCholesky(S, artesanal=True, diferencia=True):
    """
    Compara las descomposiciones Cholesky de B y B(epsilon), permitiendo
    elegir la descomposición de Cholesky preferida (artesanal o la de scipy)
    Input:  EigenValores
    Output: Distancia en supremo entre B y B(epsilon)
    """
    D = construyeDiag(S)
    B = np.dot(np.dot(Q.T,D),Q)

    #Perturbar B    
    sigma=0.01*S[0]
    S = S + np.random.normal(0,sigma,m)
    D = construyeDiag(S)
    Be = np.dot(np.dot(Q.T,D),Q)
    
    if (artesanal):
        U=Cholesky(B)
        Ue=Cholesky(Be)

    else:
        U=cholesky(B)
        Ue=cholesky(Be)

    if(diferencia):
        return normaEuc(U - Ue)
    else: 
        return (normaEuc(U - Ue)*normaEuc(B))/(normaEuc(U)*normaEuc(B-Be))

    
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

def BackwardSubst (tS, b):
    """
    Backward Substitution.
    Input:  Una matriz tS de n x n triangular superior.
            y b un vector de n x 1.
    Output: Un vector X que solucione tS*X=b.
    """
    n = len(tS)
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
    
    return X


def estimadorMC(X, Y):
    """
    Halla b tal que ||Y-Mb||_2 sea mín. Donde M es matriz de diseño obtenida 
    del vector de observaciones X y Y es el vector aleatorio a estimar. 
    Para esto llama al método que construye la matriz de diseño de X,
    y al método que calcula descomposición QR de esta, resuelve b en Rb=Q*Y 
    usando Backward sustitution.
    
    Input: X array de tamaño n= número de observaciones y Y aleatorio a estimar, 
            parametro de estimación p (grado + 1 de la matriz de diseño) y boolean
            que determina que método de descomposición QR usar. True para el descrito
            en este código y False para el propio de Numpy
    Output:  b array estimador que minimiza ||Y-Xb||
    """
    S = gs(X)        
    Q = S[0]
    R = S[1]        

    B = BackwardSubst (R, (np.dot(Q.T,Y)))

    return B
    
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

def normaEuc(A):
    """
    Halla la distancia euclideana entre A y B (matrices)
    Input: Matrices A y B
    Output:  float
    """
    m = len (A)
    n = len (A.T)
    
    dist = 0.0
    for i in xrange(0,m):
        for j in xrange(0,n):
            dist= dist + (A[i,j])**2
                
    return sqrt(dist)
            
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
    1. Sea A una matriz de tamañoo 50 × 20, aleatoria y fija, calcular su 
    descomposición QR. Sean λ1 ≥ λ2 ≥ ... ≥ λ20 ≥ 0 y
    B = Q∗diag(λ1, λ2, ..., λ20)Q y Bε = Q∗diag(λ1+ε1, λ2+ε2, ..., λ20+ε20)Q,
    con εi ∼ N(0, σ), con σ = 0.01λ1.
    
    (a) Comparar la descomposicón de Cholesky de B y de Bε usando el
    algoritmo de la tarea 1. Considerar los casos cuando B tiene un
    buen número de condición y un mal número de condición.
    """
    m=20
    n=50
    A = np.random.rand(m,n)
    Q = qr(A)[0]
    

    #Aumentar el malcondicionamiento, y tomar ||U-U_{e}||
    x = np.arange(10, 5000, 10);
    N = len(x)
    y1= np.zeros(N)
    y2= np.zeros(N)
    
    for i in xrange(0,N):
        E = np.ones(m)
        y1[i] = comparaCholesky(E, artesanal=True, diferencia=True)
        E[0] = x[i]
        y2[i] = comparaCholesky(E, artesanal=True, diferencia=True)
        
    plt.xlabel('valor del primer eigenvalor')
    plt.ylabel('norma de la diferencia')
    plt.plot(x, y1, 'go', label="Bien condicionado")
    plt.plot(x, y2, 'r^', label="Mal condicionado")    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("comparaBCvsMC")
    plt.show()
    
    #Aumentar el malcondicionamiento, y tomar log(||U-U_{e}||)

    for i in xrange(0,N):
        y1[i]=log(y1[i])
        y2[i]=log(y2[i])
    plt.xlabel('valor del primer eigenvalor')
    plt.ylabel('log de la norma de la diferencia')
    plt.plot(x, y1, 'go', label="Bien condicionado")
    plt.plot(x, y2, 'r^', label="Mal condicionado")    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("comparaBCvsMC_log")
    plt.show()
    
    #Aumentar el malcondicionamiento, y tomar mediana de ||U-U_{e}||
    rep=30
    x = np.arange(10, 1000, 10);
    N = len(x)
    y1= np.zeros(N)
    y2= np.zeros(N)
    
    for i in xrange(0,N):
        BC= [0]
        MC= [0]
        for j in xrange(1,rep):            
            E = np.ones(m)
            BC.append(comparaCholesky(E, artesanal=True, diferencia=True))
            E[0] = x[i]
            MC.append(comparaCholesky(E, artesanal=True, diferencia=True))
        
        y1[i] = mediana(BC)
        y2[i] = mediana(MC)

    plt.xlabel('valor del primer eigenvalor')
    plt.ylabel('mediana de norma de la diferencia')
    plt.plot(x, y1, 'go', label="Bien condicionado")
    plt.plot(x, y2, 'r^', label="Mal condicionado")    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("comparaMBCvsMMC_mediana")
    plt.show()
    
    #Aumentar el malcondicionamiento, y tomar log mediana de ||U-U_{e}||
    for i in xrange(0,N):
        y1[i] = log(y1[i])
        y2[i] = log(y2[i])

    plt.xlabel('valor del primer eigenvalor')
    plt.ylabel('log de norma media de la diferencia')
    plt.plot(x, y1, 'go', label="Bien condicionado")
    plt.plot(x, y2, 'r^', label="Mal condicionado")    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("comparaMBCvsMMC_mediana_log")
    plt.show()
    
    #Aumentar el malcondicionamiento, y tomar mediana de c.
    y1= np.zeros(N)
    y2= np.zeros(N)
    
    for i in xrange(0,N):
        BC= [0]
        MC= [0]
        for j in xrange(1,rep):            
            E = np.ones(m)
            BC.append(comparaCholesky(E, artesanal=True, diferencia=False))
            E[0] = x[i]
            MC.append(comparaCholesky(E, artesanal=True, diferencia=False))
        
        y1[i] = mediana(BC)
        y2[i] = mediana(MC)

    plt.xlabel('valor del primer eigenvalor')
    plt.ylabel('mediana de la alteracion al problema de Cholesky')
    plt.plot(x, y1, 'go', label="Bien condicionado")
    plt.plot(x, y2, 'r^', label="Mal condicionado")    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("comparaABCvsAMC_mediana")
    plt.show()
    
    """
    (b) Con el caso mal condicionado, comparar el resultado de su algoritmo
    con el del algoritmo de Cholesky de scipy
    """   
    pro=0
    rep=100
    for j in xrange(0,rep):
        for i in xrange(1,100):
            E = np.ones(m)
            E[0] = i
            try:
                comparaCholesky(E, artesanal=False, diferencia=True)
            except:
                pro= pro+i
                break
            
    print pro*(1.0/rep)
    """
    (c) Medir el tiempo de ejecución de su algoritmo de Cholesky con el
    de scipy.    
    """
    c = 200
    z = 1
    x = np.arange(z, c+z, 1);
    y1= np.zeros(c)
    y2= np.zeros(c)

    for k in xrange(z, c+z):
        A = np.random.rand(k,k)
        B = np.dot(A.T, A) + np.identity(k)
        
        tiempo_inicial = clock()
        Cholesky(B)
        y1[k-z] = log(clock() - tiempo_inicial)
 
        tiempo_inicial = clock()
        sc.linalg.cholesky(B)
        y2[k-z] = log(clock() - tiempo_inicial)
    
    plt.xlabel('tamano de la matriz')
    plt.ylabel('log(tiempo de ejecucion)')
    plt.plot(x, y1, 'r', label="Algortimo propio")
    plt.plot(x, y2, 'b--', label="Algoritmo scipy")    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig("Compara_tiempos_Cholesky_t3")
    plt.show()

    

    """
    2. Resolver el problema de mínimos cuadrados, y = Xβ + ε, εi ∼ N(0, σ)
    usando su implementación de la descomposición QR; β es de tamaño
    d × 1 y X de tamaño n × d.
    Sean d = 5, n = 20, β = (5, 4, 3, 2, 1)^t y σ = 0.1.
    
    (a) Hacer X con entradas aleatorias U(0, 1) y simular y. Encontrar
    βˆ y compararlo con el obtenido βˆp haciendo X + ∆X, donde las
    entradas de ∆X son N(0, σ = 0.01). Comparar a su vez con 
    βˆc = ((X + ∆X)^t (X + ∆X))−1 (X + ∆X)y usando el algoritmo
    genérico para invertir matrices scipy.linalg.inv
    """
    n=20
    d=5
    sigma = 0.1

    rep = 1000
    prom1 = 0
    prom2 = 0
    prom3 = 0

    for i in xrange(0,rep):
        X = np.random.rand(n,d)
        b = (d+1 - np.arange(1,d+1)).reshape((d,1))
        Y=np.dot(X,b)+ np.random.normal(0,sigma,n).reshape((n,1))
    
        #beta gorro
        bg = estimadorMC(X, Y)
        
        #beta perturbado
        DX= np.zeros((n,d))
        for i in xrange (0,n):
            DX[i,0:]= np.random.normal(0, 0.01, d)    
        Y=np.dot(X+DX,b)+ np.random.normal(0,sigma,n).reshape((n,1))
        bp = estimadorMC(X+DX, Y)
        #beta c
        bc = np.dot(np.dot(np.linalg.inv(np.dot((X+DX).T, (X+DX))), (X+DX).T),Y)        

        prom1 = prom1 + normaEuc(bg-bp)
        prom2 = prom2 + normaEuc(bp-bc)
        
        prom3 = prom3 + (normaEuc(bg-bp)*normaEuc(X))/(normaEuc(bg)*normaEuc(DX))
    
    print "el promedio de ||b̃ - b̃_p||_{2} es " + str(prom1*(1.0/rep)) 
    print "representa un error de " + str(prom1*(1.0/rep)/5) + " respecto a b"
    print "el promedio de ||b̃_p - b̃_c||_{2} es " + str(prom2*(1.0/rep)) 
    print "representa un error de " + str(prom2*(1.0/rep)/5) + " respecto a b"
    print "el promedio de c_{f} es " + str(prom3*(1.0/rep)) 
    
    
    
    """
    (b) Lo mismo que el anterior pero con X mal condicionada (ie. con
    casi colinealidad)
    """
    prom1 = 0
    prom2 = 0
    prom3 = 0

    for i in xrange(0,rep):
        X = np.random.rand(n,d)
        X[0:,4] = X[0:,3] + 0.01
        b = (d+1 - np.arange(1,d+1)).reshape((d,1))
        Y=np.dot(X,b)+ np.random.normal(0,sigma,n).reshape((n,1))
    
        #beta gorro
        bg = estimadorMC(X, Y)
        
        #beta perturbado
        DX= np.zeros((n,d))
        for i in xrange (0,n):
            DX[i,0:]= np.random.normal(0, 0.01, d)    
        Y=np.dot(X+DX,b)+ np.random.normal(0,sigma,n).reshape((n,1))
        bp = estimadorMC(X+DX, Y)
        #beta c
        bc = np.dot(np.dot(np.linalg.inv(np.dot((X+DX).T, (X+DX))), (X+DX).T),Y)        

        prom1 = prom1 + normaEuc(bg-bp)
        prom2 = prom2 + normaEuc(bp-bc)
        
        prom3 = prom3 + (normaEuc(bg-bp)*normaEuc(X))/(normaEuc(bg)*normaEuc(DX))
    
    print "el promedio de ||b̃ - b̃_p||_{2} es " + str(prom1*(1.0/rep)) 
    print "representa un error de " + str(prom1*(1.0/rep)/5) + " respecto a b"
    print "el promedio de ||b̃_p - b̃_c||_{2} es " + str(prom2*(1.0/rep)) 
    print "representa un error de " + str(prom2*(1.0/rep)/5) + " respecto a b"
    print "el promedio de c_{f} es " + str(prom3*(1.0/rep)) 