# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 16:11:03 2017

@author: Ricardo Chávez Cáliz
"""
import numpy as np
from time import time
import matplotlib.pyplot as plt

#Tolerencia
t = 0.001

def SistemasAleatorios(D):
    """
    Genera sistemas de ecuaciones aleatorios Vector solucion ~ U(0,1)
    Input:  Una matriz de n x n.
    Output: -
    """
    if abs(np.linalg.det(D)) < t:
        raise ValueError("D es singular")
    b = np.random.rand(n,1)
    print("Vector aleatorio:")
    print b
    print("Solución")
    x = Resuelve(D, b)
    print x
    print("Comprobación")
    print np.dot(D, x) - b

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
    VerificaTS(tS)
        
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

def ForwardSubst (tI, b):
    """
    Forward Substitution.
    Input:  Una matriz tI de n x n triangular inferior.
            y b un vector de n x 1.
    Output: Un vector X de n x 1 que solucione tI*X=b.
    """
    n = len(tI)
    if b.size != n:
        raise ValueError("Las dimensiones de tS y b no son compatibles,", n, b.size)
    
    VerificaTS(np.transpose(tI))
    X = np.zeros((n,1))    
    #Paso inicial
    X[0,0] = b[0,0]/tI[0,0]

    #Pal real
    for k in range (1, n):
        suma = 0
        for j in range (0, k):
            suma = suma + tI[k,j]*X[j,0]
            X[k,0] = (b[k,0]- suma)/tI[k,k]
    return X


def EGPP(A):
    """
    Eliminación Gaussiana con pivoteo parcial, donde se consigue la
    descomposión LUP de una matriz A.
    Input:  Una matriz U de n x n no singular.
    Output: *U matriz triangular superior a la que le fue aplicado
            el algoritmo de eliminación Gaussiana con pivoteo parcial.
            *P matriz de permutacion P
            *L matriz triangular inferior
    """
    n =  len(A)
    if abs(np.linalg.det(A)) < t:
        raise ValueError("A es singular")
    
    #Inicializamos matriz de permutación y la triangular inferior
    U = np.copy(A)
    L = np.identity(n)
    P = np.identity(n)
    
    # k representa el pivote renglón en el que estás (k tambien denotará la 
    # columna k-ésima de la diagonal).
    for k in xrange(0,n):
        """
        #Para estabilidad:
        #Selecciona el indice i (i>= k) para maximizar |u_{ik}|
        #Busca quien de los que están en la columna de u_{kk}
        #por debajo es el más grande
        """
        indMax = abs(U[k:,k]).argmax() + k
        #Verifica que la matriz sea de las buenas
        if U[indMax, k] == 0:
            raise ValueError("La matriz es singular")
            
        """    
        # Intercambia dos renglones para dejar como pivote al más grande
        # para evitar divisiones cercanas a 0. Le da estabilidad al algoritmo
        # y no aumenta la complejidad.
        #  Intercambia tambien los renglones de P y L
        
        guarda = np.copy(U[k,k:n])
        U[k,k:n] = U[indMax, k:n]        
        U[indMax, k:n] = guarda
            
        guarda = np.copy(L[k,1:k-1])
        L[k,1:k-1] = L[indMax,1:k-1]
        L[indMax,1:k-1] = guarda
            
        guarda = np.copy(P[k,0:])
        P[k,0:] = P[indMax,0:]
        P[indMax,0:] = guarda
        """

        #Para cada renglón            
        for j in xrange(k+1, n):
            L[j,k] = U[j,k]/U[k,k]
            U[j,k:n] = U[j,k:n] - (L[j,k]*U[k,k:n])
    
    return (P,L,U)

def Cholesky(A):
    """
    Factorización Cholesky.
    Input:  Una matriz A de n x n hermitiana definida postiva.
    Output: Una matriz R tal que A=R*R^{*} 
    """
    n =  len(A)
    R = np.copy(A)
    if abs(np.linalg.det(A)) < t:
        raise ValueError("A es singular")
 
    for k in xrange(0,n):
        for j in xrange(k+1, n):
            R[j,j:n] = R[j,j:n] - R[k,j:n]*R[k,j]/R[k,k]
        R[k,k:n] = R[k,k:n]/np.sqrt(abs(R[k,k]))
    
    for k in xrange(1,n):
        for j in xrange(0, k):
            R[k,j]= 0.0
    
    return R

def ConstruyeA(n):
    """
    Función constructora para una matriz con -1's abajo de la diagonal,
    1's en diagonal y 1's en la última columna)
    Input:  n entero tamaño de la matriz.
    Output: Una matriz de n*n como se indicó 
    """
    #ifn entero other wise no hags nada
    A = np.zeros((n,n))
    
    for k in xrange(0,n):
        for j in xrange(0,k):
            A[k,j] = -1
        A[k,k] = 1
        A[k,n-1] = 1
    return A

def ConstruyeTS(n):
    """
    Función constructora para una matriz aleatoria triangular superior
    de entradas U(0,1) de tamaño nxn
    Input:  n entero tamaño de la matriz.
    Output: Una matriz de nxn como se indicó 
    """
    #ifn entero other wise no hags nada
    A = np.zeros((n,n))
    
    for k in xrange(0,n):
        for j in xrange(k,n):
            A[k,j] = np.random.rand()
    return A
    
def Imprime(T):
    """
    Función que imprime información sobre la descomposición LUP de la
    matriz M.
    Input:  Arreglo proviniente de EGPP(A)=[P,L,U]
    Output: - 
    """
    print ("P matriz de permutación")
    print T[0]
    print ("L matriz triangular inferior")
    print T[1]
    print ("U matriz triangular superior")
    print T[2]
    
def Verificar (L,U,P):
    """
    Input: Matrices L,U,P provinientes de descomposicón de A
    Output: Producto de matrices LUP
    """
    return np.dot(P, np.dot(L, U))
    

def Resuelve(D, b):
    """
    Resuelva sistema de ecuaciones usando la descomposición LUP y los
    métodos de Backward and Forward substitution.
    Input:  Una matriz D de n*n no singular y b un vector de solución 
            de tamaño n*1.
    Output: Una matriz x de n*n con los valores que satisfacen el
            el sistema Dx=b.
    """
    #Revise si D es singular, si b tiene las mismas entradas
    n = len(D)
    if b.size != n:
        raise ValueError("Las dimensiones de la matriz y b no son compatibles,", n, b.size)
    if abs(np.linalg.det(A)) < t:
        raise ValueError("A es singular")
    #Obtener descomposición LUP
    T = EGPP(D)
    L = T[1]
    U = T[2]
    P = T[0]
    
    # Sistema es Ax=b => P'Ax=P'b (P' inversa de P, su transpuesta)
    # como A = LUP ent P'A=LU. Ahora tenemos el sistema LUx=P'b = c
    # si y=Ux entonces Ly=c (podemos aplicar FS) ahora y será conocida
    #y podemos resolver para Ux=y usando BS.
    c = np.dot(np.transpose(P),b)
    y= ForwardSubst(L,c)
    x= BackwardSubst (U, y)
    return x
    
if __name__ == "__main__":
    
    """
    1. Implementación de algoritmos Backward y Forward substitution.
    """
    n=5
    tS = ConstruyeTS(n)
    tI = tS.T
    b = np.random.rand(n,1)
    print ("1. Verificación de algoritmos BS y FS. Las siguientes dos entradas deben aproximarse a 0 en R^{n}")
    x = BackwardSubst(tS,b) #Debe resolver el sistema tS*x=b
    print np.dot(tS,x)-b
    x = ForwardSubst(tI,b) #Debe resolver el sistema tI*x=b
    print np.dot(tI,x)-b
    
    """
    2. Implementar el algoritmo de eliminación gaussiana con pivotteo
    parcial LUP, 21.1 del Trefethen (p. 160).
    
    Dicho algoritmo aparece definido en EGPP. Salvo que la parte del pivoteo
    parcial aparece comentada porque por alguna razón no modificaba de manera
    adecuada a L, y resultaba una descomposión incorrecta.
    """
    
    #Matrices que serán usadas en el problema 3
    n = 5 
    D = np.random.rand(n,n)
    A = ConstruyeA(n)
    S = EGPP(A)
    T = EGPP(D)

    print ("3. Verificación del algoritmo LUP. Las siguientes dos matrices deben aproximarse a 0")
    print Verificar(S[1],S[2],S[0]) - A
    print Verificar(T[1],T[2],T[0]) - D


    """
    3. Dar la descomposición LUP para una matriz aleatoria de entradas
    U(0, 1) de tamaño 5 × 5, y para la matriz específica A.
    """
    
    print("Para la matriz A,")
    print (A)
    Imprime(S)
    
    print("Para la matriz aleatoria D,")
    print (D)
    Imprime(T)
    
        
    """
    4. Usando la descomposición LUP anterior, resolver el sistema de la forma
    Dx = b
    donde D son las matrices del problema 3, para 5 diferentes b aleatorios
    con entradas U(0, 1). Verificando si es o no posible resolver el sistema.
    """ 
    for i in xrange(1,6):
        print("Prueba "+ str(i))
        SistemasAleatorios(D)
        
    """
    5. Implementar el algoritmo de descomposición de Cholesky 23.1 del Trefethen
    (p. 175).
    """
    n = 5 
    A = np.random.rand(n,n)
    B = np.dot(A.T, A) + np.identity(n)
    print ("5. Verificación del algoritmo Choleski. Las siguiente matriz debe aproximarse a 0")
    L = Cholesky(B)
    print np.dot(L.T.conj(), L)-B
    
    
                
    """
    6. Comparar la complejidad de su implementación de los algoritmos de
    factorización de Cholesky y LUP mediante la medición de los tiempos
    que tardan con respecto a la descomposición de una matriz aleatoria
    hermitiana definida positiva. Graficar la comparación.
    c = 200
    x = np.arange(0, c, 1);
    y1= np.zeros(c)
    y2= np.zeros(c)
    
    for k in xrange (0, c):
        A = np.random.rand(k,k)
        B = np.dot(A.T, A) + np.identity(k)
        
        tiempo_inicial = time()  
        Cholesky(B)    
        tiempo_final = time() 
        y1[k] = tiempo_final - tiempo_inicial
        
        tiempo_inicial = time()  
        EGPP(B)    
        tiempo_final = time() 
    
        y2[k] = tiempo_final - tiempo_inicial
        
    plt.plot(x, y1, 'r', x, y2, 'b')
    plt.savefig("Graf1")
    """