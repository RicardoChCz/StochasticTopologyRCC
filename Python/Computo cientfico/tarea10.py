# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 21:04:46 2017

@author: Ricardo Chávez Cáliz
"""
import numpy as np
from numpy import zeros
from numpy import log
from numpy.random import rand
from scipy.stats import norm
from numpy.random import normal
from time import clock
import matplotlib.pyplot as plt

def evalua(M):
    n=len(M)
    x = np.arange(0., n, 1)
    y = log(M)
    
    plt.figure(figsize=(5, 5), dpi=80)
    plt.plot(x, y, color='#990033')
    plt.xlabel('numero de iteracion (n)')
    plt.ylabel('log(media simulada)')
    plt.title('Burn-in')
    plt.savefig('BurnNC.png')
    plt.show()
    
def propuesta(numP,m,z1,z2,z3,z4,z5,datos):
    """
    Devuelve distintas propuestas de acuerdo a la opción selecccionada.
    Input: int (opción deseada)
    Output: float, float, float (pp,Np,ro)
    """
    a=21
    q=norm.cdf(a-m)
    u=rand(1)*(1-q) + q
    
    if numP>6 or numP<1:
        raise ValueError("No conozco esa propuesta") 

    #Condicionales totales de Gibs
    
    #Mueve a teta
    if numP==1:
        todos = np.concatenate((datos,np.array([z1,z2,z3,z4,z5])),axis=0)
        mp=np.mean(todos)
        return mp,z1,z2,z3,z4,z5   
 
    #Mueve a z1
    elif numP==2:
        z1p=norm.ppf(u)+m
        return m,z1p,z2,z3,z4,z5
    
    #Mueve a z2
    elif numP==3:
        z2p=norm.ppf(u)+m
        return m,z1,z2p,z3,z4,z5
    
   #Mueve a z3
    elif numP==4:
        z3p=norm.ppf(u)+m
        return m,z1,z2,z3p,z4,z5
 
   #Mueve a z4
    elif numP==5:
        z4p=norm.ppf(u)+m
        return m,z1,z2,z3,z4p,z5
 
   #Mueve a z5
    else:
        z5p=norm.ppf(u)+m
        return m,z1,z2,z3,z4,z5p

def GibsNormal(datos,tam):
    """
    Completa la muestra con un Gibs sampler con parámetros artificiales en
    la posterior  f(m, z_{1}, .... z_{n-m}|x) y s implementa el algoritmo MCMC
    
    Input:
    Output: array(tam) (muestra)
    """
    M = zeros((tam,6))#m,z1,z2,z3,z4,z5
    datos=datos[0:15]    
    a=21
     
    #Inical.
    m=normal(18,1)
    q=norm.cdf(a-m)
    
    M[0,0] = m #m inicial
    for i in xrange(1,6):
        u=rand(1)*(1-q) + q
        M[0,i] = norm.ppf(u)+m
        
    for i in xrange(1,tam):
        m=M[i-1][0]
        z1=M[i-1][1]
        z2=M[i-1][2]
        z3=M[i-1][3]
        z4=M[i-1][4]
        z5=M[i-1][5]
        
        numP = int(6*rand(1)[0])+1
        
        R = propuesta(numP,m,z1,z2,z3,z4,z5,datos)
        for j in xrange(0,6):
            M[i,j] = R[j]

    return M
#------------------------------------------------------------------------------
def dx(f, x):
    return abs(0-f(x))
 
def newton(f, df, x0, e):
    delta = dx(f, x0)
    i=0
    while delta > e:
        i=i+1
        x0 = x0 - f(x0)/df(x0)
        delta = dx(f, x0)
    print 'Raíz es: ', x0
    print 'se llegó a la solución en ', i , ' iteraciones'
    	
def f(t):
    """
    derivada de log de verosimilitud 
    """
    X = np.array([ 18.221753,  18.418174,  18.720224,  19.067637,  19.128777,  19.402623,
           19.507782,  19.580571, 19.632340,  19.930952,  20.116566,  20.142095,  20.445327,
           20.461254,  20.646856])
    
    m=15.0
    n=20.0
    
    return sum(X-t) + (n-m)*(norm.pdf(a-t))/(1-norm.cdf(a-t))
 
def df(t):
    """
    segunda derivada de log de verosimiitud
    """
    n = 20.0
    m = 15.0
    a = 21.0
    return -m - ((n-m)/(1-norm.cdf(a-t))**2)*(norm.pdf(a-t)*(a-t)*(norm.cdf(a-t)-1) + (norm.pdf(a-t))**2)
#------------------------------------------------------------------------------
def algoritmoEM(r,m,a):    
    n=len(r)
    r=np.sort(r)
    B=r[0:m]
    X=np.mean(B)
    teta = np.mean(B)
    
    m=float(m)
    n=float(n)
    a=float(a)
    
    for i in xrange(100):
        tetaCh= (m/n)*X + ((n-m)/n)*(teta + ( (norm.pdf(a-teta)) / (1-norm.cdf(a-teta)) ) )
        if np.abs(tetaCh - teta) < 1e-6:
            break
        teta=tetaCh
    print 'en ', i , ' iteraciones'
    return tetaCh

#------------------------------------------------------------------------------
if __name__ == "__main__":
    
    """
    1.EM: Considere que se tienen los datos de distribución Normal truncada.
    Conocemos el logaritmo de la distribución conjuta y la función de 
    verosimilitud. Implementar el algoritmo EM para obtener el estimador máximo
    verosímil de θ.
    """
    np.random.seed(0)
    n=1000
    a=1
    r = norm.rvs(size=n)
    c=0
    for i in xrange(0,len(r)):
        if r[i]>a:
            c=c+1
            r[i]=a
            
    b=np.sort(r)
    b=b[0:n-c]
    s, bins, patches = plt.hist(b, 20, normed=1, facecolor='#006666', alpha=0.75,edgecolor='#003366', linewidth=1.2)
    plt.title('Histograma de prueba de control')
    plt.savefig("EM-Prueba")
    plt.show()    
    
    print "Estimador de la media ignorando censura es " , np.mean(b)
    print "Estimador de la media con algoritmoEM ", algoritmoEM(r,n-c,a)

    #--------------------------------------------------------------------------
    datos = np.array([ 18.221753,  18.418174,  18.720224,  19.067637,  19.128777,  19.402623,
           19.507782,  19.580571, 19.632340,  19.930952,  20.116566,  20.142095,  20.445327,
           20.461254,  20.646856, 21.000000, 21.000000,  21.000000,  21.000000,  21.000000 ])
    
    n, bins, patches = plt.hist(datos[0:15], 5, normed=1, facecolor='#006666', alpha=0.75,edgecolor='#003366', linewidth=1.2)
    plt.title('Histograma de datos censurados')
    plt.savefig("EM-Datos")
    plt.show()

    m = 15
    a = 21
    print "Estimador de la media ignorando censura es " , np.mean(datos)
    tiempo_inicial = clock()
    print "Estimador de la media con algoritmoEM ", algoritmoEM(datos,m,a)
    print "Tiempo de ejecución ", (clock() - tiempo_inicial)
    print "-------------------------------------------------------------------"   
    """"
    2. MCMC: Bajo el mismo esquema de datos del ejercicio anterior, poner
    la distribuci´on a priori para θ ∼ N(18, 1) e implementar un MCMC
    cuya función objetivo es f(θ, z_m+1, . . . , z_n|x).
    
    Es decir, se agregan parámetros espurios z para salvar el problema de
    la censura; la misma idea de EM. Utilice solamente kérneles Gibbs.
    Para simular de las condicionales totales f(zi|θ, x) utilice el método de
    la Transformada Inversa.
    """
    N=5000
    tiempo_inicial = clock()
    M = GibsNormal(datos,N)
    print "Tiempo de ejecución ", (clock() - tiempo_inicial) , " para N=", N
    
    print "Promedio de las medias simuladas " , np.mean(M[0:,0])

    evalua(M[0:,0])
    
    k=20
    y=M[k:,0]
    x=M[k:,1]

    #Scatter    
    colors = np.arange(0, 1, 1.0/(N-k))
    area = 50*np.ones(N)
    plt.figure(figsize=(5, 5), dpi=80)
    plt.scatter(x, y, s=area, c=colors, alpha=0.5) 
    plt.xlabel('Valores de media')
    plt.ylabel('valores de z1 artificial')
    plt.savefig('ScatterNC.png')
    plt.show()

    #Trayectoria
    N=300
    M = GibsNormal(datos,N)

    k=10
    y=M[k:,0]
    x=M[k:,1]

    area = 50*np.ones(N)
    colors = np.arange(0, 1, 1.0/(N-k))
    plt.figure(figsize=(5, 5), dpi=80)
    plt.scatter(x, y, s=area, c=colors, alpha=0.8)
    plt.plot(x, y, color='#993366')
    plt.savefig('TrayectoriaNC.png')
    plt.xlabel('Valores de media')
    plt.ylabel('valores de z1 artificial')
    plt.title('Trayectoria de la cadena')
    plt.show()
    
    """
    3. Método de Newton-Raphson: Bajo el mismo esquema de EM, implementar
    el m´etodo de Newton-Raphson para encontrar a ˆθ que maximice
    la log-verosimilitud log L(x|θ) (MLE)
    """
    print '-------------------------------------------------'
    print("Aproximación por Newton-Rapson")
    tiempo_inicial = clock()
    x0 = np.mean(datos)
    newton(f, df, 15, 1e-6)
    print "Tiempo de ejecución", (clock() - tiempo_inicial)