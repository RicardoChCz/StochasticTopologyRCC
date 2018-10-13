# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 05:14:51 2017

@author: Ricardo Chávez Cáliz
"""
import numpy as np
from matplotlib.mlab import bivariate_normal
from matplotlib import pylab as plt

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
