# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 16:39:13 2017

@author: Pc
"""

from numpy.random import uniform
from numpy import sin
from numpy import exp
from numpy import sqrt


if __name__ == "__main__":
    
    #Número de repeticiones
    n=1000
    
    #Muestra
    X = uniform(1,2,n)    
    #Función
    G = (X)**3/(1+ sqrt(X))    
    #Calcula la suma
    suma = 0.0
    for g in G:
        suma = suma + g        
    #Calcula valor esperado empírico
    I = (1.0/n)*suma    
    #Calcula la integral
    print I
    
    #Muestra
    X = uniform(0,3,n)    
    #Función
    G = (exp(X)*sin(X))/(1+ X**2)    
    #Calcula la suma
    suma = 0.0
    for g in G:
        suma = suma + g        
    #Calcula valor esperado empírico
    E = (1.0/n)*suma    
    #Calcula la integral
    I = 3.0*E
    print I