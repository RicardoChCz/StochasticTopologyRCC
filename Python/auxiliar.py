# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 15:03:11 2018

@author: rechavez
"""
from numpy.random import randint
from numpy import arange

def randomSet(k,n):
    """
    Returns a random set of size k from the set [0,1,...,n]
    Input: int, int
    Output: List
    """
    if k > n:
        raise ValueError("You must pick k smaller than n")
    S= set()
    j=0
    while j<k:
        S.add(randint(n))
        j = len(S)
    return S
    
def powerset(seq):
    """
    Returns all the possible subsets that can be made with the input set.
    Input: Set
    Output: list
    """
    seq = list(seq)
    
    #Empty set or one element sets
    if len(seq) <= 1:
        yield seq
        yield []
        
    else:
        for item in powerset(seq[1:]):
            yield [seq[0]]+item
            yield item

def powersetEfective(seq,kMin,kMax):
    """
    Returns all the possible subsets that can be made with the input set.
    Input: Set
    Output: list
    """
    
    seq = list(seq)
    
    #Empty set or one element sets
    if len(seq) <= 1:
        yield seq
        yield []

    else:
        for item in powerset(seq[1:]):
            if (len([seq[0]]+item) <= kMax and len([seq[0]]+item) >= kMin):
                yield [seq[0]]+item
            if (len(item) <= kMax and len(item) >= kMin):    
                yield item
            
def intervalo (x0, x1, d):
    """
    Devuelve un intervalo discretizado que empieza en x0 y termina en x1
    con separación d entre cada punto de I.
    Input:  3 números: extremo izquierdo x0, extremo derecho x1,
            separación d.
    Output: Un arreglo con valores que empiezan en x0, terminan en x1 y
            difieren en d
    """
    I = arange(x0, x1, d)
    return I
    
def is_notempty(any_structure):
    if any_structure:
        return True
    else:
        return False
    