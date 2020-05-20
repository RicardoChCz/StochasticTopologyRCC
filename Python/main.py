#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 23:30:28 2019

@author: ricardochcz
"""
import networkx as nx
from _aux import random_set

#Computational experiments:
from _rigidity_experiments import individual_experiments
from _rigidity_experiments import set_of_experiments

#Simulation of Rigid expansion as a Markov Chain
from _rigid_expansion_markov_chain

#Erdös-Renyi thresholds
from _thresholds_erdos_renyi

#Random Graphs Simulations
from _random_graphs

#Rigid Expansion Visualization
from _aux_graphic

#Optimization Analisis
from 

if __name__ == "__main__":
    """
    -------------------------------------------------------------------------
      Computational experiments related to rigid expansions.
    -------------------------------------------------------------------------    
      Event 1: Fixed vertex
      Event 2: Any vertex
      Event 3: Expansion
    -------------------------------------------------------------------------
    """
    n,p=15,0.5
    print("EMPIRICAL VS ESTIMATED PROBABILITIES. n="+ str(n)+ " , p=" + str(p))
    print("Event 1. Prob. of a fixed vertex being uniquely det.")
    individual_experiments(15,0.5,event=1)

    print("Event 2. Prob. of any vertex being uniquely det.")
    individual_experiments(15,0.5,event=2)
    
    print("Event 3. Prob. that a subset generates a rigid expansion.")
    individual_experiments(15,0.5,event=3)

    #------------------------------------------------------------------------
    N,P=[5,7,10],[0.2,0.5,0.8]
    print("EMPIRICAL VS ESTIMATED PROBABILITIES. N="+ str(n)+ " , P=" + str(p))
    print("Event 1. Prob. of a fixed vertex being uniquely det.")
    set_of_experiments(N,P,1)

    print("Event 2. Prob. of any vertex being uniquely det.")
    set_of_experiments(N,P,2) 

    print("Event 3. Prob. that a subset generates a rigid expansion.")
    set_of_experiments(N,P,3)
    
    """
    -------------------------------------------------------------------------
      Simulation of Rigid expansion as a Markov Chain
    -------------------------------------------------------------------------    
    """
    
    """
    -------------------------------------------------------------------------
      Erdös-Renyi thresholds
    -------------------------------------------------------------------------    
    """

    """
    -------------------------------------------------------------------------
      Random Graphs Simulations
    -------------------------------------------------------------------------    
    """

    """
    -------------------------------------------------------------------------
      Rigid Expansion Visualization
    -------------------------------------------------------------------------    
    """

    """
    -------------------------------------------------------------------------
      Optimization Analysis
    -------------------------------------------------------------------------    
    """

    """
    -------------------------------------------------------------------------
      Generate algorithm
    -------------------------------------------------------------------------    
    """