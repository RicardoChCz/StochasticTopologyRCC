# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 11:40:57 2017

@author: Ricardo Esteban Chávez Cáliz
"""

import numpy as np
import random
import networkx as nx
import matplotlib.pyplot as plt
from time import clock
from numpy import log
from numpy.random import rand
from numpy.random import randint
from __aux_graphic import color_template_6
from __aux_graphic import plot_lines_styles


def _bernoulli(p):
    """
    Simulates a Bernoulli(p) r.v
    Input: float
    Output: int (success or failure)
    """
    M = rand(1)
    if p > 1 or p<0:
        raise ValueError("Wrong paramater value")
    if M < p:
        return True
    else:
        return False

def _empiric_connectedness(n,p,M):
    """
    Calculates empiric probabilty of obtainign a coneccted random graph in G(n,p)
    Input: int n, float 0<p<1, int M - number of experiments
    Output: float (empiric probability)
    """
    r=0.0
    for i in range(0,M):
        if nx.is_connected(er_graph(n,p))==True:
            r=r+1
    return r/M

def er_graph(n,p):
    """
    Return a Erdös- Rényi random graph with parameters n and p
    Input: int n, float p en [0,1]
    Output: object (networkX)
    """
    G=nx.Graph()
    G.add_nodes_from(list(range(0, n)))
    for i in range(0,n):
        for j in range(i+1,n):
            if _bernoulli(p):
                G.add_edge(i,j)
    return G


def generate_algorithm(G):
    """
    Sample a random tree using Generate algorithm
    Input: object (networkX graph)
    Output: object (networkX tree)
    """
    V=G.nodes()
    n=len(V)
    #punto inicial
    u = randint(n)
    visited = [u]
    T=nx.empty_graph(n)
    
    while len(visited) < n:
        N = G.neighbors(u)
        v = random.choice(list(N))
        if (v in visited)==False:
            visited.append(v)
            T.add_edge(u,v)
        #change standing-point
        u=v
    
    return T
    
def draw_er_graphs(n):
    """
    Graph a set of Erdös-Rényi random graphs (varing p)
    Input: int n
    Output: .png file
    """
    r = c =3
    k=1.0
    f, axarr = plt.subplots(r, c,figsize=(6, 6.3), dpi=80)
    for i in range(0,3):
        for j in range(0,3):
            G = er_graph(n,(k/10))
            nx.draw_shell(G,node_size=30, node_color='#0d3c55',edge_color='#595959', ax=axarr[i, j])
            axarr[i, j].set_title('p = ' + str(int(k)) + "/" + str(n))
            k=k+1
    plt.savefig('Figures/Set-of-ER-graphs.png')
    f.subplots_adjust(hspace=0.4)
    plt.show()

def draw_generate_trees(n):
    """
    Graph a sample of random trees with Generate algorithm
    Input: int n
    Output: .png file
    """
    r = c = 4
    f, axarr = plt.subplots(r, c,figsize=(6, 6), dpi=80)
    for i in range(0,4):
        for j in range(0,4):
            T = generate_algorithm(nx.complete_graph(n))
            nx.draw_spectral(T,node_size=20, node_color='#0d3c55', edge_color='#595959', ax=axarr[i, j])
    plt.savefig('Figures/Tree-sample'+str(n)+'.png')
    plt.show()

def time_execution_performance(experimental_interval, number_of_experiments, xlabel,algorithms,label):
    plt.gcf().clear()
    fig, axarr = plt.subplots(1, 2, figsize=(6,3),dpi=80)    
    colors = color_template_6()
    styles = plot_lines_styles()
    c=0
    for algorithm in algorithms:
        y1 = np.zeros(len(experimental_interval))
        p = algorithm["p"]
        algorithm_name = algorithm["name"]
        j=0
        for n in experimental_interval:
            sumation = 0
            for e in range(number_of_experiments):
                starting_time = clock()
                if(algorithm_name == "er_graph"):
                    er_graph(n,p)
                elif(algorithm_name == "generate_algorithm"):
                    generate_algorithm(nx.complete_graph(n))
                elif(algorithm_name == "er_fast"):
                    m = np.random.binomial(n*(n-1)/2, p)
                    nx.gnm_random_graph(n, m)
                sumation += clock() - starting_time
                
            y1[j] = sumation/number_of_experiments
            j += 1
            
        axarr[0].plot(experimental_interval, y1, 
                     linestyle=styles[c%len(algorithms)], color=colors[c])
        axarr[1].plot(experimental_interval, log(y1), 
                     linestyle=styles[c%len(algorithms)], color=colors[c], 
                     label=algorithm_name+", p="+str(p))
        
        axarr[0].set_title('Linear scale')
        axarr[1].set_title('Logarithmic scale')
        axarr[0].set_ylabel('Time')
        axarr[0].set_xlabel(xlabel)
        axarr[1].set_xlabel(xlabel)
        c += 1
        
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.savefig("Figures/Time-execution-"+label+".png", bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    """
    Sample Erdös-Rényi graphs
    """
    print ("Sample Erdös-Rényi graphs")
    draw_er_graphs(10)
    """
    Sample random trees with Generate algorithm
    """
    print ("Sample random trees with Generate algorithm")
    draw_generate_trees(10)
    """
    Time execution for ER graphs
    """
    print ("Time execution test for ER-graphs samples. Might take a while...")

    test_1 = {"name": "er_fast" ,"p": 0.1}
    test_2 = {"name": "er_fast" ,"p": 0.9}
    test_3 = {"name": "er_graph","p": 0.1}
    test_4 = {"name": "er_graph","p": 0.9}
    
    upper_bundle = 500
    experimental_interval = np.arange(5, upper_bundle+55, 50)
    
    time_execution_performance(experimental_interval, 5, xlabel="Number of vertices", 
                               algorithms=[test_1, test_2, test_3, test_4], 
                               label="er-generation-algoriths")

    """
    Time execution for RGT graphs
    """
    print ("Time execution test for RGT. Might take a while...")

    test_1 = {"name": "generate_algorithm","p": 0.1}
    test_2 = {"name": "generate_algorithm","p": 0.9}
    
    upper_bundle = 500
    experimental_interval = np.arange(5, upper_bundle+55, 50)
    
    time_execution_performance(experimental_interval, 30, xlabel="Number of vertices", 
                               algorithms=[test_1, test_2], 
                               label="generate-algorithm")