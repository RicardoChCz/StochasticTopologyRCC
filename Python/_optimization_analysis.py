# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:55:08 2018
@author: Ricardo Chávez Cáliz - RCC

Optimization analysis for rigid expansions
"""
from _execution_tests import er_graph
from __aux import random_set
from __aux import intervalo
from __aux_graphic import color_template_3
from __aux_graphic import plot_lines_styles
from __rigid_expansion import rigid_expansion

from numpy import log
from time import clock
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from math import floor

def calculate_mean_expansion_time(n, p, k, case, number_of_experiments):
    """
    Input: n (int), p (float), k(int), case(dict), number_of_experiments(int)
    Output: Mean time, float
    """
    sumation = 0
    for e in range(number_of_experiments):
        A = random_set(k, n)
        starting_time = clock()
        if(case['optimizations'] == "none"):
            G = er_graph(n,p)
            rigid_expansion(G,A,p,False,False)
        else:
            G = nx.fast_gnp_random_graph(n, p)
            rigid_expansion(G,A,p,False,True)
        sumation += clock() - starting_time
    return 1000 * sumation/number_of_experiments

def time_execution_performance(n, p, samples_size, log_scale=False):
    """
    Plots a mean time rigid expansions experiments varing n, p  and k.
    Input: n (max graph size), p (list), samples_size(int)
    Output: Mean time, float
    """

    n_interval = intervalo(3,n+2,2)
    k_proportion = [0.25, 0.5, 0.75]
    # For the plotting
    rows = len(k_proportion)
    cols = len(p)
    fig, axarr = plt.subplots(rows, cols, figsize=(3*cols+3,3*rows),dpi=80)
    colors = color_template_3()
    styles = plot_lines_styles()

    cases = [
      {"name":"Whitout optimizations", "optimizations":'none'},
      {"name":"Fully optimized", "optimizations":"all"}]

    for r in range(len(k_proportion)):
        _proportion = k_proportion[r]
        for c in range(len(p)):
            _p = p[c]
            for case in cases:
                case_index = cases.index(case)
                t = np.zeros(len(n_interval))
                for i in range(len(n_interval)):
                    _n = n_interval[i]
                    _k = floor(_n * _proportion)
                    t[i] = calculate_mean_expansion_time(_n, 
                      _p, _k, case, sample_size)
                if log_scale:
                    t = log(t)
                axarr[r,c].plot(n_interval, t, label=case['name'],
                    linestyle=styles[case_index%len(cases)], color=colors[case_index])
            # Set labels for plot in row r and col c
            print('Just finished', r,c)
            axarr[r, c].spines['top'].set_color('white')
            axarr[r, c].spines['right'].set_color('white')
            axarr[r, c].set_title('p = ' + str(_p) + ', k_prop = ' + str(_proportion), backgroundcolor='silver')
            axarr[r, c].set_ylabel('Time (ms)' + (' (log scale)' if log_scale else ''))
            axarr[r, c].set_xlabel('Number of vertices')
            if r == 0 and c == len(p)-1:
                axarr[r,c].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    fig.subplots_adjust(bottom=0.2)
    plt.tight_layout()
    plt.savefig('Figures/Time-execution-rigid-expansions'+('-log-scale' if log_scale else '')+'.png', bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    """
    Time execution for rigid expansions
    """
    print ("Time execution test for Rigid expansions. Might take a while...")
    sample_size = 30
    n = 20
    p = [0.1, 0.5, 0.9]

    time_execution_performance(n, p, sample_size)
    time_execution_performance(n, p, sample_size, True)