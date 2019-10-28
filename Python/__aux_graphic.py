#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 13:21:40 2019
@author: Ricardo Chávez Cáliz - RCC

Methods to generate graphical visualization of rigid expansions
"""
import imageio
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from numpy.random import randint

from __aux import random_set

def layout_graph(G, S):
    #Position of vertices
    pos=nx.circular_layout(G)
    #Size
    fig = plt.figure(figsize=(7, 7))
    #Background
    fig.set_facecolor('white')
    #Colors
    colors=['#020A18']*len(G.nodes())
    #Paint S
    for j in S:    
        colors[j]='#FFA500'

    #Edges
    nx.draw_networkx_edges(G, pos, alpha=0.8, edge_color='#020a18')
    #Vertices
    nx.draw_networkx_nodes(G, pos,
                           node_size=200, #Tamaño de nodos
                           node_color=colors) #Define los colores
    nx.draw_networkx_labels(G, pos)

    plt.xlim(-1.05, 1.2)
    plt.ylim(-1.05, 1.05)
    plt.axis('off')

def visual_rigid_exp(G,S,filenames,i):
    layout_graph(G,S)
    plt.savefig('Figures/rigid_expansion_gif/rigid_exp_'+str(i)+'.png')
    filenames.append('Figures/rigid_expansion_gif/rigid_exp_'+str(i)+'.png')
    plt.show()

def color_template_3():
    return(["#0d3c55","#f16c20","#1395ba"])

def color_template_6():
    return(["#0d3c55","#1395ba","#a2b86c","#ebc844","#f16c20","#c02e1d"])

def color_template_12():
    return(["#0d3c55","#0f5b78","#117899","#1395ba","#5ca793","#a2b86c","#ebc844","#ecaa38","#ef8b2c","#f16c20","#d94e1f","#c02e1d"])

def plot_lines_styles():
    return(["-","--","-.",":"])