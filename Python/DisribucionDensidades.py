# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 02:38:37 2017

@author: Ricardo Chávez Cáliz
"""
#Paquetería necesaria
import numpy as np
from numpy.random import normal
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV

#Datos dados
datos = normal(0,1,150)

#Tamaño de muestra
N=150
#Ajuste de formato de datos
X= np.asarray(datos).reshape(-1,1)

#Gráficar
X_plot = np.linspace(-4, 4, 1000)[:, np.newaxis]
bins = np.linspace(-4, 4, 40)
fig, ax = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(8, 6), dpi=80)
fig.subplots_adjust(hspace=0.05, wspace=0.05)
# Histograma 1
ax[0, 0].hist(X[:, 0], bins=bins, fc='#006666', normed=True, alpha=0.75,edgecolor='#003366', linewidth=1.2)
# Histograma 2
ax[0, 1].hist(X[:, 0], np.linspace(-4, 4, 8) , fc='#006666', normed=True, alpha=0.75,edgecolor='#003366', linewidth=1.2)
# Tophat KDE
kde = KernelDensity(kernel='tophat', bandwidth=0.88).fit(X)
log_dens = kde.score_samples(X_plot)
ax[1, 0].fill(X_plot[:, 0], np.exp(log_dens), fc='#660033', alpha=0.5,edgecolor='#4d0026', linewidth=1.2)
# Gaussiana KDE
kde = KernelDensity(kernel='gaussian', bandwidth=0.88).fit(X)
log_dens = kde.score_samples(X_plot)
ax[1, 1].fill(X_plot[:, 0], np.exp(log_dens), fc='#660033', alpha=0.5,edgecolor='#4d0026', linewidth=1.2)

#Ajustar limites de graficas
for axi in ax.ravel():
    axi.plot(X[:, 0], np.zeros(X.shape[0]) - 0.01, '+k')
    axi.set_xlim(-4, 4)
    axi.set_ylim(-0.02, 0.5)
#Agregar etiquetas
for axi in ax[:, 0]:
    axi.set_ylabel('Densidad normalizada')
for axi in ax[1, :]:
    axi.set_xlabel('x')
plt.savefig("histoYkernel")

#------------------------------------------------------------------------------
# Grafica con distintos kerneles 1D
X_plot = np.linspace(-4, 4, 1000)[:, np.newaxis]
fig, ax = plt.subplots(figsize=(8, 6), dpi=80)

#Elegir los kerneles a usar
for kernel in ['gaussian', 'tophat', 'epanechnikov']:
    kde = KernelDensity(kernel=kernel, bandwidth=4).fit(X)
    log_dens = kde.score_samples(X_plot)
    ax.plot(X_plot[:, 0], np.exp(log_dens), '-',
            label="kernel = '{0}'".format(kernel))

#Caja con leyenda
ax.legend(loc='upper left')
ax.plot(X[:, 0], -0.005 - 0.01 * np.random.random(X.shape[0]), '+k')

#Límites
ax.set_xlim(-4, 4)
ax.set_ylim(-0.02, 0.5)
plt.savefig("bandaChafa")
plt.show()

#------------------------------------------------------------------------------
#Calcula la amplitud, revisa entre -1 y 1
params = {'bandwidth': np.logspace(-1, 1, 20)}
grid = GridSearchCV(KernelDensity(), params)
grid.fit(X)
print("Mejor amplitud de banda: {0}".format(grid.best_estimator_.bandwidth))

#Segundo intento con ventana apropiada
X_plot = np.linspace(-4, 4, 1000)[:, np.newaxis]
fig, ax = plt.subplots(figsize=(8, 6), dpi=80)

for kernel in ['gaussian', 'tophat', 'epanechnikov']:
    kde = KernelDensity(kernel=kernel, bandwidth=0.88).fit(X)
    log_dens = kde.score_samples(X_plot)
    ax.plot(X_plot[:, 0], np.exp(log_dens), '-',
            label="kernel = '{0}'".format(kernel))

ax.legend(loc='upper left')
ax.plot(X[:, 0], -0.005 - 0.01 * np.random.random(X.shape[0]), '+k')

ax.set_xlim(-4, 4)
ax.set_ylim(-0.02, 0.5)
plt.savefig("bandaChida")
plt.show()
