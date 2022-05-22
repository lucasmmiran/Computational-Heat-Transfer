# Autor: Lucas Mendes Miranda
# Data: 20/05/2022

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.cm as cm

# parametros da simulacao
Q = 0.3
rho = 1.0
cv = 1.0
k = 1.0
alpha = k/(rho*cv)
dt = 0.005

# parametros da malha
Lx = 1.0
Ly = 1.0
nx = 8
ny = 8
npoints = nx*ny
ne = (nx-1)*(ny-1)
dx = Lx/(nx-1)
dy = Ly/(ny-1)

# Formando os pontos da malha em X e Y
Xv = np.linspace(0, Lx, nx )
Yv = np.linspace(0, Ly, ny )

# Criando o Grid
x, y =np.meshgrid(Xv, Yv)

#Colocando no formato desejado (Transformando a matriz para um vetor)
X = np.reshape(x, npoints)
Y = np.reshape(y, npoints)

# Criando os vetores de CC
cc1, cc2, cc3, cc4, miolo  = [], [], [], [], []

# Preenchendo os vetores 1 e 3
for i in range(0, nx):
    cc1 = cc1 + [i]
    cc3 = cc3 + [(npoints - i - 1)]

# Preenchendo os vetores 2 e 4
for i in range(0, ny-2):
    cc2 = cc2 + [((i + 2) * nx) -1]
    cc4 = cc4 + [(i + 1) * nx]

cc = cc1 + cc2 + cc3 + cc4

for i in range(0, npoints-1):
    if not i in cc:
        miolo = miolo + [i]

#Vetor de condições de contorno
Tcc = np.zeros( (npoints), dtype = 'float')

# Preenchendo o vetor de CCs
#CC1
for i in cc1:
    Tcc[i] = X[i]

# CC2
for i in cc2:
    Tcc[i] = Y[i]*Y[i] + 1.0

# CC3
for i in cc3:
    Tcc[i] = X[i]*X[i] + 1.0

# CC4
for i in cc4:
    Tcc[i] = Y[i]

# Vetor de condicao inicial
T = np.zeros( (npoints),dtype='float' )

# Aplicando as c.c.s na condicao inicial
for i in cc:
    T[i] = Tcc[i]

# Vetor fonte de calor
Qvec = Q*np.ones( (npoints),dtype='float' )
Z = T.reshape(ny,nx)

# plot 3D da condicao incial + condicao de contorno
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, Z, cmap=cm.turbo)
ax.set_xlabel('Eixo X')
ax.set_ylabel('Eixo Y')
ax.set_zlabel('Temperatura')

plt.show()

for n in range(0,10):
    for i in miolo:
        T[i] = T[i] + \
            dt*alpha*(T[i+1]-2*T[i]+T[i-1])/(dx*dx) + \
            dt*alpha*(T[i+nx]-2*T[i]+T[i-nx])/(dy*dy) + \
        Qvec[i]/(rho*cv) 

    Z = T.reshape(ny,nx)

    # Plot 3D 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, Z, cmap=cm.turbo)
    ax.set_xlabel('Eixo X')
    ax.set_ylabel('Eixo Y')
    ax.set_zlabel('Temperatura')

    plt.show()