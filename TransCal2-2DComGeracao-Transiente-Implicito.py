# Autor: Lucas Mendes Miranda
# Data: 22/05/2022

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.cm as cm

# Parâmetros da simulacao
Q = 10
rho = 1.0
cv = 1.0
k = 1.0
alpha = k/(rho*cv)
dt = 1.0

# Parâmetros da malha
Lx = 1.0
Ly = 1.0
nx = 4
ny = 4
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

# Vetor de condições de contorno da Temperatura
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
Qvec = (Q/(rho*cv))*np.ones( (npoints),dtype='float' )


# Inicializando matriz e vetor A e b
A = np.zeros( (npoints,npoints),dtype='float' )

for i in cc:
    A[i,i] = 1.0

for i in miolo:
    A[i,i+1] = - alpha/(dx*dx)
    A[i,i] = 1/dt +2*alpha/(dx*dx) + 2*alpha/(dy*dy)
    A[i,i-1] = - alpha/(dx*dx)
    A[i,i+nx] = - alpha/(dy*dy)
    A[i,i-nx] = - alpha/(dy*dy)

# Calculando a inversa de A
Ainv = np.linalg.inv(A)

# # plot 3D da condicao incial + condicao de contorno
# Z = T.reshape(ny,nx)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(x, y, Z, cmap=cm.turbo)
# ax.set_xlabel('Eixo X')
# ax.set_ylabel('Eixo Y')
# ax.set_zlabel('Temperatura')

# plt.show()

b = np.zeros( (npoints),dtype='float' )
# contorno
for i in cc:
    b[i] = Tcc[i]

for n in range(0,10):
    for i in miolo:
        b[i] = T[i]/dt + Qvec[i]

for n in range(0,10):
    for i in miolo:
        b[i] = T[i]/dt + Qvec[i]
    
    T = Ainv@b

    Z = T.reshape(ny,nx)

    # Plot 3D 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, Z, cmap=cm.jet)
    ax.set_xlabel('Eixo X')
    ax.set_ylabel('Eixo Y')
    ax.set_zlabel('Temperatura')

    plt.show()