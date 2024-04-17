import numpy as np

'''
    Config file to set Simulation Parameters
    u_n[k,j] -> k = y_index, j = x_index
'''

# X Domain
xmin, xmax = 0, 2*np.pi     # Domain limits
Lx = abs(xmax - xmin)       # Domain Length
x_ghost_cells = 0           # Number of ghost cells
Nx = 100 + x_ghost_cells    # Number of spatial points
dx = Lx / (Nx-1)    # Cell width
# dx = 0.01 # Cell width
# Nx = int(L/dx) + 1    # Number of spatial points
x = np.linspace(xmin + dx/2, xmax - dx/2, Nx)  # Spatial grid

# Y Domain
ymin, ymax = 0, 1
Ly = abs(xmax - xmin)
y_ghost_cells = 0
Ny = 100 + y_ghost_cells
dy = Ly / (Ny-1)
# dy = 0.01
# Nx = int(Ly/dy) + 1
y = np.linspace(ymin + dy/2, ymax - dy/2, Nx)

X, Y = np.meshgrid(x, y)

# First Row, Col
x0 = X[0, :]
y0 = Y[:, 0]

# x Indices
j = np.arange(0, Nx)  # Indices 'j'
j_plus1 = np.roll(j, -1)  # Indices 'j+1'
j_min1 = np.roll(j, 1)   # Indices 'j-1'

midpoint = y[Ny//2]

cfl = 0.35  # Stability Parameter - CFL Number

Tf = 0.5  # Final time / Total Time

x_first_cell, x_last_cell = 0, Nx-1  # j domain Limits
y_first_cell, y_last_cell = 0, Ny-1  # k domain Limits


K = 1 # Coupling Strength


# Returns Non-local flux (1D array), f(u) = L[Mu]Mu
def flux(L_U, u, xindices, y_index,):
   return -K*L_U[y_index, xindices] + y0[y_index]*u[y_index, xindices]
