import numpy as np

'''
    Config file to set Simulation Parameters
'''

# X Domain
xmin, xmax = 0, 2*np.pi     # Domain limits
Lx = abs(xmax - xmin)       # Domain Length
x_ghost_cells = 0           # Number of ghost cells
Nx = 512 + x_ghost_cells    # Number of spatial points
dx = Lx / (Nx-1)    # Cell width
# dx = 0.01 # Cell width
# Nx = int(L/dx) + 1    # Number of spatial points
x = np.linspace(xmin + dx/2, xmax - dx/2, Nx+1)  # Spatial grid


# Y Domain
ymin, ymax = 0, 1
Ly = abs(xmax - xmin)
y_ghost_cells = 0
Ny = 512 + y_ghost_cells
dy = Ly / (Ny-1)
# dy = 0.01
# Nx = int(Ly/dy) + 1
y = np.linspace(ymin + dy/2, ymax - dy/2, Nx+1)


X, Y = np.meshgrid(x, y)

midpoint = y[Ny//2]


cfl = 0.35  # Stability Parameter - CFL Number

Tf = 1.0  # Final time / Total Time

first_cell, last_cell = 0, Nx-1  # j domain Limits


om = None # Natural Frequency
K = 1 # Coupling Strength


# Returns 
def flux(V_U, u, x_index):
    pass

