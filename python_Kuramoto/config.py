import numpy as np

'''
    Config file to set Simulation Parameters
'''

equation = 2  # Problem Statement

xmin, xmax = -2 * np.pi, 2 * np.pi  # Domain limits
L = abs(xmax - xmin)  # Domain Length
ghost_cells = 2  # Number of ghost cells
Nx = 1000 + ghost_cells  # Number of spatial points
dx = L / (Nx - 1)  # Cell width
# dx = 0.01 # Cell width
# Nx = L/dx + ghost_cells # Number of spatial points

cfl = 0.75  # Stability Parameter - CFL Number
c = 1.0  # Wave Velocity

dt = cfl * dx / c  # Time step
Tf = 2.0  # Final time / Total Time
Nt = int(Tf / dt)  # No. of time steps

first_cell, last_cell = 1, Nx - 2  # j domain Limits


# Returns flux, f(u)
a = 1.0  # Constant Flux Multiplier
def flux(u):
    if equation == 1:
        return a * u
    elif equation == 2:
        return 0.5 * u * u
