from scipy.fft import fft, ifft
from .functions import *

'''
    Finite Volume Solver for Kuramoto Equation with NON-IDENTICAL natural fruequencies  YET TO BE IMPLEMENTED 
'''

# Returns Lax-Friedrich Numerical Flux 
def num_flux(V_U, u, x_index, dt):
    if x_index == Nx-1: x_index = 0 # Right Boundary
    return 0.5 * (flux(V_U, u, x_index) + flux(V_U, u, x_index+1)) - ((0.5 / (dt/dx)) * (u[x_index+1] - u[x_index]))


# Initial condition
def initialise(u, condition):

    print("Initial Condition: U_0(x_j) = ")

    # Something
    if condition == 1:
        pass

    # Something 1
    elif condition == 2:
        pass
    
    # Something 2
    elif condition == 3:
       pass 

    u_ex(condition) # Compute exact solution based on initial condition and problem statement


# Simulating the time stepping of PDE using Finite Volume Method and Flux Scheme
def simulate(u_n):
    pass

# Calculates the TV bound
def calculate_tv(u):
    tv_norm = 0.0
    for j in range(first_cell, last_cell):
        tv_norm += abs(u[j] - u[j-1])
    
    return tv_norm


# Calculates exact Solution
def u_ex(condition):
    ex_file = open("uex.txt", "w")  # Exact Solution data dump file

    ex_file.close()

