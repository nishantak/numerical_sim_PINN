from scipy.fft import fft, ifft
from .functions import *

'''
    Finite Volume Solver for Kuramoto Equation with IDENTICAL natural fruequencies
'''

# Returns Lax-Friedrich Numerical Flux 
def num_flux(V_U, u, x_index, dt):
    if x_index == Nx-1: x_index = 0 # Right Boundary
    return 0.5 * (flux(V_U, u, x_index) + flux(V_U, u, x_index+1)) - ((0.5 / (dt/dx)) * (u[x_index+1] - u[x_index]))


# Initial condition
def initialise(u, condition):

    print("Initial Condition: U_0(x_j) = ")

    # U_0(x_j) = 1/4 * ((x_j >= 3*pi/4) && (x_j <= 5*pi/4)) + 1/2 * ((x_j >= pi/2) && (x_j <= 3*pi/2)
    if condition == 1:
        print("1/4 * ((x_j >= 3*pi/4) && (x_j <= 5*pi/4)) + 1/2 * ((x_j >= pi/2) && (x_j <= 3*pi/2)")
        for j in range(Nx):
            u[j] = 0.25*((x[j]>=3*np.pi/4.0) & (x[j]<=5*np.pi/4.0)) + 0.5*((x[j]>=np.pi/2.0) & (x[j]<=3*np.pi/2.0))

    # Something
    elif condition == 2:
        pass
    
    # Something 2
    elif condition == 3:
       pass 

    u_ex(condition) # Compute exact solution based on initial condition and problem statement


# Simulating the time stepping of PDE using Finite Volume Method and Flux Scheme
def simulate(u_n):
    # Output dump files
    out_file = open("simulation_data.txt", "w")
    fin_file = open("U_final.txt", "w")
    # debug_file = open("all_data.txt", "w")

    t = 0  # Current Time

    write_data(out_file, u_n, first_cell, last_cell)  # Write initial data
    # write_data(debug_file, u_n, 0, Nx-1)

    # Time stepping loop
    while t < Tf:

        u_n_plus1 = np.copy(u_n)  # Next Time Step, U^n+1_j | initialised with copy of current time step data, U^n_j
        
        # Convolution term, V[mu^n]_j
        V_U = np.real(ifft(fft(np.sin(x)) * fft(u_n))) * dx

        # Setting time step
        dt = cfl * dx / max(abs(V_U))
        if(t+dt > Tf): dt = Tf-t

        for j in range(first_cell, last_cell+1):
            
            # Numerical Flux
            F_j_plus_half = num_flux(V_U, u_n, j, dt)
            F_j_min_half = num_flux(V_U, u_n, j-1, dt)

            # Update using Numerical Scheme
            u_n_plus1[j] -= (dt/dx) * (F_j_plus_half - F_j_min_half)

        u_n_plus1[0] = u_n_plus1[Nx-1]  # LEFT Boundary

        # Store u^n+1_j in u^n_j for next time step
        u_n = u_n_plus1

        write_data(out_file, u_n, first_cell, last_cell)  # Write Simulation Data for THIS time step
        # write_data(debug_file, u_n, 0, Nx-1)
        t += dt

    write_data(fin_file, u_n, first_cell, last_cell)  # Write Simulation Data for FINAL time step
    out_file.close(), fin_file.close()


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

