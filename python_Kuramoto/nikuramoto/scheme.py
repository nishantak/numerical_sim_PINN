from scipy.fft import fft, ifft
from .functions import *

'''
    Finite Volume Solver for Kuramoto Equation with NON-IDENTICAL oscillators
'''
# Returns Lax-Friedrich Numerical Flux 
def num_flux(L_U, u, y_index, dt):
    return 0.5 * (flux(L_U, u, j, y_index,) + flux(L_U, u, j_plus1, y_index)) - ((0.5 / (dt/dx)) * (u[y_index, j_plus1] - u[y_index, j]))


# Initial condition
def initialise(u, condition):

    print("Initial Condition: U_0(x_j) = ")

    # Polynomial Initial data
    if condition == 1:
        print("((thet >= pi/4) && (thet < pi/2) && (om >= 0) && (om <= 1)) * (64/3*pi^2) * thet*om; \n\n 0 , else")
        for k in range(Ny):
            for j in range(Nx):
                u[k, j] = ((X[k, j] >= np.pi/4) & (X[k, j] < np.pi/2) & (Y[k, j] >= 0) & (Y[k, j] <=1)) * (64/(3*(np.pi**2))) * X[k, j]*Y[k, j]

    # Something 1
    elif condition == 2:
        pass

    u_ex(condition) # Compute exact solution based on initial condition and problem statement


# Simulating the time stepping of PDE using Finite Volume Method and Flux Scheme
def simulate(u_n):
    # Output dump files
    out_file = open("simulation_data.txt", "w")
    fin_file = open("U_final.txt", "w")
    # debug_file = open("all_data.txt", "w")
    
    t = 0  # Current Time

    write_data(out_file, u_n, x_first_cell, x_last_cell, y_first_cell, y_last_cell)  # Write initial data
    # write_data(debug_file, u_n, 0, Nx-1, 0, Ny-1)

    # Time stepping loop
    while t < Tf:
        # Compute convolution (integration)
        rho = dy * np.sum(u_n, axis=0)
        L = np.real(ifft(fft(np.sin(x0)) * fft(rho))) * dx
        L_V = np.tile(L, (Ny, 1))
        L_U = L_V*u_n

        # Setting time step
        dt = cfl * dx / np.max(np.abs(K*L_V-Y))
        if(t+dt > Tf): dt = Tf-t
        
        # Numerical Flux
        F= np.zeros_like(u_n)
        for k in range(Ny):
            F[k, j] = num_flux(L_U, u_n, k, dt)

        # Update using Numerical Scheme (Vectorized operations)
        u_n[:, j_plus1] = u_n[:, j] - (dt/dx) * (F[:, j] - F[:, j_min1])
        
        write_data(out_file, u_n, x_first_cell, x_last_cell, y_first_cell, y_last_cell)  # Write Simulation Data for THIS time step
        # write_data(debug_file, u_n, 0, Nx-1, 0, Ny-1)

        t += dt

    write_data(fin_file, u_n, x_first_cell, x_last_cell, y_first_cell, y_last_cell)  # Write Simulation Data for FINAL time step
    out_file.close(), fin_file.close()


# Calculates the TV bound
def calculate_tv(u):
    tv_norm = 0.0
    for k in range(y_first_cell, y_last_cell+1):
        for i in range(x_first_cell+1, x_last_cell+1):
            tv_norm += abs(u[k, i] - u[k, i-1])
    
    return tv_norm


# Calculates exact Solution
def u_ex(condition):
    ex_file = open("uex.txt", "w")  # Exact Solution data dump file

    ex_file.close()
