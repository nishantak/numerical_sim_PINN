from scipy.fft import fft, ifft
from .functions import *

'''
    Finite Volume Solver for Kuramoto Equation with IDENTICAL oscillators
'''

# Returns Lax-Friedrich Numerical Flux 
def num_flux(V_U, u, x_index, dt):
    if x_index == Nx-1: x_index = 0 # Right Boundary, F^n_Nx-1/2 = F^n_1/2
    return 0.5 * (flux(V_U, u, x_index) + flux(V_U, u, x_index+1)) - ((0.5 / (dt/dx)) * (u[x_index+1] - u[x_index]))


# Initial condition
def initialise(u, condition):

    print("Initial Condition: U_0(x_j) = ")

    # Singular Initial data
    if condition == 1:
        print("1/4 * (delta_3pi/4 + delta_5pi/4) + 1/2 * X_[pi/2, 3pi/2] ")
        def dirac_delta(x, a, epsilon=0.01):
            return np.where(np.abs(x-a) < epsilon, 1/(2*epsilon), 0)

        def heaviside(x):
            return np.where(x>=0, 1, 0)

        u[:] = ( 1/4 * (dirac_delta(x, 3*np.pi/4) + dirac_delta(x, 5*np.pi/4)) +
                1/2 * heaviside(x - np.pi/2) * (1 - heaviside(x - 3*np.pi/2)) )


    # Polynomial Initial data
    elif condition == 2:
        print("(6/pi^3) * (3*pi/2 - x) * (x - pi/2) , if pi/2 <= x < 3*pi/2; \n\n 0 , else")
        for j in range(Nx):
            if (x[j] >= np.pi/2) and (x[j] < 3*np.pi/2):
                u[j] = (6/np.pi**3) * (3*np.pi/2 - x[j]) * (x[j] - np.pi/2)
            else: 0
    

    # Piecewise Initial Data
    elif condition == 3:
        print("2/(3pi) if pi/2 <= x <= 3pi/2; 1/(3pi), else\n")
        mask = (x >= (np.pi/2)) & (x <= (3*np.pi/2))
        u[mask]  = 2.0 / (3.0*np.pi)
        u[~mask] = 1.0 / (3.0*np.pi)


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

            # Update using Numerical Scheme (Two array method)
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
    for j in range(first_cell+1, last_cell+1):
        tv_norm += abs(u[j] - u[j-1])
    
    return tv_norm


# Calculates exact Solution
def u_ex(condition):
    ex_file = open("uex.txt", "w")  # Exact Solution data dump file

    ex_file.close()

