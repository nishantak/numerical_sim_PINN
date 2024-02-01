from functions import *


'''
    Finite Volume Solver for Kuramoto Equation
'''


# Returns Lax-Friedrich Numerical Flux 
def num_flux(u, u_next):
    return 0.5 * (flux(u) + flux(u_next)) - ((0.5 / (dt/dx)) * (u_next-u))


# initialise with initial condition
def initialise(u, condition):

    print("Initial Condition: U_0(x_j) = ")

    # U_0(x_j) = sin(x_j+1/2)
    if condition == 1:
        print("sin(x_j+1/2)")

        for j in range(Nx):
            u[j] = np.sin(xmin + (j+0.5)*dx)

    # Discrete initial data, x_i > 0 ? 0 : 1
    elif condition == 2:
        print("x_i > 0 ? 0 : 1")

        for j in range(Nx):
            u[j] = 0 if xmin + (j+0.5)*dx > 0 else 1

    # U_0(x_j) = 0.25 * ( sech(sqrt(0.5)/2 * x -7) )^2 
    elif condition == 3:
        print("0.25 * ( sech(sqrt(0.5)/2 * x -7) )^2 ")

        for j in range(Nx):
            u[j] = 0.25 * (1.0 / np.cosh(np.sqrt(0.5) / 2 * (xmin + (j+0.5)*dx) - 7))**2

    u_ex(condition) # Compute exact solution based on initial condition and problem statement


# Simulating the time stepping of PDE using Finite Volume Method and Flux Scheme
def simulate(u_n, boundary_condition):
    # Output dump files
    out_file = open("simulation_data.txt", "w")
    fin_file = open("U_final.txt", "w")
    # debug_file = open("all_data.txt", "w")

    t = 0  # Current Time

    write_data(out_file, u_n, first_cell, last_cell)  # Write initial data
    # write_data(debug_file, u_n, 0, Nx-1)

    # Time stepping loop
    while t < Tf:

        u_n_plus1 = np.copy(u_n)  # Next Time Step, U^n+1_j, initialised with U^n_j

        for j in range(first_cell, last_cell+1):
            
            # Numerical Flux
            F_j_plus_half = num_flux(u_n[j], u_n[j+1])
            F_j_min_half = num_flux(u_n[j-1], u_n[j])

            # Update using Numerical Scheme
            u_n_plus1[j] -= (dt/dx) * (F_j_plus_half - F_j_min_half)

        # Boundary Conditions
        if boundary_condition == 2:
            u_n_plus1[Nx - 1] = u_n_plus1[Nx - 2]  # RIGHT Boundary
            u_n_plus1[0] = u_n_plus1[1]  # LEFT Boundary
        else:
            u_n_plus1[Nx - 1] = u_n_plus1[Nx - 2]
            u_n_plus1[0] = u_n_plus1[Nx - 1]

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
    for j in range(first_cell, last_cell + 1):
        tv_norm += abs(u[j] - u[j - 1])
    
    return tv_norm


# Calculates exact Solution
def u_ex(condition):
    
    ex_file = open("uex.txt", "w")  # Exact Solution data dump file
    # Transport equation
    if equation == 1:
        if condition == 1:
            # Writing exact solution to file
            for j in range(first_cell, last_cell + 1):
                ex_file.write(str(np.sin(xmin + (j+0.5)*dx - a*Tf)) + " ")

        elif condition == 2:
            for j in range(first_cell, last_cell + 1):
                u = 0 if xmin + (j+0.5)*dx - a*Tf > 0 else 1
                ex_file.write(str(u) + " ")
    
    # Burger's Equation
    if equation == 2:
        if condition == 3:
            pass

    ex_file.close()

