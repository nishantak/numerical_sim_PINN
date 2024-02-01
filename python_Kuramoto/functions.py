import scipy as sp
import matplotlib.pyplot as plt
import os
from config import *


# Plots graph
def plot():
    # Read simulation data from output dump file
    sim_data = np.loadtxt('simulation_data.txt')
    
    ex_data_exists = 0
    if os.path.getsize("uex.txt"):
        ex_data = np.loadtxt('uex.txt')
        ex_data_exists = 1

    # Plot
    x = np.linspace(xmin, xmax, Nx-ghost_cells)
    time_steps = len(sim_data)
    plt.xlabel('x')
    plt.ylabel('u(x)')

    for t in range(time_steps):
        plt.clf()
        plt.title(f'FVM Simulation after Time Step {t}')

        # Exact Soltuion | Will give error if ex_data file is not according to initial data; in such case comment it
        if ex_data_exists:
            plt.plot(x, ex_data, linestyle=":", marker="o", markersize=1, markerfacecolor='none', label='Exact Solution') 
        
        # Numerical Solution
        plt.plot(x, sim_data[t, :], linestyle=":", marker="o", markersize=1, markerfacecolor='none', label='Numerical Solution') 
        
        plt.legend()
        plt.draw()
        plt.pause(0.05)

    plt.show()


# write data to file
def write_data(filename, u, start, end):
    filename.write(" ".join(map(str, u[start : end+1])) + "\n")


# Print Simulation Parameters
def get_param():
    print("\nDomain Limits (xmin, xmax):", xmin, ", ", xmax)
    print("Domain Length (L):", L)
    print("Number of Spatial Points (Nx):", Nx - ghost_cells)
    print("Cell Width (dx):", dx, "\n")
    print("Stability Parameter (CFL Number):", cfl, "\n")
    print("Wave Velocity (c):", c, "\n")
    print("Final Time (Tf):", Tf)
    print("Time Step (dt):", dt)
    print("Number of Time Steps (Nt):", Nt, "\n")

