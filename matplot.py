import matplotlib.pyplot as plt
import numpy as np
import os

# Read simulation data from output dump file
sim_data = np.loadtxt('simulation_data.txt')
ex_data = np.loadtxt('uex.txt')

#Read environment variables set by C++ code
xmin = float(os.getenv("xmin"))
xmax = float(os.getenv("xmax"))
Nx = int(os.getenv("Nx"))

# Plot
x = np.linspace(xmin, xmax, Nx)
time_steps = len(sim_data)
plt.xlabel('x')
plt.ylabel('u(x)')

for t in range(time_steps):
    plt.clf()
    plt.title(f'FVM Simulation after Time Step {t}')

    # Exact Soltuion
    plt.plot(x, ex_data, linestyle=":", marker="o", markersize=1, markerfacecolor='none', label='Exact Solution') 
    # Numerical Solution
    plt.plot(x, sim_data[t, :], linestyle=":", marker="o", markersize=1, markerfacecolor='none', label='Numerical Solution') 
    
    plt.legend()
    plt.draw()
    plt.pause(0.03)

plt.show()