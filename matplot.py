import matplotlib.pyplot as plt
import numpy as np
import os

# Read simulation data from output dump file
data = np.loadtxt('simulation_data.txt')

#Read environment variables set by C++ code
xmin = float(os.getenv("xmin"))
xmax = float(os.getenv("xmax"))
Nx = int(os.getenv("Nx"))

# Plotting
x_values = np.linspace(xmin, xmax, Nx)
time_steps = len(data)
plt.xlabel('x')
plt.ylabel('U(x)')

for t in range(time_steps):
    plt.clf()
    plt.title(f'FVM Simulation after Time Step {t}')
    plt.plot(x_values, data[t, :], linestyle="solid", marker="o", markersize=1, markerfacecolor='none')
    plt.draw()
    plt.pause(0.05)

plt.show()