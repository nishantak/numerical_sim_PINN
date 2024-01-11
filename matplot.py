import matplotlib.pyplot as plt
import numpy as np

# Read simulation data from the file
data = np.loadtxt('simulation_data.txt')

# Plotting
x_values = np.linspace(0, 2.0, 200)
time_steps = len(data)
plt.xlabel('x')
plt.ylabel('u(x)')

for t in range(time_steps):
    plt.clf()
    plt.title(f'FVM KdV Simulation at Time Step {t}')
    plt.plot(x_values, data[t, :])
    plt.draw()
    plt.pause(0.1)  # Adjust the pause time as needed

plt.show()