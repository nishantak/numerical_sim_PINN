import matplotlib.pyplot as plt
import os
from .config import *

# Plot 3D graph
def plot():
    # Read simulation data from output dump file
    sim_data = np.loadtxt('simulation_data.txt')

    ex_data_exists = os.path.getsize("uex.txt") > 0
    if ex_data_exists:
        ex_data = np.loadtxt('uex.txt')

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x = np.linspace(xmin, xmax, Nx - x_ghost_cells)
    y = np.linspace(ymin, ymax, Ny - y_ghost_cells)
    X, Y = np.meshgrid(x, y)

    time_steps = len(sim_data)
    Dt = Tf / (time_steps-1)

    sim_data = sim_data.reshape(time_steps, Ny - y_ghost_cells, Nx - x_ghost_cells)

    for t in range(time_steps):
        ax.clear()
        ax.set_title(f'FVM Simulation after Time Step {t}, t={round(t*Dt, 3)}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('u(x, y)')
        # Exact solution
        if ex_data_exists:
            ax.plot_surface(X, Y, ex_data[t, :, :])
        # Numerical solution
        ax.plot_surface(X, Y, sim_data[t, :, :])

        plt.draw()
        plt.pause(0.08)

    plt.show()


# write data to file
def write_data(filename, u, start_x, end_x, start_y, end_y):
    for k in range(start_y, end_y+1):
        for j in range(start_x, end_x+1):
            filename.write(str(u[k, j]) + " ")
    filename.write("\n")


# Print Simulation Parameters
def get_param():
    print("\nNon-Identical Kuramoto Equation\n")
    print("X Domain Limits (xmin, xmax):", xmin, ", ", xmax)
    print("X Domain Length (Lx):", Lx)
    print("Number of X Points (Nx):", Nx - x_ghost_cells)
    print("Cell Width (dx):", dx, "\n")
    print("Y Domain Limits (ymin, ymax):", ymin, ", ", ymax)
    print("Y Domain Length (Ly):", Ly)
    print("Number of Y Points (Ny):", Ny - y_ghost_cells)
    print("Cell Height (dy):", dy, "\n")
    print("Stability Parameter (CFL Number):", cfl, "\n")
    print("Coupling Strength (K):", K, "\n")
    print("Final Time (Tf):", Tf, "\n")
    print(X, Y)
