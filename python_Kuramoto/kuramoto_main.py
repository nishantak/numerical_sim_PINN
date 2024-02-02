from scheme import *
from config import *

U = np.zeros(Nx)  # U(x);

get_param()

initialise(U, 1)
simulate(U)

print("Total Variation:", calculate_tv(U), "\n")

plot()




