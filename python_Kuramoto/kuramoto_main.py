from scheme import *
from config import *

U = np.zeros(Nx)  # U(x);

get_param()

initialise(U, 2)
simulate(U, 2)

print("Total Variation:", calculate_tv(U), "\n")

plot()




