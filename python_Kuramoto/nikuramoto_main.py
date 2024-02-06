from nikuramoto import *

U = np.zeros(Nx)  # U(x);

get_param()

'''
    initialise() function takes 2 inputs, the vector U and and an integer denoting the initial condition. 
    The Initial Conditions are as follows: -

        1 : 

'''
initialise(U, 1)
simulate(U)

print("Total Variation:", calculate_tv(U), "\n")

plot()