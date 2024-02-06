from idkuramoto import *

U = np.zeros(Nx)  # U(x);

get_param()

'''
    initialise() function takes 2 inputs, the vector U and and an integer denoting the initial condition. 
    The Initial Conditions are as follows: -

        1 : U_0(x_j) = 1/4 * ((x_j >= 3*pi/4) && (x_j <= 5*pi/4)) + 1/2 * ((x_j >= pi/2) && (x_j <= 3*pi/2)
        
        2 : Something

        3 : Something 2

'''
initialise(U, 1)
simulate(U)

print("Total Variation:", calculate_tv(U), "\n")

plot()




