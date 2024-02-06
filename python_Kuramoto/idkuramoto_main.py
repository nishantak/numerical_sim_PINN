import idkuramoto as ik
import nikuramoto as nk

U = np.zeros(Nx)  # U(x);

equation = 1

# Identical Natural Frequencies
if(equation == 1):
    ik.get_param()

'''
    ik.initialise() function takes 2 inputs, the vector U and and an integer denoting the initial condition. 
    The Initial Conditions are as follows: -

        1 : U_0(x_j) = 1/4 * ((x_j >= 3*pi/4) && (x_j <= 5*pi/4)) + 1/2 * ((x_j >= pi/2) && (x_j <= 3*pi/2)
        
        2 : Something

        3 : Something 2

'''
    ik.initialise(U, 1)
    ik.simulate(U)

# NON-Identical Natural Frequencies 
else if(equation == 2):
    nk.get_param()

'''
    nk.initialise() function takes 2 inputs, the vector U and and an integer denoting the initial condition. 
    The Initial Conditions are as follows: -

        1 : Something 1

'''
    nk.initialise(U, 1)
    nk.simulate(U)



print("Total Variation:", ik.calculate_tv(U), "\n")

ik.plot()

# functions ik.plot() and nk.plot() are identical to each other, so are ik.calculate_tv() and nk.calculate-tv(). Hence, it becomes irrelevant from which module they are imported.



