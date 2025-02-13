'''
    'equation' is a choice variable denoting the equation to simulate.
        
        1 : idkuramoto

        2 : nikuramoto
''' 
equation = 1

# Identical Natural Frequencies
if(equation == 1): 
    from idkuramoto import * 
    U = np.zeros(Nx)  # U(x);
    
# NON-Identical Natural Frequencies 
elif(equation == 2): 
    from nikuramoto import *
    U = np.zeros((Ny, Nx))  # U(y, x);
    

get_param()

'''
    initialise() function takes 2 inputs, the vector U and and an integer denoting the initial condition. 

    For IDkuramoto the Initial Conditions are as follows: -
        
        # Singular Initial data
        1 : U_0(x_j) = 1/4 * ((x_j >= 3*pi/4) && (x_j <= 5*pi/4)) + 1/2 * ((x_j >= pi/2) && (x_j <= 3*pi/2)
        
        # Polynomial Initial data
        2 : U_0(x_j) = (6/pi^3) * (3*pi/2 - x) * (x - pi/2) , if pi/2 <= x < 3*pi/2;
                     = 0 , else

        
    For NIkuramoto the Initial Conditions are as follows: -
        # Polynomial initial data
        1 : U_0(x_j) = ((thet >= pi/4) && (thet < pi/2) && (om >= 0) && (om <= 1)) * (64/3*pi^2) * thet*om ;
                     = 0 ,  else
'''
initialise(U, 2)
simulate(U)

print("\nTotal Variation:", calculate_tv(U), "\n")

plot()
