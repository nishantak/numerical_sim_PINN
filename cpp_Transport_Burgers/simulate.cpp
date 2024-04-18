#include <iostream>
#include <vector>
#include <math.h>
#include "FV_sim.h" // Custom finite volume PDE solver library

using namespace std;

// Simulation Parameters 
int equation; // Problem Statement

double xmin = -2*M_PI, xmax = 2*M_PI ;  // Domain limits
double L = abs(xmax- xmin);   //Domain Length
int ghost_cells = 0;    // Number of ghost cells 
int Nx = 512 + ghost_cells;   // Number of spatial points
long double dx = L/(Nx-1);  // Cell width 
// long double dx = 0.2; // Cell Width
// int Nx = (int)(L/dx) + 1 + ghost_cells; // Number of Spatial Points

double cfl = 0.5;  // Stability Parameter - CFL Number 
long double c = 1.0;  // Wave Velocity

long double dt;    // Time step
double Tf = 2.0;         // Final time / Total Time

int first_cell = 0, last_cell = Nx-1;   // j domain Limits


// Returns Flux, f(u)
float a;  // Constant Flux Multiplier
long double flux(long double u){
    switch (equation){
        // flux = a*u | Transport Equation
        case 1:
            return a*u;
        break;

        // flux = u^2 / 2 | Burger's Equation
        case 2:
         return 0.5*u*u;

        default: break;
    }   
}


//Driver Code 
int main(){

    vector<long double> U(Nx, 0); // U(x);
    
    get_param();

    /*
        "equation" is a choice variable for user to choose which equation to simulate. The choices are 
        as follows: -

            1 : Transport Equation

            2 : Burger's Equation


        initialise() function takes 2 inputs, the vector U and and an integer denoting the 
        initial condition. The Initial Conditions are as follows: -

            1 : U_0(x_j) = sin(x_j+1/2) || U_0(x_j) = -cos(x_j+1/2)

            2 : Discrete initial data, U_0(x_j) = (x_i > 0) ? 0 : 1

            # 3 : U_0(x_j) = 0.25 * ( sech(sqrt(0.5)/2 * x -7) )^2  |  KdV stuff (nvm)


        Similarly, simulate() function takes 3 inputs, the vector U, an integer denoting the flux scheme,
        and an integer denoting the boundary conditons (to be set according to the initial conditions). 
        
        The Flux Schemes are as follows: -

            1: Lax-Friedrich

            2: Lax-Wendroff
        
        The boundary conditions are as follows: -
        
            2 :  u_n_plus1[Nx-1] = u_n_plus1[Nx-2]; // RIGHT Boundary
                 u_n_plus1[0] = u_n_plus1[1];   // LEFT Boundary
            
            1 (default) :   u_n_plus1[Nx-1] = u_n_plus1[Nx-2]; 
                            u_n_plus1[0] = u_n_plus1[Nx-1]; 
    */ 
    
    equation = 1; 

    a = 1;  // Constant Flux Multiplier

    initialise(U, 1);
    simulate(U, 1, 1);

    cout << "Total Variation: " << calculate_tv(U) << endl << endl;
    
    plot();

    return 0;
}


