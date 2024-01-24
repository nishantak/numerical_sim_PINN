#include <iostream>
#include <vector>
#include <math.h>
#include "numSim.h"

using namespace std;

// Simulation Parameters 
double xmin = 0, xmax = 40;  // Domain limits
double L = abs(xmax- xmin);   //Domain Length
int ghost_cells = 4;    // Number of ghost cells 
int Nx = 100 + ghost_cells;   // Number of spatial points
long double dx = L/(Nx-1);  // Cell width 
// long double dx = 0.2; // Cell Width
// int Nx = L/dx + ghost_cells; // Number of Spatial Points

double cfl = 0.01;  // Stability Parameter - CFL Number 
long double c = 1.0;  // Wave Velocity

long double dt = cfl * dx / c;  // Time step
double Tf = 2.0;         // Final time / Total Time
int Nt = (int)(Tf/dt);  // No. of time steps

int first_cell = 3, last_cell = Nx-2;   // j domain Limits


// Returns Flux, u^2 / 2
long double flux(long double u){
    return 0.5*u*u;
};


//Driver Code 
int main(){

    vector<long double> U(Nx, 0); // U(x);
    //vector<long double> U_0(U); // Copy of intial u_0, For some reason?

    get_param();

    /*
        initialise() function takes 2 inputs, the vector U and and an integer denoting the 
        initial condition. The Initial Conditions are as follows: -

            1 : U_0(x_j) = sin(x_j+1/2) || U_0(x_j) = -cos(x_j+1/2)

            2 : Discrete initial data, U_0(x_j) = (x_i > 0) ? 0 : 1

            3 : U_0(x_j) = 0.25 * ( sech(sqrt(0.5)/2 * x -7) )^2 

        Similarly, simulate() function takes 2 inputs, the vector U and an integer denoting the
        flux scheme. The Flux Schemes are as follows: -

            1: Lax-Friedrich

            2: Lax-Wendroff
    */ 
    initialise(U, 3);
    simulate(U, 1);

    cout << "Total Variation: " << calculate_tv(U) << endl << endl;
    
    plot();

    return 0;
}


