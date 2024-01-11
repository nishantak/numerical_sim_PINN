#include<bits/stdc++.h>
using namespace std;

//Function Signatures
void intialise(vector<double>&);
void simulate(vector<double>&);
void plot();
// Function Signatures

// Simulation Parameters
int Nx = 201;      // Number of points
double xmin = 0.0, xmax = 2.0;  // Domain limits
double L = abs(xmax- xmin);   // Domain Length
double dx = L/(Nx-1);  // Cell width

double cfl = 0.75; // Staility Parameter - CFL Number 
double c=1; // Wave Velocity

double dt = cfl * dx / c;     // Time step
double Tf = 2.0;         // Final time / Total Time
int Nt = (int)(Tf/dt);  // No. of time steps

vector<double> U(Nx, 0); // U(x);
vector<double> rhs(Nx, 0); // RHS of PDE

double flux(double u){
    return 2*u;
};

// Numerical Flux Scheme (Lax-Friedrich)
double num_flux(double u, double u_next){
    return ( (0.5 * (flux(u) + flux(u_next))) - (0.5 / (dt/dx) * (u_next - u)) );
}


//Driver Code 
int main(){
    intialise(U);
    simulate(U);
    plot();

    return 0;
}

/// @brief initialsies with intial condition
void intialise(vector<double> &u){
    for(int j=0; j<Nx; j++)
        u[j] = cos(M_PI*( xmin + j*dx));
}

/// @brief Simulating the time stepping of PDE using Finite Volume Method and Flux Scheme
void simulate(vector<double> &u){
    ofstream outFile("simulation_data.txt");
    double t=0; //Current Time

    while(t<Tf){
        for(int j=0; j<=Nx-2; j++){
            double flux = num_flux(u[j], u[j+1]);
            
            rhs[j] += flux;
            rhs[j+1] -= flux;
            
            if (j!=0) u[j] -= (dt/dx) * rhs[j]; 

            outFile << u[j] << " "; // Write the simulation data for this time step to the file
        }outFile << endl;
        t+=dt;
    }outFile.close();
}

void plot(){
    system("python matplot.py");
}