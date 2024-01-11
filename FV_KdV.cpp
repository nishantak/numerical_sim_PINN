#include<bits/stdc++.h>
using namespace std;

//Function Signatures
void intialise(vector<double>&);
void simulate(vector<double>&);
void get_param();
void plot();
void write_data(string, vector<double>);
double calculate_tv(vector<double>&);
// Function Signatures


// Simulation Parameters
int Nx = 201;      // Number of points
double xmin = 0.0, xmax = 2.0;  // Domain limits
double L = abs(xmax- xmin);   // Domain Length
double dx = L/(Nx-1);  // Cell width

double cfl = 0.5; // Stability Parameter - CFL Number 
double c=1; // Wave Velocity

double dt = cfl * dx / c;     // Time step
double Tf = 2.0;         // Final time / Total Time
int Nt = (int)(Tf/dt);  // No. of time steps


vector<double> U(Nx, 0); // U(x);
vector<double> rhs(Nx, 0); // RHS of PDE


// Returns Flux, u^2 / 2
double flux(double u){
    return 0.5*u*u;
};


// Numerical Flux Scheme (Lax-Friedrich)
double num_flux(double u, double u_next){
    return ( (0.5 * (flux(u) + flux(u_next))) - (0.5 / (dt/dx) * (u_next - u)) );
}


//Driver Code 
int main(){
    intialise(U);
    simulate(U);

    get_param();

    cout << "Total Variation: " << calculate_tv(U) << endl;
    
    plot();

    return 0;
}


/// @brief initialsies with intial condition
void intialise(vector<double> &u){
    for(int j=0; j<Nx; j++)
        u[j] = cos(M_PI*(xmin + j*dx));
}


/// @brief Simulating the time stepping of PDE using Finite Volume Method and Flux Scheme
void simulate(vector<double> &u){

    double t=0; //Current Time
    while(t<=Tf){
        write_data("simulation.txt", U); // Write initial data, Time Step 0
        
        for(int j=1; j<=Nx-2; j++){
            
            double third_derivative = (u[j] - 3.0*u[j-1] + 3.0*u[j-2] - u[j-3]) / (dx*dx*dx); 

            // Update using Numerical Scheme
            u[j] -= ( (dt/dx)*(num_flux(u[j], u[j+1]) - num_flux(u[j-1], u[j])) + dt*third_derivative );
        }
        
        // Boundary Conditions
        u[0] = u[1]; // Left Boundary
        u[Nx-1] = u[Nx-2]; // Right Boundary

        write_data("simulation.txt", U); // Write Simulation Data for THIS time step

        t+=dt;

    } write_data("U_final.txt", U); // Write Simulation Data for FINAL time step

}


/// @brief Plots the graph using python matplotlib script
void plot(){
    system("python matplot.py");
}


/// @brief Calculates the TV bound
double calculate_tv(vector<double>& u) {
    double tv_norm = 0.0;
    for (int j = 1; j < Nx; j++) {
        tv_norm += abs(u[j] - u[j - 1]);
    }
    
    return tv_norm;
}


/// @brief write data to file
void write_data(string filename, vector<double> u){
    ofstream out_file(filename);
    for(int i=0; i<Nx; i++)
        out_file << u[i] << " ";
    out_file << endl;
}


/// @brief Print Simulation Parameters
void get_param(){
    cout << "Number of Spatial Points (Nx): " << Nx << endl;
    cout << "Domain Limits (xmin, xmax): " << xmin << ", " << xmax << endl;
    cout << "Domain Length (L): " << L << endl;
    cout << "Cell Width (dx): " << dx << endl << endl;
    cout << "Stability Parameter (CFL Number): " << cfl << endl << endl;
    cout << "Wave Velocity (c): " << c << endl << endl;
    cout << "Time Step (dt): " << dt << endl;
    cout << "Final Time (Tf): " << Tf << endl;
    cout << "Number of Time Steps (Nt): " << Nt << endl << endl;
}