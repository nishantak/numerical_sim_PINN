#include<bits/stdc++.h>
using namespace std;

//Function Signatures
void simulate(vector<long double>&);
void get_param();
void plot();
void write_data(ofstream&, vector<long double>, int);
long double calculate_tv(vector<long double>&);
// Function Signatures


// Simulation Parameters
int Nx = 204;   // Number of points | Number of ghost cells = 4
double xmin = 0.0, xmax = 2*M_PI;  // Domain limits
double L = abs(xmax- xmin);   // Domain Length
long double dx = L/(Nx-1);  // Cell width

double cfl = 0.5; // Stability Parameter - CFL Number 
double c=1; // Wave Velocity

long double dt = cfl * dx / c;     // Time step
double Tf = 2.0;         // Final time / Total Time
int Nt = (int)(Tf/dt);  // No. of time steps


vector<long double> U(Nx, 0); // U(x);


// Returns Flux, 2u
long double flux(long double u){
    return 2*u;
};


// Returns Numerical Flux, Numerical Flux Scheme (Lax-Friedrich)
long double num_flux(long double u, long double u_next){
    return 0.5 * (flux(u) + flux(u_next)) - (0.5 * (dt/dx) * (u_next - u));
}


/// @brief initialise with intial condition, U_0(x) = sin(x)
void intialise(vector<long double> &u){
    for(int j=0; j<Nx; j++)
        u[j] = sin((xmin + j*dx));
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


/// @brief Simulating the time stepping of PDE using Finite Volume Method and Flux Scheme
void simulate(vector<long double> &u){
    // Output dump files
    ofstream out_file("simulation_data.txt");
    ofstream fin_file("U_final.txt");

    double t=0; //Current Time
    
    write_data(out_file, u, 1); // Write initial data
    // Time Stepping Loop
    while(t<=Tf){
        for(int j=1; j<=Nx-2; j++){
            
            // Numerical Flux
            long double F_j_plus_half = num_flux(u[j], u[j+1]);
            long double F_j_min_half = num_flux(u[j-1], u[j]);

            // Update using Numerical Scheme
            u[j] -= (dt/dx) * (F_j_plus_half - F_j_min_half);
        }
        
        // Boundary conditions
        u[Nx - 1] = u[Nx-2]; // Right Boundary
        u[0] = u[Nx - 1]; // Left Boundary

        write_data(out_file, u, 1); // Write Simulation Data for THIS time step

        t+=dt;

    } write_data(fin_file, u, 1); // Write Simulation Data for FINAL time step

out_file.close(); fin_file.close();

}


/// @brief Pass plot data to python matplotlib script and plots graph
void plot(){
    // Set environment variables to pass to python script
    setenv("xmin", to_string(xmin).c_str(), 1);
    setenv("xmax", to_string(xmax).c_str(), 1);
    setenv("Nx", to_string(Nx).c_str(), 1);
    // CLI run python script
    system("python matplot.py");
}


/// @brief Calculates the TV bound
long double calculate_tv(vector<long double>& u) {
    double tv_norm = 0.0;
    for (int j = 1; j < Nx; j++) {
        tv_norm += abs(u[j] - u[j - 1]);
    }
    
    return tv_norm;
}


/// @brief write data to file using file stream
void write_data(ofstream& filename, vector<long double> u, int start){
    for(int i=start; i<Nx-1; i++)
        filename << u[i] << " ";
    filename << endl;
}


/// @brief Print Simulation Parameters
void get_param(){
    cout << endl << "Number of Spatial Points (Nx): " << Nx << endl;
    cout << "Domain Limits (xmin, xmax): " << xmin << ", " << xmax << endl;
    cout << "Domain Length (L): " << L << endl;
    cout << "Cell Width (dx): " << dx << endl << endl;
    cout << "Stability Parameter (CFL Number): " << cfl << endl << endl;
    cout << "Wave Velocity (c): " << c << endl << endl;
    cout << "Time Step (dt): " << dt << endl;
    cout << "Final Time (Tf): " << Tf << endl;
    cout << "Number of Time Steps (Nt): " << Nt << endl << endl;
}