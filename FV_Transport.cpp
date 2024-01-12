#include<bits/stdc++.h>
using namespace std;

//Function Signatures
void simulate(vector<long double>&);
void get_param();
void plot();
void write_data(ofstream&, vector<long double>, int, int);
long double calculate_tv(vector<long double>&);
// Function Signatures


// Simulation Parameters
int Nx = 500;   // Number of spatial points 
int ghost_cells = 2; // Number of ghost cells 
double xmin = 0, xmax = 2*M_PI;  // Domain limits
double L = abs(xmax- xmin);   // Domain Length
long double dx = L/(Nx-1);  // Cell width

double cfl = 0.5; // Stability Parameter - CFL Number 
double c=1; // Wave Velocity

long double dt = cfl * dx / c;     // Time step
double Tf = 2.0;         // Final time / Total Time
int Nt = (int)(Tf/dt);  // No. of time steps


// Returns Flux, 2u
long double flux(long double u){
    return 0.5*u*u;
};


// Returns Numerical Flux, Numerical Flux Scheme (Lax-Friedrich)
long double num_flux(long double u, long double u_next){
    return 0.5 * (flux(u) + flux(u_next)) - ((0.5 / (dt/dx)) * (u_next - u));
}


/// @brief initialise with intial condition, U_0(x) = sin(x)
void intialise(vector<long double> &u){
    for(int j=0; j<Nx; j++)
        u[j] = sin((xmin + (j+0.5)*dx));
}


//Driver Code 
int main(){

    Nx += ghost_cells;
    vector<long double> U(Nx, 0); // U(x);
    vector<long double> U_0(U.begin(), U.end()); // Copy of intial u_0 | For some reason

    intialise(U);
    simulate(U);

    get_param();

    cout << "Total Variation: " << calculate_tv(U) << endl;
    
    plot();

    return 0;
}


/// @brief Simulating the time stepping of PDE using Finite Volume Method and Flux Scheme
/// @param u_n This Current Time Step data
void simulate(vector<long double> &u_n){
    // Output dump files
    ofstream out_file("simulation_data.txt");
    ofstream fin_file("U_final.txt");
    ofstream debug_file("all_data.txt");

    int first_cell = 1, last_cell = Nx-2; // Spatial domain Limits (Discretized)
    double t=0; //Current Time
    
    write_data(out_file, u_n, first_cell, last_cell); // Write initial data
    write_data(debug_file, u_n, 0, Nx-1);
    
    // Time Stepping Loop
    while(t<=Tf){
        
        vector<long double> u_n_plus1(u_n.begin(), u_n.end()); // Next Time Step

        for(int j=first_cell; j<=last_cell; j++){
            
            // Numerical Flux
            long double F_j_plus_half = num_flux(u_n[j], u_n[j+1]);
            long double F_j_min_half = num_flux(u_n[j-1], u_n[j]);
            
            // Update using Numerical Scheme
            u_n_plus1[j] -= (dt/dx) * (F_j_plus_half - F_j_min_half);
        }

        // Boundary conditions        
        u_n_plus1[Nx-1] = u_n_plus1[Nx-2]; // RIGHT Boundary
        u_n_plus1[0] = u_n_plus1[Nx-1]; // LEFT Boundary

        for(int i=0; i<=Nx-1; i++)
            u_n[i] = u_n_plus1[i]; // Store u^n+1_j in u^n_j for next time step

        write_data(out_file, u_n, first_cell, last_cell); // Write Simulation Data for THIS time step
        write_data(debug_file, u_n, 0, Nx-1);

        t+=dt;

    } write_data(fin_file, u_n, 1, Nx-2); // Write Simulation Data for FINAL time step

out_file.close(); fin_file.close();

}


/// @brief Calculates the TV bound
long double calculate_tv(vector<long double>& u) {
    double tv_norm = 0.0;
    for (int j = 1; j < Nx; j++) {
        tv_norm += abs(u[j] - u[j - 1]);
    }
    
    return tv_norm;
}


/// @brief Pass plot data to python matplotlib script and plots graph
void plot(){
    // Set environment variables to pass to python script
    setenv("xmin", to_string(xmin).c_str(), 1);
    setenv("xmax", to_string(xmax).c_str(), 1);
    setenv("Nx", to_string(Nx - ghost_cells).c_str(), 1);
    // CLI command to run python script
    system("python matplot.py");
}


/// @brief write data to file using file stream
void write_data(ofstream& filename, vector<long double> u, int start, int end){
    for(int i=start; i<=end; i++)
        filename << u[i] << " ";
    filename << endl;
}


/// @brief Print Simulation Parameters
void get_param(){
    cout << endl << "Number of Spatial Points (Nx): " << Nx - ghost_cells << endl;
    cout << "Domain Limits (xmin, xmax): " << xmin << ", " << xmax << endl;
    cout << "Domain Length (L): " << L << endl;
    cout << "Cell Width (dx): " << dx << endl << endl;
    cout << "Stability Parameter (CFL Number): " << cfl << endl << endl;
    cout << "Wave Velocity (c): " << c << endl << endl;
    cout << "Time Step (dt): " << dt << endl;
    cout << "Final Time (Tf): " << Tf << endl;
    cout << "Number of Time Steps (Nt): " << Nt << endl << endl;
}