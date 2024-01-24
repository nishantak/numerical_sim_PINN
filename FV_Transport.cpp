#include<bits/stdc++.h>
using namespace std;

// Function Signatures
void simulate(vector<long double>&, int);
void uex(int);
void get_param();
void plot();
void write_data(ofstream&, vector<long double>, int, int);
long double calculate_tv(vector<long double>);
// Function Signatures


// Simulation Parameters
double xmin = 0, xmax = 40;  // Domain limits
double L = abs(xmax- xmin);   // Domain Length
int ghost_cells = 2;    // Number of ghost cells 
// int Nx = 200 + ghost_cells;   // Number of spatial points 
// long double dx = L/(Nx-1);  // Cell width
long double dx = 0.01; // Cell Width
int Nx = L/dx + ghost_cells; // Number of Spatial Points

double cfl = 0.5;   // Stability Parameter - CFL Number 
long double c = 1.0;  // Wave Velocity

long double dt = cfl * dx / c;  // Time step
double Tf = 2.0;         // Final time / Total Time
int Nt = (int)(Tf/dt);  // No. of time steps

int first_cell = 1, last_cell = Nx-2;   // j domain Limits


// Returns Flux, 2u
long double flux(long double u){
    return 2*u;
};


// Returns Numerical Flux, Numerical Flux Scheme 
long double num_flux(long double u, long double u_next, int scheme){
    switch(scheme){
        // Lax-Friedrich
        case 1:
            return 0.5 * (flux(u) + flux(u_next)) - ((0.5 / (dt/dx)) * (u_next - u));
            break;

        // Lax-Wendroff
        case 2:
            return u;
            break;
        
        default: break;
    
    }
    
}


/// @brief initialise with intial condition
void intialise(vector<long double> &u, int condition){
    
    cout << "Initial Condition: U_0(x_j) = ";
    switch (condition){

            // U_0(x_j) = sin(x_j+1/2) || U_0(x_j) = -cos(x_j+1/2)
            case 1:
                cout << "sin(x_j+1/2)" << endl << endl;
                
                for(int j=0; j<Nx; j++)
                    u[j] = sin((xmin + (j+0.5)*dx)); 
                    //u[j] = -cos((xmin + (j+0.5)*dx));
                
                break;

            // Discrete initial data, x_i > 0 ? 0 : 1
            case 2:
                cout << "x_i > 0 ? 0 : 1" << endl << endl;
                
                for(int j=0; j<Nx; j++)
                    u[j] = xmin + (j+0.5)*dx > 0 ? 0 : 1;
                
                break;
            
            // U_0(x_j) = 0.25 * ( sech(sqrt(0.5)/2 * x -7) )^2 
            case 3:
                cout << "0.25 * (sech(sqrt(0.5)/2 * x -7))^2" << endl << endl;
                
                for(int j=0; j<Nx; j++)
                    u[j] = 0.25 * (pow(1.0/cosh(sqrt(0.5)/2 * (xmin + (j+0.5)*dx) - 7), 2));
                
                break;

            default: break;
    } 
    uex(condition); // Compute exact solution based on initial condition
}


//Driver Code 
int main(){

    vector<long double> U(Nx, 0); // U(x);
    //vector<long double> U_0(U); // Copy of intial u_0, For some reason?

    get_param();

    intialise(U, 1);
    simulate(U, 1);

    cout << "Total Variation: " << calculate_tv(U) << endl << endl;
    
    plot();

    return 0;
}


/// @brief Simulating the time stepping of PDE using Finite Volume Method and Flux Scheme
/// @param u_n Current Time Step data U^n_j
void simulate(vector<long double> &u_n, int flux_scheme){
    // Output dump files
    ofstream out_file("simulation_data.txt");
    ofstream fin_file("U_final.txt");
    //ofstream debug_file("all_data.txt");

    double t=0; //Current Time
    
    write_data(out_file, u_n, first_cell, last_cell); // Write initial data
    //write_data(debug_file, u_n, 0, Nx-1);
    
    // Time Stepping Loop
    while(t<Tf){
        
        vector<long double> u_n_plus1(u_n); // Next Time Step, U^n+1_j, initialised with U^n_j
        for(int j=first_cell; j<=last_cell; j++){
            
            // Numerical Flux
            long double F_j_plus_half = num_flux(u_n[j], u_n[j+1], flux_scheme);
            long double F_j_min_half = num_flux(u_n[j-1], u_n[j], flux_scheme);
            
            // Update using Numerical Scheme
            u_n_plus1[j] -= (dt/dx) * (F_j_plus_half - F_j_min_half);
        }

        // Boundary conditions        
        u_n_plus1[Nx-1] = u_n_plus1[Nx-2]; // RIGHT Boundary
        u_n_plus1[0] = u_n_plus1[Nx-1];  // LEFT Boundary
        
        // Store u^n+1_j in u^n_j for next time step
        u_n = u_n_plus1; 

        write_data(out_file, u_n, first_cell, last_cell); // Write Simulation Data for THIS time step
        //write_data(debug_file, u_n, 0, Nx-1);

        t+=dt;

    } write_data(fin_file, u_n, 1, Nx-2); // Write Simulation Data for FINAL time step

out_file.close(); fin_file.close();

}


///@brief Calcuates exact Solution
///@param condition Initial Condition
void uex(int condition){
    ofstream ex_file("uex.txt"); // Exact Solution data dump file
    switch (condition){
        case 3:
            // Writing exact solution to file
            for(int j=first_cell; j<=last_cell; j++)
                ex_file << 0.25 * (pow(1.0/cosh(sqrt(0.5)/2.0 * (xmin + (j+0.5)*dx - 2.5) - 7.0), 2)) << " ";

            break;

        default: break;
    } ex_file.close();
}


/// @brief Calculates the TV bound
long double calculate_tv(vector<long double> u) {
    double tv_norm = 0.0;
    for (int j = first_cell; j <= last_cell + 1; j++) {
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
    filename.precision(8);
    for(int i=start; i<=end; i++)
        filename << u[i] << " ";
    filename << endl;
}


/// @brief Print Simulation Parameters
void get_param(){
    cout << endl << "Domain Limits (xmin, xmax): " << xmin << ", " << xmax << endl;
    cout << "Domain Length (L): " << L << endl;
    cout << "Number of Spatial Points (Nx): " << Nx - ghost_cells << endl;
    cout << "Cell Width (dx): " << dx << endl << endl;
    cout << "Stability Parameter (CFL Number): " << cfl << endl << endl;
    cout << "Wave Velocity (c): " << c << endl << endl;
    cout << "Final Time (Tf): " << Tf << endl;
    cout << "Time Step (dt): " << dt << endl;
    cout << "Number of Time Steps (Nt): " << Nt << endl << endl;
}