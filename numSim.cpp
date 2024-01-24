#include "numSim.h"
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

/// @brief Returns Numerical Flux
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
    } uex(condition); // Compute exact solution based on initial condition
}


/// @brief Simulating the time stepping of PDE using Finite Volume Method and Flux Scheme
/// @param u_n Current Time Step data U^n_j
void simulate(vector<long double> &u_n, int flux_scheme){
    // Output dump files
    ofstream out_file("dump_files/simulation_data.txt");
    ofstream fin_file("dump_files/U_final.txt");
    //ofstream debug_file("dump_files/all_data.txt");

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

            // U_xxx Discretization 1: (u^n_j - 3u^n_j-1 + 3u^n_j-2 - u^n_j-3) / dx^3
            long double u_xxx = (u_n[j] - 3.0*u_n[j-1] + 3.0*u_n[j-2] - u_n[j-3]) / (dx*dx*dx);
            
            // U_xxx Discretization 2: (u^n_j+2 - 2u^n_j+1 + 2u^n_j-1 - u^n_j-2) / 2dx^3
            // long double u_xxx = (u_n[j+2] - 2.0*u_n[j+1] + 2.0*u_n[j-1] - u_n[j-2]) / (2.0 * (dx*dx*dx));
            
            // U_xxx Calculation 3: Newton forward difference method
            // vector<long double> u_xxx(Nx,0);
            // u_xxx = third_derivative(u_n);

            // Update using Numerical Scheme
            u_n_plus1[j] = dt*u_xxx;
        }

        //Boundary conditions      
        // long double dU = u_n_plus1[first_cell] - u_n_plus1[first_cell+1];
        // for(int i=first_cell-1; i>=0; i--)
        //     u_n_plus1[i] = u_n_plus1[i+1] - dU; 
        // u_n_plus1[Nx-1] = u_n_plus1[Nx-2] + dU;
        u_n_plus1[Nx-1] = 0; // RIGHT Boundary 
        u_n_plus1[0] = 0; // LEFT Boundary
        u_n_plus1[1] = u_n_plus1[0]; u_n_plus1[2] = u_n_plus1[1]; // Cells -1, -2
        
        // Store u^n+1_j in u^n_j for next time step
        u_n = u_n_plus1; 

        write_data(out_file, u_n, first_cell, last_cell); // Write Simulation Data for THIS time step
        //write_data(debug_file, u_n, 0, Nx-1);

        t+=dt;

    } write_data(fin_file, u_n, 1, Nx-2); // Write Simulation Data for FINAL time step

out_file.close(); fin_file.close();

}


// Calculates Third Derivative using Newton Difference Interpolation Derivative
vector<long double> third_derivative(vector<long double> &u){
    vector<long double> u_der(u.size(), 0);
    // First Derivative 
    for (int i = first_cell; i <= last_cell; i++)
        u_der[i] = derivative(u, i);
    u_der[Nx - 1] = 1;

    // // Second Derivative 
    // for (int i = first_cell; i <= last_cell; i++)
    //     u_der[i] = derivative(u_der, i);
    // u_der[Nx - 1] = 0;

    // // Third Derivative 
    // for (int i = first_cell; i = last_cell; i++)
    //     u_der[i] = derivative(u_der, i);
    // u_der[Nx - 1] = -1;

    return u_der;
}


/// @brief Newton Difference Interpolation Method for Derivative
/// @param u Spatial Array u^n_j
/// @return First Derivative of u^n_{x_index}
long double derivative(vector<long double> &u, int x_index) {
    long double u_x;
    
    // Calculate forward differnce del^n_u
    if(x_index <= 50) return (u[x_index+1] - u[x_index]) / dx;

    // Calculate backward differnce del^n_u
    else if(x_index >= 150) return (u[x_index] - u[x_index-1]) / dx; 
    
    // Calculate Stirling Interpoaltion
    else if(x_index > 50 && x_index < 150) return ((u[x_index+1] - u[x_index-1])) / (2*dx);

}


/// @brief Calculates the TV bound
long double calculate_tv(vector<long double> u) {
    double tv_norm = 0.0;
    for (int j = first_cell; j <= last_cell + 1; j++) {
        tv_norm += abs(u[j] - u[j - 1]);
    }
    
    return tv_norm;
}


///@brief Calcuates exact Solution
///@param condition Initial Condition
void uex(int condition){
    ofstream ex_file("dump_files/uex.txt"); // Exact Solution data dump file
    switch (condition){
        case 3:
            // Writing exact solution to file
            for(int j=first_cell; j<=last_cell; j++)
                ex_file << 0.25 * (pow(1.0/cosh(sqrt(0.5)/2.0 * (xmin + (j+0.5)*dx - 2.5) - 7.0), 2)) << " ";
            break;

        default: break;
    } ex_file.close();
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