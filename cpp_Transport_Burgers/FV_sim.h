#ifndef FV_SIM_H
#define FV_SIM_H

#include <vector>
#include <fstream>

using namespace std;

// Global Variables | Simulation Parameters
extern int equation; // Problem Statement 

extern double xmin, xmax, L;    // Spatial Domain Parameters
extern int Nx, ghost_cells;     // Spatail Domain Points
extern long double dx;   // Cell Width

extern double cfl;  // Stability Parameter - CFL Number
extern long double c;   // Wave Velocity

extern long double dt;  // Time Step
extern double Tf;   // Final time / Total Time
extern int Nt;  // No. of time steps

extern int first_cell, last_cell;   // j domain Limits

extern float a; // Constant multiplier of flux


// Function Prototypes
long double flux(long double u);
long double num_flux(long double u, long double u_next, int scheme);
void initialise(vector<long double> &u, int condition);
void simulate(vector<long double> &u_n, int flux_scheme, int boundary_condition);
void u_ex(int condition);
void get_param();
void plot();
void write_data(ofstream &filename, vector<long double> u, int first_cell, int last_cell);
long double calculate_tv(vector<long double> u);
//Unused functions
vector<long double> third_derivative(vector<long double> &);
long double derivative(vector<long double>&, int);


#endif 
