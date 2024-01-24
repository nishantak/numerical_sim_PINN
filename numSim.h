#ifndef NUMSIM_H
#define NUMSIM_H

#include <vector>
#include <fstream>

using namespace std;

// Global Variables | Simulation Parameters
extern double xmin, xmax, L;    // Spatial Domain Parameters
extern int Nx, ghost_cells;     // Spatail Domain Points
extern long double dx;   // Cell Width

extern double cfl;  // Stability Parameter - CFL Number
extern long double c;   // Wave Velocity

extern long double dt;  // Time Step
extern double Tf;   // Final time / Total Time
extern int Nt;  // No. of time steps

extern int first_cell, last_cell;   // j domain Limits

// Function Prototypes
long double flux(long double);
long double num_flux(long double, long double, int);
void intialise(vector<long double>&, int);
vector<long double> third_derivative(vector<long double> &);
void simulate(vector<long double>&, int);
void uex(int);
void get_param();
void plot();
void write_data(ofstream&, vector<long double>, int, int);
long double calculate_tv(vector<long double>);
long double derivative(vector<long double>&, int);
// Function Prototypes

#endif 
