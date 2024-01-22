#include<bits/stdc++.h>
using namespace std;

// Simulation Parameters 
double xmin = 0, xmax = 2*M_PI;  // Domain limits
double L = abs(xmax- xmin);   //Domain Length
int ghost_cells = 1;    // Number of ghost cells 
int Nx = 200 + ghost_cells;   // Number of spatial points
long double dx = L/(Nx-1);  // Cell width 

int first_cell = 0, last_cell = Nx-2;   // j domain Limits

/// @brief Newton Difference Interpolation Method for Derivative
/// @param u Spatial Array u^n_j
/// @return First Derivative of u^n_{x_index}
long double derivative(vector<long double> &u, int x_index, int count, long double u_x) {
    // Base case when del^n_y = 0 or only 1 differnce 
    if (u.size()<=1 || u[1]-u[0] == 0){
        u_x /= (dx);
        return u_x; 
    }

    vector<long double> del_u(u.size() - x_index- 1, 0);
    if(x_index <= 50){
        // Calculate forward differnce del^n_u
        for (int i = x_index; i < u.size() - 1; i++) {
            del_u[i-x_index] = u[i + 1] - u[i];

        } u_x += pow(-1, count++) * (1.0/count) * del_u[0];
            
        return derivative(del_u, 0, count, u_x); // Recursive call to caculate forward difference of del^n_u
    }

    else if(x_index >= 150){
        // Calculate backward differnce del^n_u
        for (int i = x_index; i < 0; i++) {
            del_u[i-x_index] = u[i] - u[i-1];

        } u_x += pow(-1, count++) * (1.0/count) * del_u[0];
            
        return derivative(del_u, 0, count, u_x); // Recursive call to caculate forward difference of del^n_u
    }
    
    else if(x_index > 50 && x_index < 150){
        //Implement Stirling Interpoaltion
    }

}


// Calculates Third Derivative using Newton Difference Interpolation Derivative
vector<long double> third_derivative(vector<long double> &u){
    
    vector<long double> u_der(u.size(), 0);

    // First Derivative 
    for (int i = first_cell; i < last_cell; i++)
        u_der[i] = derivative(u, i, 0, 0);
    u_der[Nx - 1] = 1;

    // Second Derivative 
    for (int i = first_cell; i < last_cell; i++)
        u_der[i] = derivative(u_der, i, 0, 0);
    u_der[Nx - 1] = 0;

    // Third Derivative 
    for (int i = first_cell; i < last_cell; i++)
        u_der[i] = derivative(u_der, i, 0, 0);
    u_der[Nx - 1] = -1;

    return u_der;
}


int main(){
    vector<long double>u={0.015707317, 0.047106451, 0.078459096, 0.10973431, 0.14090123, 0.1719291, 0.2027873, 0.23344536, 0.26387305, 0.29404033, 0.32391742, 0.35347484, 0.38268343, 0.41151436, 0.43993917, 0.46792981, 0.49545867, 0.52249856, 0.54902282, 0.57500525, 0.60042023, 0.62524266, 0.64944805, 0.67301251, 0.6959128, 0.7181263, 0.73963109, 0.76040597, 0.78043041, 0.79968466, 0.81814972, 0.83580736, 0.85264016, 0.86863151, 0.88376563, 0.89802758, 0.91140328, 0.92387953, 0.93544403, 0.94608536, 0.95579301, 0.96455742, 0.97236992, 0.97922281, 0.98510933, 0.99002366, 0.99396096, 0.99691733, 0.99888987, 0.99987663, 0.99987663, 0.99888987, 0.99691733, 0.99396096, 0.99002366, 0.98510933, 0.97922281, 0.97236992, 0.96455742, 0.95579301, 0.94608536, 0.93544403, 0.92387953, 0.91140328, 0.89802758, 0.88376563, 0.86863151, 0.85264016, 0.83580736, 0.81814972, 0.79968466, 0.78043041, 0.76040597, 0.73963109, 0.7181263, 0.6959128, 0.67301251, 0.64944805, 0.62524266, 0.60042023, 0.57500525, 0.54902282, 0.52249856, 0.49545867, 0.46792981, 0.43993917, 0.41151436, 0.38268343, 0.35347484, 0.32391742, 0.29404033, 0.26387305, 0.23344536, 0.2027873, 0.1719291, 0.14090123, 0.10973431, 0.078459096, 0.047106451, 0.015707317, -0.015707317, -0.047106451, -0.078459096, -0.10973431, -0.14090123, -0.1719291, -0.2027873, -0.23344536, -0.26387305, -0.29404033, -0.32391742, -0.35347484, -0.38268343, -0.41151436, -0.43993917, -0.46792981, -0.49545867, -0.52249856, -0.54902282, -0.57500525, -0.60042023, -0.62524266, -0.64944805, -0.67301251, -0.6959128, -0.7181263, -0.73963109, -0.76040597, -0.78043041, -0.79968466, -0.81814972, -0.83580736, -0.85264016, -0.86863151, -0.88376563, -0.89802758, -0.91140328, -0.92387953, -0.93544403, -0.94608536, -0.95579301, -0.96455742, -0.97236992, -0.97922281, -0.98510933, -0.99002366, -0.99396096, -0.99691733, -0.99888987, -0.99987663, -0.99987663, -0.99888987, -0.99691733, -0.99396096, -0.99002366, -0.98510933, -0.97922281, -0.97236992, -0.96455742, -0.95579301, -0.94608536, -0.93544403, -0.92387953, -0.91140328, -0.89802758, -0.88376563, -0.86863151, -0.85264016, -0.83580736, -0.81814972, -0.79968466, -0.78043041, -0.76040597, -0.73963109, -0.7181263, -0.6959128, -0.67301251, -0.64944805, -0.62524266, -0.60042023, -0.57500525, -0.54902282, -0.52249856, -0.49545867, -0.46792981, -0.43993917, -0.41151436, -0.38268343, -0.35347484, -0.32391742, -0.29404033, -0.26387305, -0.23344536, -0.2027873, -0.1719291, -0.14090123, -0.10973431, -0.078459096, -0.047106451, -0.015707317};
    
    vector<long double> u_xxx(Nx,0);
    u_xxx = third_derivative(u);

    return 0;
}