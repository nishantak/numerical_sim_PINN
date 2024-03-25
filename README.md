# Numerical Simulation
Finite Volume Numerical analysis scheme <br>

# Steps to run
Ensure you have the *dependencies* installed and follow the steps for the problem which you want to simulate. 
<br><br> (All shell commands are for a windows environment)
### Dependencies: -
1. Python (in PATH)
2. C++ (in PATH)
2. matplotlib
3. NumPy
4. SciPy (for kuramoto)
```bash
pip install matplotlib numpy scipy
```
```bash
git clone 
```
## Transport_Burgers C++ Code: -
***FV_sim.h*** and ***FV_sim.cpp*** constitute A CUSTOM NUMERICAL ANALYSIS LIBRARY that contains _all the numerical analysis functions and scheme implementations._   
  1. Clone this git repository into your _working directory_;

  2. Run the ***simulate.cpp*** file

```bash
cd [YOUR_WORKING_DIRECTORY]
g++ simulate.cpp FV_sim.cpp -o simulate && ./simulate
```

## Kuramoto Python Code: -
In **each directory**, ***scheme.py*** is the module that contains _the respective numerical analysis scheme implementation._ ***functions.py*** is the module that contains functionality functions *(plot, write_data, and get_param)* 
<br><br> **config.py** contains all simulation parameters and the flux definition (f(u)). Changes can be made there as per simulation need.
  1. Clone / Download this git repository into your _working directory_ ;

  2. Run the ***kuramoto_main.py*** script
```bash
cd [YOUR_WORKING_DIRECTORY]
python kuramoto_main.py
```