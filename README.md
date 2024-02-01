# Numerical Simulation
Finite Volume Numerical analysis scheme <br>

# Steps to run
Make sure you have the *dependencies* installed and follow the steps for the code which you want to run / problem which you want to simulate. 
<br><br> (All shell commands are according to a windows environment)
## Dependencies: -
1. Python (and add to PATH)
2. C++ (and add to PATH)
2. matplotlib
3. NumPy
4. SciPy
```bash
pip install matplotlib numpy scipy
```

## Transport_Burgers CPP Code: -
***FV_sim.h*** and ***FV_sim.cpp*** constitute A CUSTOM NUMERICAL ANALYSIS LIBRARY that contains _all the numerical analysis functions and scheme implementations._   
  1. Clone / Download this git repository into your _working directory_; <br><br> OR <br><br> Download the _**simulate.cpp, FV_sim.h, FV_sim.cpp, and matplot.py**_  files into a (same) directory;
    
  2. Run the **simulate.cpp** file

```bash
cd [YOUR_WORKING_DIRECTORY]
g++ simulate.cpp FV_sim.cpp -o simulate && ./simulate
```

## Kuramoto Python Code: -
***scheme.py*** is the module that contains _the numerical analysis scheme implementation._ ***functions.py*** is the module that contains functionality functions *(plot, write_data, and get_param)* 
<br><br> **config.py** contains all simulation parameters and the flux definition (f(u)). Changes can be made there as per simulation need.
  1. Clone / Download this git repository into your _working directory_ ; <br><br> OR <br><br> Download the _**scheme.py, functions.py, config.py, and kuramoto_main.py**_  files into a (same) directory;
    
  2. Run the *kuramoto_main.py* script
```bash
cd [YOUR_WORKING_DIRECTORY]
python kuramoto_main.py
```