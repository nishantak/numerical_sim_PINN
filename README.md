# Numerical Simulation
Finite Volume Numerical analysis scheme <br>

# Steps to run
Make sure you have the *dependencies* installed and follow the steps for the code which you want to run / problem which you want to simulate.
## Dependencies: -
1. Python
2. C++ 
2. matplotlib
3. NumPy
4. SciPy
```bash
pip install matplotlib numpy scipy
```

## Transport_Burgers CPP Code: -
***FV_sim.h*** and ***FV_sim.cpp*** constitute A CUSTOM NUMERICAL ANALYSIS LIBRARY that contains _all the numerical analysis functions and scheme implementations._   
  1. Clone / Download this git repository into your _working directory_; <br><br> OR <br><br> Download the _**simulate.cpp, FV_sim.h, FV_sim.cpp, and matplot.py**_ files into a (same) directory;
    
  2. Run the **C++** file

```bash
cd [YOUR_WORKING_DIRECTORY]
g++ simulate.cpp FV_sim.cpp -o simulate && ./simulate
```

## Kuramoto Python Code: -
***scheme.py*** is the module that contains _the numerical analysis scheme implementation._ **functions.py*** is the module that contains functionality functions like *plot, write_data, and get_param.* 
<br> **config.py** contains all simulation parameters and the flux definition (f(u)). Changes can be made there as per need.
  1. Clone / Download this git repository into your _working directory_ ; <br><br> OR <br><br> Download the _**simulate.cpp, FV_sim.h, FV_sim.cpp, and matplot.py**_ files into a (same) directory;
    
  2. Run the *kuramoto_main.py* script
```bash
cd [YOUR_WORKING_DIRECTORY]
python kuramoto_main.py
```