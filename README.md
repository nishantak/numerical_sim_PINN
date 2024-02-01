# Numerical Simulation
Finite Volume Numerical analysis scheme <br>

# Steps to run
<br>
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
  1. Clone / Download this git repository into your _working directory_ (Make sure the Python script and C++ file are in the **same** directory) ; <br><br> OR <br><br> Download the _**simulate.cpp, FV_sim.h, FV_sim.cpp, and matplot.py**_ files into a (same) directory;
    
  2. Run the **C++** file (Make sure that if you change simulation parameters then the matplot parameters are tweaked accordingly)

```bash
cd [YOUR_WORKING_DIRECTORY]
g++ simulate.cpp FV_sim.cpp -o simulate && ./simulate
```

## Kuramoto Python Code: -

