# Numerical Simulation
Finite Volume Numerical analysis scheme <br><br>

# Steps to run Transport_Burgers CPP Code: -
***numSim.h*** and ***numSim.cpp*** constitute A CUSTOM NUMERICAL ANALYSIS LIBRARY that contains _all the numerical analysis functions and scheme implementations._
  1. Make sure you have **Python**, **matplotlib** and **C/C++** installed;
     
  2. Clone / Download this git repository into your _working directory_ (Make sure the Python script and C++ file are in the **same** directory) ; <br><br> OR <br><br> Download the _**simulate.cpp, FV_sim.h, FV_sim.cpp, and matplot.py**_ files into a (same) directory;
    
  3. Run the **C++** file (Make sure that if you change simulation parameters then the matplot parameters are tweaked accordingly)

```bash
cd [YOUR_WORKING_DIRECTORY]
g++ simulate.cpp FV_sim.cpp -o simulate && ./simulate

# Steps to run Kuramoto Python Code: -