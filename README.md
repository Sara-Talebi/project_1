# project_1
**Overview**
This project involves solving two physics problems using numerical methods:

Solving the Hill, second-order differential equation using Euler's and Runge-Kutta 4th order (RK4) methods.
Calculating the definite integral of a Gaussian function using various numerical integration methods.
Both problems are compared against exact solutions or Scipy-provided algorithms for validation.

**File Structure**

   **main.py:** This file contains the main function where variables and constants are defined, and functions are called for solving the ODE and the integral. To modify parameters for the ODE or the integration, adjust them directly in the script.
   
  **ode_method.py** Contains implementations of Euler's method and Runge-Kutta 4th order (RK4) method for solving the Hill differential equation. It also compares the results of these methods to Scipy's RK45 solver. Additionally, the code computes the total energy (kinetic + potential) to verify energy conservation in the system.
   
  **integ_method.py:** Implements the following numerical integration methods for solving the Gaussian function:
   
         * Riemann sum (left and right)
         * Composite Simpson’s rule (two versions)
         * Simpson's rule (1/3)
         * Trapezoidal rule (two versions)
         * Scipy’s integration methods (Trapezoidal and Simpson)
**Error Analysis:**
   For both the ODE and the integration problems, errors are computed and printed:
   
   __ODE Problem:__ Errors in the Euler and RK4 methods are compared to Scipy’s RK45 solution in different number of steps.
   
   __Integral Problem:__ Errors are computed by comparing the numerical results with the exact analytic solution of the Gaussian function.

**How to Use:**
  
   By running the *main.py*, the ode_solver() and int_calculator() functions are called within main.py. The plots are shown, and the results are printed.
