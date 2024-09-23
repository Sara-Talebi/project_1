import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import solve_ivp,trapezoid,simpson
from integ_method import int_calculator
from ode_method import  hill_solv_scipy, hill_eq_solv


##-----------------------------------------------------------------##
#main function to define all variables and constant and calling functions. To observe any changes for ODE and integration just change parameters in this function. 
def main_function():
    #Parameters for the Hill eq. ODE
    n  = int(1e4)	#number of steps
    t0 = 0.0		#initial time
    dt = 0.1		#time step size
    x0 = 0.0		#initial displacement
    v0 = 1.0		#initial velocity
    amp   = 0.5		#amplitude (A)
    omega = 5.0		#frequency 
    
    #solve using Euler's method and Runge-Kutta 4th order method and returning a plot of both method results
    x_rk4, v_rk4, t_rk4 = hill_eq_solv(n,t0,dt,x0,v0,amp,omega)
    #solve Hill differential eq using scipy 
    sol = hill_solv_scipy(amp,omega)   
    
    #ploting the result of scipy/sara based RK4 for comarison
    plt.figure()
    plt.plot(sol.t,sol.y[0], label="scipy rk4")
    plt.plot(t_rk4,x_rk4, label='sara rk4')
    plt.xlabel('Time')
    plt.ylabel('Displacement')
    plt.title('Solution to the Hill eq')
    plt.legend()
    plt.show() 	
  
    
    #Set up integration limits and parameters
    up_lim  = 1e3	#upper limit 
    low_lim = -up_lim	#lower_limit
    n       = int(1e6)	#number of divisions for integration
    
    int_calculator(n,up_lim,low_lim)
    
    
##-----------------------------------------------------------------##
main_function()
##-----------------------------------------------------------------##
 
