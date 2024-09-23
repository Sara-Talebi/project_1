import numpy as np
import matplotlib.pyplot as plt
import math
from integ_method import int_calculator
from ode_method import  ode_solver


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
    
    #solve Hill eq using Euler's method and Runge-Kutta 4th order method and returning a plot of both developed bethod
    #solve Hill eq using scipy package and ploting the results + RK4 of mine developement for comparison 
    ode_solver(n,t0,dt,x0,v0,amp,omega)
    
  
    
    #Set up integration limits and parameters
    up_lim  = 1e3	#upper limit 
    low_lim = -up_lim	#lower_limit
    n       = int(1e6)	#number of divisions for integration
    
    #This function prints the result of gaussian integration based on trapz, simpson, riemann, scipy methods, and analytically
    int_calculator(n,up_lim,low_lim)
    
    
##-----------------------------------------------------------------##
main_function()
##-----------------------------------------------------------------##
 
