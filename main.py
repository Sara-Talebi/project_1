import numpy as np
import matplotlib.pyplot as plt
import math
#import scipy as sp
import integ_method
from ode_method import euler_calculator, rk4_calculator

##-----------------------------------------------------------------##
#Define the ODE problem (Hill Differential Eqation)
def hill_eq_solv(n,t0,dt,x0,v0,amp,omega):    
    #Implementation of Euler's method and Runge-Kutta 4th order method
    x_em, v_em, t_em = euler_calculator(n,t0,dt,x0,v0,amp,omega)
    x_rk4, v_rk4, t_rk4 = rk4_calculator(n,t0,dt,x0,v0,amp,omega)
    
    return x_em, v_em, t_em, x_rk4, v_rk4, t_rk4
##-----------------------------------------------------------------##
def func(mass,a):
    f = m*a
    return f 
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
    
    #solve using Euler's method and Runge-Kutta 4th order method
    x_em, v_em, t_em, x_rk4, v_rk4, t_rk4 = hill_eq_solv(n,t0,dt,x0,v0,amp,omega)
    
     	
 	
    plt.figure()
    plt.plot(t_em,x_em, label="Euler's method")
    plt.plot(t_rk4,x_rk4, label='Runge-Kutta')
    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.title('ODE(Hill eq)')
    plt.legend()
    plt.show()
##-----------------------------------------------------------------##
main_function()
##-----------------------------------------------------------------##
 
