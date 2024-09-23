import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import solve_ivp,trapezoid,simpson
from integ_method import trapz_rule_v1,composite_simpson_v2,riemann_left,riemann_right
from ode_method import euler_calculator, rk4_calculator

##-----------------------------------------------------------------##
#solve Hill differential eq using scipy
def hill_solv_scipy(amp,omega):
    #parameters
    t_span = [0.0, 1000.0]	#time range from 0 to 1000
    y0     = [0.0,1.0] 
    sol = solve_ivp(hill_func,t_span,y0,args=(amp,omega),method='RK45')
    return sol
##-----------------------------------------------------------------##
#Defining Hill eg as a system of first-order ODEs
def hill_func(t,y,amp,omega):
    x, v = y 	#A vector of [x,v]
    dxdt = v
    dvdt = -amp * np.sin(omega*t) * x
    return [dxdt, dvdt]
##-----------------------------------------------------------------##
#Define the ODE problem (Hill Differential Eqation) as a system of first-order ODEs
def hill_eq_solv(n,t0,dt,x0,v0,amp,omega):    
    #Implementation of Euler's method and Runge-Kutta 4th order method
    x_em, v_em, t_em = euler_calculator(n,t0,dt,x0,v0,amp,omega)
    x_rk4, v_rk4, t_rk4 = rk4_calculator(n,t0,dt,x0,v0,amp,omega)
    
    #Plot to compare the Hill equation's solutions obtained by Euler's method and RK4	
    plt.figure()
    plt.plot(t_em,x_em, label="Euler's method")
    plt.plot(t_rk4,x_rk4, label='Runge-Kutta method')
    plt.xlabel('Time')
    plt.ylabel('Displacement')
    plt.title('Solution to the Hill eq')
    plt.legend()
    plt.show()
    
    #Returning RK4 calculated solutions 
    return x_rk4, v_rk4, t_rk4
##-----------------------------------------------------------------##
#Define the gaussian function to integrate 
def func(x):
    a = np.exp(-(x**2))
    return a
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
    x_rk4, v_rk4, t_rk4 = hill_eq_solv(n,t0,dt,x0,v0,amp,omega)
    #solve Hill differential eq using scipy
    sol = hill_solv_scipy(amp,omega)   
    
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
    
    print("solution to the gaussian function by sara-developed methods :D")
    #Compute the integral using Trapezoidal rule
    ans1 = trapz_rule_v1(n,up_lim,low_lim,func)
    print("Trapezoidal= \t", ans1)
    
    #Compute the integral using Composite-Simpson's rule
    val2 = composite_simpson_v2(n,up_lim,low_lim,func)
    print("Simpson= \t", val2)
    
    #Compute the integral using Riemann left and right rule
    val1 = riemann_left(n,up_lim,low_lim,func)
    val2 = riemann_right(n,up_lim,low_lim,func)
    print("Riemann_left= \t", val1)
    print("Riemann_right= \t", val2)
    
    #Integrating the Gaussian function using scipy
    print("solution to the gaussian function using scipy")
    samp_point = np.linspace(low_lim, up_lim, n)
    ans = trapezoid(func(samp_point),samp_point)
    print("Trapezoidal= \t", ans)
    val = simpson(func(samp_point),samp_point)
    print("Simpson= \t", val)
    
    #The analytical answer to the qaussian function
    print("The analytical answer to the qaussian function")
    print("Analytical= \t",math.sqrt(math.pi))
    
    
    
##-----------------------------------------------------------------##
main_function()
##-----------------------------------------------------------------##
 
