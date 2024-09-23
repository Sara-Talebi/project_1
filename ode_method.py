import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

##-----------------------------------------------------------------##
#solve Hill differential eq using scipy
def hill_solv_scipy(amp,omega):
    #parameters
    t_span = [0.0, 1000.0]	#time range from 0 to 1000
    y0     = [0.0,1.0] 
    sol = solve_ivp(hill_func,t_span,y0,args=(amp,omega),method='RK45')
    return sol
##-----------------------------------------------------------------##
#Defining Hill eg as a system of first-order ODEs to calculate by scipy
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
#using the rk4_method function to loop over n for calculating next points
def rk4_calculator(n,t,dt,x,v,amp,omega):
    x_list = [x]	#creating a list to store x values
    v_list = [v]	#creating a list to store v values
    t_list = [t]	#creating a list to store t values

    	
    for i in range(n):
        x_last = x_list[-1]
        v_last = v_list[-1]
        t_last = t_list[-1]
        new_x, new_v = rk4_method(n,t_last,dt,x_last,v_last,amp,omega)
        x_list.append(new_x)
        v_list.append(new_v)
        new_t = t_last + dt
        t_list.append(new_t)
  
    return x_list, v_list, t_list
##-----------------------------------------------------------------##
def rk4_method(n,t,dt,x,v,amp,omega):
    #k1=f(tn,yn)
    k1x = v
    k1v = -amp * np.sin(omega*t) * x
    
    #k2=f(tn+dt/2, yn+dt.k1/2)
    k2x = v + (dt/2 * k1v)
    k2v = -amp * np.sin(omega*(t + (dt/2))) * (x + (dt/2 * k1x))
    
    #k3=ftn+dt/2, yn+dt.k2/2)
    k3x = v + (dt/2 * k2v)
    k3v = -amp * np.sin(omega*(t + (dt/2))) * (x + (dt/2 * k2x))
    
    #k4=ftn+dt/2, yn+dt.k3/2)
    k4x = v + (dt * k3v)
    k4v = -amp * np.sin(omega*(t + (dt/2))) * (x + (dt * k3x))
    
    #yn+1 = yn + dt/6(k1+2k2+2k3+k4)
    new_x = x + (dt/6)*(k1x + (2 * k2x) + (2 * k3x) + k4x)
    new_v = v + (dt/6)*(k1v + (2 * k2v) + (2 * k3v) + k4v)
    
    return new_x, new_v
##-----------------------------------------------------------------##
#using the euker_method function to loop over n for calculating next points
def euler_calculator(n,t,dt,x,v,amp,omega):
    x_list = [x]	#creating a list to store x values
    v_list = [v]	#creating a list to store v values
    t_list = [t]	#creating a list to store t values

    	
    for i in range(n):
        x_last = x_list[-1]
        v_last = v_list[-1]
        t_last = t_list[-1]
        new_x, new_v = euler_method(n,t_last,dt,x_last,v_last,amp,omega)
        x_list.append(new_x)
        v_list.append(new_v)
        new_t = t_last + dt
        t_list.append(new_t)
  
    return x_list, v_list, t_list
##-----------------------------------------------------------------##
def euler_method(n,t,dt,x,v,amp,omega):
    new_x = x + (dt * v)
    new_v = v + (dt * -amp * np.sin(omega*t) * x)
    return new_x, new_v
##-----------------------------------------------------------------##
 
