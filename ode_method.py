import numpy as np
import math
import matplotlib.pyplot as plt


##-----------------------------------------------------------------##
def rk4_calculator(n,t,dt,x,v,amp,omega):
    x_list = [x]
    v_list = [v]
    t_list = [t]

    	
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
def euler_calculator(n,t,dt,x,v,amp,omega):
    x_list = [x]
    v_list = [v]
    t_list = [t]

    	
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
#This part of the code is added for the testing purpose independent from main code. make sure to uncomment #test_function_1B() down below.
def test_function_1B():
    n  = int(1e4)
    t0 = 0.0
    dt = 0.1
    x0 = 0.0
    v0 = 1.0
    amp   = 0.5
    omega = 5.0
    
    x_list1, v_list1, t_list1 = euler_calculator(n,t0,dt,x0,v0,amp,omega)
    x_list2, v_list2, t_list2 = rk4_calculator(n,t0,dt,x0,v0,amp,omega)
    
    plt.figure()
    plt.plot(t_list1,x_list1, label='Euler')
    plt.plot(t_list2,x_list2, label='RK4')
    plt.legend()
    plt.show()
    
##-----------------------------------------------------------------##
#test_function_1B()
##-----------------------------------------------------------------##
 
