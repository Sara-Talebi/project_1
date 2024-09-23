import numpy as np
import math
import matplotlib.pyplot as plt

##-----------------------------------------------------------------##
def rk4_method(n,t,dt,x,v,amp,omega):

    return 
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
#This part of the code is added for the testing purpose independent from main code.
def test_function_1B():
    n  = int(1e4)
    t0 = 0.0
    dt = 0.1
    x0 = 0.0
    v0 = 1.0
    amp   = 0.5
    omega = 5.0
    
    #x_list, v_list, t_list = euler_calculator(n,t0,dt,x0,v0,amp,omega)
    
    plt.figure()
    plt.plot(t_list,x_list)
    plt.show()
    
##-----------------------------------------------------------------##
test_function_1B()
##-----------------------------------------------------------------##
 
