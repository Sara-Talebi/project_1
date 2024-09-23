import numpy as np
import math

##-----------------------------------------------------------------##
def euler_method(n,t,dt,x,g,amp,omega):
    x_list = []
    g_list = []
    
    x_list.append(x0)
    g_list.append(g0)
    	
    for i in range(n):
        x_last = x_list[-1]
        g_last = g_list[-1]
        new_x, new_g = euler_method(n,t,dt,x_last,g_last,amp,omega)
        x_list.append(new_x)
        g_list.append(new_g)
  
    return x_list, g_list
##-----------------------------------------------------------------##
def euler_method(n,t,dt,x,g,amp,omega):
    new_x = x + (dt * g)
    new_g = g + (dt * -amp * np.sin(omega*t) * x)
    return new_x, new_g
##-----------------------------------------------------------------##
def test_function():
    n  = int(1e2)
    t  = 0.0
    dt = 0.1
    x0 = 0.0
    g0 = 0.0
    amp   = 1.0
    omega = 1.0
    
    
##-----------------------------------------------------------------##
test_function()
##-----------------------------------------------------------------##
 
