import numpy as np
import math

##-----------------------------------------------------------------##
def composite_simpson(n,up_lim,low_lim,func):
    h = (up_lim-low_lim)/n
    k1 = func(low_lim)
    k2 = 0.0 
    b = int((n+1)/2)
    for i in range(1, b):
        k2 += func(low_lim + (2 * i * h) - h)
    k3 = 0.0
    b = int(((n+1)/2)-1)
    for j in range(1, b):
    	k2 += func(low_lim + (2 * i * h))
    k4 = func(up_lim)
    final_val = (h/3) * (k1 + (4 * k2) + (2 * k3) + k4)
    return final_val
##-----------------------------------------------------------------##
def simpson_rule(up_lim,low_lim,func):
    h = (up_lim-low_lim)/6
    final_val = h * (func(low_lim) + (4 * func((up_lim+low_lim)/2)) + func(up_lim))
    return final_val
##-----------------------------------------------------------------##
def trapz_rule_v1(n,up_lim,low_lim,func):
    delta_x = (up_lim - low_lim)/n
    k_sum   = 0.0
    for i  in range(n-1):
        temp_val = func(low_lim + (i * delta_x))
        k_sum   += temp_val
        
    temp_val2 = (func(up_lim) + func(low_lim))/2
    final_val = delta_x * (temp_val2 + temp_val)
    return final_val
##-----------------------------------------------------------------##
def trapz_rule_v2(n,up_lim,low_lim,func):
    vol = up_lim - low_lim
    
    delta_x = (up_lim - low_lim)/n
    f_sum   = 0.0
    for i  in range(n-1):
        temp_val = func(low_lim + (i * delta_x))
        f_sum   += temp_val
    temp_val2 = ((func(up_lim) + func(low_lim))/2) + f_sum
    f_avg = temp_val2/n
    final_val = vol * f_avg
    return final_val
##-----------------------------------------------------------------##
def test_function():
    def func(x):
        a = math.exp(-(x**2))
        return a
    
    up_lim  = 1e3
    low_lim = 0.0
    n       = 1000000
    
    ans = trapz_rule_v2(n,up_lim,low_lim,func)
    print("trapizoidal= \t", ans)
    
    val2 = composite_simpson(n,up_lim,low_lim,func)
    print("simpson= \t", val2)
    
    print("actual= \t", math.sqrt(math.pi)/2)
##-----------------------------------------------------------------##
test_function()
##-----------------------------------------------------------------##
 
