import numpy as np
import math

##-----------------------------------------------------------------##
def trapz_rule(n,up_lim,low_lim,func):
    delta_x = (up_lim - low_lim)/n
    k_sum   = 0.0
    for i  in range(n-1):
        temp_val = func(low_lim + delta_x)
        k_sum   += temp_val
        
    temp_val2 = (func(up_lim) + func(low_lim))/2
    final_val = delta_x * (temp_val2 + temp_val)
    return final_val
##-----------------------------------------------------------------##
def test_function():
    def func(x):
        a = math.exp(-(x**2))
        return a
    
    up_lim  = 1e3
    low_lim = 0
    n       = 1000
    
    ans = trapz_rule(n,up_lim,low_lim,func)
    print(ans)
    print(math.sqrt(math.pi)/2)
 	
 	
##-----------------------------------------------------------------##
test_function()
##-----------------------------------------------------------------##
 
