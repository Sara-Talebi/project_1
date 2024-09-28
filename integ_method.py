import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid, simpson


##-----------------------------------------------------------------##
def func(x):
    """
    Define the gaussian function to integrate using Trapezoidal, Simpson, and Riemann methods.
    """
    a = np.exp(-(x**2))
    return a
##-----------------------------------------------------------------##
def riemann_right(n, up_lim, low_lim, func):
    """
    For the left rule, the function is approximated by its values at the left endpoints of the subintervals. This gives multiple rectangles with base Δx and height f(a + iΔx). Doing this for i = 0, 1, ..., n − 1, and summing the resulting areas gives the approximate answer.
    """
    delta_x = (up_lim - low_lim) / n
    r_sum = 0.0
    for i in range(1, n + 1):
        temp_val = func(low_lim + (i * delta_x)) * delta_x
        r_sum += temp_val
    return r_sum, delta_x
##-----------------------------------------------------------------##
def riemann_left(n, up_lim, low_lim, func):
    """
    For the right rule, the function is approximated by its values at the right endpoints of the subintervals. This gives multiple rectangles with base Δx and height f(a + iΔx). Doing this for i = 1, ..., n, and summing the resulting areas gives the approximate answer.
    """
    delta_x = (up_lim - low_lim) / n
    r_sum = 0.0
    for i in range(n):
        temp_val = func(low_lim + (i * delta_x)) * delta_x
        r_sum += temp_val
    return r_sum, delta_x
##-----------------------------------------------------------------##
def composite_simpson_v1(n, up_lim, low_lim, func):
    h = (up_lim - low_lim) / n
    k1 = func(low_lim)
    k2 = 0.0
    b = int((n) / 2)
    for i in range(1, b):
        k2 += func(low_lim + (2 * i * h) - h)
    k3 = 0.0
    b = int((n / 2) - 1)
    for j in range(1, b):
        k2 += func(low_lim + (2 * i * h))
    k4 = func(up_lim)
    final_val = (h / 3) * (k1 + (4 * k2) + (2 * k3) + k4)
    return final_val, delta_x
##-----------------------------------------------------------------##
def composite_simpson_v2(n, up_lim, low_lim, func):
    delta_x = (up_lim - low_lim) / n
    s_sum = 0.0
    for i in range(n + 1):
        if i == 0:
            s_sum += func(low_lim)
        elif i % 2 == 0:
            s_sum += 2 * func(low_lim + (i * delta_x))
        elif i % 2 != 0:
            s_sum += 4 * func(low_lim + (i * delta_x))
        elif i == n:
            s_sum += func(up_lim)
        final_val = s_sum * (delta_x / 3)
    return final_val, delta_x
##-----------------------------------------------------------------##
def simpson_rule(up_lim, low_lim, func):
    h = (up_lim - low_lim) / 6
    final_val = h * (func(low_lim) + (4 * func((up_lim + low_lim) / 2)) + func(up_lim))
    return final_val
##-----------------------------------------------------------------##
def trapz_rule_v1(n, up_lim, low_lim, func):
    """
    Uniform grid trapezoidal rule, for a domain discretized into N equally spaced panels.
    """
    delta_x = (up_lim - low_lim) / n
    k_sum = 0.0
    for i in range(1, n):
        temp_val = func(low_lim + (i * delta_x))
        k_sum += temp_val
    temp_val2 = (func(up_lim) + func(low_lim)) / 2
    final_val = delta_x * (temp_val2 + k_sum)
    return final_val, delta_x
##-----------------------------------------------------------------##
def trapz_rule_v2(n, up_lim, low_lim, func):
    vol = up_lim - low_lim
    delta_x = (up_lim - low_lim) / n
    f_sum = 0.0
    for i in range(n):
        temp_val = func(low_lim + (i * delta_x))
        f_sum += temp_val
    temp_val2 = ((func(up_lim) + func(low_lim)) / 2) + f_sum
    f_avg = temp_val2 / n
    final_val = vol * f_avg
    return final_val, dleta_x
##-----------------------------------------------------------------##
def int_calculator(n, up_lim, low_lim):
    """
    Solving the Gaussian function using Trapezoidal, Simpson, Riemann method and also calculating by scipy library (trapz and simpson). At the last line, printing the exact analytical solution too.
    """
    
    analytical_solution = math.sqrt(math.pi)  # The exact analytical answer to the Gaussian function
    
    error_trap = []
    error_simp = []
    error_riemann_left = []
    error_riemann_right = []
    error_scipy_trap = []
    error_scipy_simp = []
    n_step = []
    
    dx_trapz = []
    dx_simp = []
    dx_riem_left = []
    
    
    #Increasing the number of steps to increase the acuracy
    
    while n >= 100:
        n = int(n)
        print("Results of number of steps equals to ", n)
        print("Solution to the Gaussian function by Sara-developed methods :D")
	    
        # Compute the integral using Trapezoidal rule
        ans1, dx = trapz_rule_v1(n, up_lim, low_lim, func)
        e_trap = np.abs(analytical_solution - ans1)
        error_trap.append(e_trap)
        dx_trapz.append(dx)
        print(f"Trapezoidal= \t {ans1} Error= {e_trap}")

        # Compute the integral using Composite-Simpson's rule
        val2, dx = composite_simpson_v2(n, up_lim, low_lim, func)
        e_simp = np.abs(analytical_solution - val2)
        error_simp.append(e_simp)
        dx_simp.append(dx)
        print(f"Simpson= \t {val2} Error= {e_simp}")

        # Compute the integral using Riemann left and right rule
        val1, dx = riemann_left(n, up_lim, low_lim, func)
        val2, dx = riemann_right(n, up_lim, low_lim, func)
        e_riemann_left = np.abs(analytical_solution - val1)
        error_riemann_left.append(e_riemann_left)
        dx_riem_left.append(dx)
        error_riemann_right.append(np.abs(analytical_solution - val2))
        print(f"Riemann_left= \t {val1} Error= {error_riemann_left}")
        print(f"Riemann_right= \t {val2} Error= {error_riemann_right}")

        # Integrating the Gaussian function using Scipy
        samp_point = np.linspace(low_lim, up_lim, n)
        ans = trapezoid(func(samp_point), samp_point)
        val = simpson(func(samp_point), samp_point)
        error_scipy_trap.append(np.abs(analytical_solution - ans))
        error_scipy_simp.append(np.abs(analytical_solution - val))
        print(f"Scipy Trapezoidal= \t {ans} Error= {error_scipy_trap}")
        print(f"Scipy Simpson= \t {val} Error= {error_scipy_simp}")
	
        n_step.append(n)
        n -= 100
       

        # Print the exact analytical solution for comparison
        #print(f"Analytical= \t {analytical_solution}")
    
    # Plot error
    plt.figure()
    plt.plot(n_step, error_trap, label="Trapezoidal Error")
    plt.plot(n_step, error_simp, label="Simpson Error")
    plt.plot(n_step, error_riemann_left, 'D', label="Riem_left Error")
    #plt.plot(n_step, error_riemann_right, 'D', label="Riem_right Error")
    plt.plot(n_step, error_scipy_trap, label="Scipy_trapz Error")
    plt.plot(n_step, error_scipy_simp, '*', label="Scipy_simp Error")
    #plt.gca().invert_xaxis()
    plt.xlabel("Number of steps (n)")
    plt.ylabel("Absolute Error (|Analytical - Numerical|)")
    plt.legend()
    plt.show()

    #print("n_step=", n_step)
    #print("error_trap=", error_trap)
    #print("error_simp=", error_simp)
    #print("error_riemann_left=", error_riemann_left)
    #print("error_riemann_right=", error_riemann_right)
    #print("error_scipy_trap=", error_scipy_trap)
    #print("error_scipy_simp=", error_scipy_simp)
##-----------------------------------------------------------------##
