import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


##-----------------------------------------------------------------##
def hill_func(t, y, amp, omega):
    """
    Defining Hill eg, as a system of first-order ODEs to calculate by Scipy library
    """
    x, v = y  # A vector of [x,v]
    dxdt = v
    dvdt = -amp * np.sin(omega * t) * x
    return [dxdt, dvdt]

##-----------------------------------------------------------------##
def compute_energy(x, v):
    """
    Computes the total energy (kinetic + potential) for validation of energy conservation.
    Assuming mass m = 1 and spring constant k = 1 for simplicity.
    """
    m = 1.0  # Mass
    k = 1.0  # Spring constant
    
    kinetic_energy = 0.5 * m * v**2
    potential_energy = 0.5 * k * x**2
    total_energy = kinetic_energy + potential_energy
    
    return total_energy
##-----------------------------------------------------------------##
def rk4_calculator(n, t, dt, x, v, amp, omega):
    """
    Using the rk4_method function to loop over n_step for calculating the next n points
    """
    x_list = [x]  # creating a list to store x values
    v_list = [v]  # creating a list to store v values
    t_list = [t]  # creating a list to store t values

    for i in range(n):
        x_last = x_list[-1]
        v_last = v_list[-1]
        t_last = t_list[-1]
        new_x, new_v = rk4_method(n, t_last, dt, x_last, v_last, amp, omega)
        x_list.append(new_x)
        v_list.append(new_v)
        new_t = t_last + dt
        t_list.append(new_t)

    return x_list, v_list, t_list


##-----------------------------------------------------------------##
def rk4_method(n, t, dt, x, v, amp, omega):
    # k1=f(tn,yn)
    k1x = v
    k1v = -amp * np.sin(omega * t) * x

    # k2=f(tn+dt/2, yn+dt.k1/2)
    k2x = v + (dt / 2 * k1v)
    k2v = -amp * np.sin(omega * (t + (dt / 2))) * (x + (dt / 2 * k1x))

    # k3=ftn+dt/2, yn+dt.k2/2)
    k3x = v + (dt / 2 * k2v)
    k3v = -amp * np.sin(omega * (t + (dt / 2))) * (x + (dt / 2 * k2x))

    # k4=ftn+dt/2, yn+dt.k3/2)
    k4x = v + (dt * k3v)
    k4v = -amp * np.sin(omega * (t + (dt / 2))) * (x + (dt * k3x))

    # yn+1 = yn + dt/6(k1+2k2+2k3+k4)
    new_x = x + (dt / 6) * (k1x + (2 * k2x) + (2 * k3x) + k4x)
    new_v = v + (dt / 6) * (k1v + (2 * k2v) + (2 * k3v) + k4v)

    return new_x, new_v


##-----------------------------------------------------------------##
def euler_calculator(n, t, dt, x, v, amp, omega):
    """
    using the euler_method function to loop over n for calculating the next n points
    """
    x_list = [x]  # creating a list to store x values
    v_list = [v]  # creating a list to store v values
    t_list = [t]  # creating a list to store t values

    for i in range(n):
        x_last = x_list[-1]
        v_last = v_list[-1]
        t_last = t_list[-1]
        new_x, new_v = euler_method(n, t_last, dt, x_last, v_last, amp, omega)
        x_list.append(new_x)
        v_list.append(new_v)
        new_t = t_last + dt
        t_list.append(new_t)

    return x_list, v_list, t_list


##-----------------------------------------------------------------##
def euler_method(n, t, dt, x, v, amp, omega):
    new_x = x + (dt * v)
    new_v = v + (dt * -amp * np.sin(omega * t) * x)
    return new_x, new_v


##-----------------------------------------------------------------##
def ode_solver(n, t0, dt, x0, v0, amp, omega):
    """
    Implementing Euler's method and Runge-Kutta's 4th order method on the Hill differential equation and plotting them. also, plotting the RK4 result versus scipy provided RK45.
    """
    # Implementation of Euler's method and Runge-Kutta 4th order method
    x_em, v_em, t_em = euler_calculator(n, t0, dt, x0, v0, amp, omega)
    x_rk4, v_rk4, t_rk4 = rk4_calculator(n, t0, dt, x0, v0, amp, omega)

    # solve Hill differential eq using Scipy
    # parameters
    t_span = [t0, n*dt]  # time range from 0 to 1000
    y0 = [x0, v0]
    sol = solve_ivp(hill_func, t_span, y0, args=(amp, omega), method="RK45", t_eval=np.linspace(t0, n*dt, n+1))
    # Scipy solution
    t_scipy = sol.t
    x_scipy = sol.y[0]
    v_scipy = sol.y[1]
    
    # Energy Conservation Check for RK4 Method
    energy_rk4 = [compute_energy(x, v) for x, v in zip(x_rk4, v_rk4)]
    energy_scipy = [compute_energy(x, v) for x, v in zip(x_scipy, v_scipy)]
    # Error computation
    error_em = np.abs(np.array(x_em) - np.array(x_scipy))
    error_rk4 = np.abs(np.array(x_rk4) - np.array(x_scipy))
    
    # Plot solutions
    plt.figure()
    plt.plot(t_em, x_em, label="Euler's Method")
    plt.plot(t_rk4, x_rk4, label="Runge-Kutta (RK4)")
    plt.plot(t_scipy, x_scipy, label="Scipy RK45")
    plt.xlabel("Time")
    plt.ylabel("Displacement")
    plt.title("Solution to Hill's Equation (Euler, RK4, Scipy RK45)")
    plt.legend()
    plt.show()

    # Plot error
    plt.figure()
    plt.plot(t_em, error_em, label="Euler Error")
    plt.plot(t_rk4, error_rk4, label="RK4 Error")
    plt.xlabel("Time")
    plt.ylabel("Error (|Numerical - Scipy|)")
    plt.title("Error in RK4 vs Scipy RK45")
    plt.legend()
    plt.show()

    # Plot total energy over time for RK4 method to check energy conservation
    plt.figure()
    plt.plot(t_scipy, energy_scipy, label="Total Energy (scipy)")
    plt.plot(t_rk4, energy_rk4, label="Total Energy (RK4)")
    plt.xlabel("Time")
    plt.ylabel("Total Energy")
    plt.title("Energy Conservation in RK4 and Scipy RK45 Method")
    plt.legend()
    plt.show()
    
##-----------------------------------------------------------------##
