from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


R = 2         # Amount of rainfall [mm/day]
r_m = 0.3     # Maximum growth rate of tree cover [1/year]
h_R = 0.5     # Half-saturation constant for rainfall effect on growth [mm/day]
m_n = 0.15    # Natural mortality rate of trees [1/year]
h_n = 0.1     # Half-saturation constant for natural mortality
m_f = 0.11    # Mortality rate of trees due to external factors [1/year]
h_f = 0.6     # Half-saturation constant for external mortality factors
p = 7         # Hill coefficient for external mortality response
k = 0.9       # Carrying capacity of tree cover
b = 2         # Baseline tree reproduction rate [mm/day]
r_R = 1       # Rainfall regeneration rate [1/year]




def new_model(mu, s_n, T0, R0, S0):
    """
    Simulates the coupled tree cover-rainfall dynamics over time.
    
    Returns:
    T_arr : Array containing tree cover levels over time.
    """
    # Constant rainfall rate for the simulation
    R_const = 2

    # Growth rate of tree cover as a function of rainfall
    r = lambda x: r_m * x / (h_R + x)

    # Differential equation for tree cover growth
    def dTdt(t, y1, y2, y3):
        return (r(y2) * y1 * (1 - y1 / k) - 
                m_n * y1 * h_n / (y1 + h_n) -
                m_f * y1 * (h_f ** p) / (y1 ** p + h_f ** p) -
                mu*y1*s_n/(y3+s_n))

    # Differential equation for rainfall dynamics
    def dRdt(t, y1, y2, y3):
        return r_R * ((R_const + b * y1 / k) - y2)


    def dSdt(t, y1, y2, y3):
        return -dTdt(t, y1, y2, y3)

    # System of differential equations for tree cover and rainfall
    def dzdt(t, y):
        return [dTdt(t, *y), dRdt(t, *y), dSdt(t, *y)]

    # Array of initial rainfall levels for simulation


    # Time span for the simulation
    t_span = (0, 600)
    # Time points at which to solve the system
    t_eval = np.linspace(t_span[0], t_span[1], 600)



    sol = solve_ivp(fun=lambda t, y: dzdt(t, y), 
                            y0=[T0, R0, S0], t_span=t_span, t_eval=t_eval, 
                            method='RK45')


    return t_eval, sol.y[0], sol.y[1], sol.y[2]

t, T, R, S = new_model(0.2, 0.1, 1, 10, 0)

plt.figure()
plt.plot(t, T, label='Tree cover')
plt.plot(t, R, label='Rainfall')
plt.plot(t, S, label='Cutted Tree Cover')
plt.xlabel('Time [years]')
plt.ylabel('Level')
plt.legend()
plt.show()

