from scipy.integrate import solve_ivp
import numpy as np

# Constants representing the parameters for the model
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

def run_model1():
    """
    Simulates the tree cover dynamics without considering the rainfall dynamics.
    
    Returns:
    T_arr : Array containing tree cover levels over time for different rainfall rates.
    R_arr : Array of rainfall rates [mm/day].
    t_span : Tuple representing the time interval of the simulation [days].
    """
    # Growth rate of tree cover as a function of rainfall
    r = lambda R: r_m * R / (h_R + R)

    # Differential equation for tree cover growth
    def dTdt(t, T, R):
        return (r(R) * T * (1 - T / k) - 
                m_n * T * h_n / (T + h_n) -
                m_f * T * (h_f ** p) / (T ** p + h_f ** p))

    # Array of rainfall rates for simulation
    R_arr = np.linspace(0, 5, 100)
    # Initial tree cover levels
    T0 = np.array([0.05, 0.25, 0.50, 0.75, 1])

    # Time span for the simulation
    t_span = (0, 600)
    # Time points at which to solve the system
    t_eval = np.linspace(t_span[0], t_span[1], 600)

    # Array to store the tree cover levels
    T_arr = np.zeros((len(t_eval), len(T0), len(R_arr)))

    # Solve the differential equation for each rainfall value
    for i, R_value in enumerate(R_arr):
        sol = solve_ivp(fun=lambda t, y: dTdt(t, y, R_value), 
                        y0=T0, t_span=t_span, t_eval=t_eval, 
                        method='RK45')
        T_arr[:, :, i] = sol.y.T

    return T_arr, R_arr, t_span

def run_model2():
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
    def dTdt(t, y1, y2):
        return (r(y2) * y1 * (1 - y1 / k) - 
                m_n * y1 * h_n / (y1 + h_n) -
                m_f * y1 * (h_f ** p) / (y1 ** p + h_f ** p))

    # Differential equation for rainfall dynamics
    def dRdt(t, y1, y2):
        return r_R * ((R_const + b * y1 / k) - y2)

    # System of differential equations for tree cover and rainfall
    def dzdt(t, y):
        return [dTdt(t, *y), dRdt(t, *y)]

    # Array of initial rainfall levels for simulation
    R0 = np.linspace(0, 5, 100)
    # Array of initial tree cover levels
    T0 = np.linspace(0, 1, 30)

    # Time span for the simulation
    t_span = (0, 600)
    # Time points at which to solve the system
    t_eval = np.linspace(t_span[0], t_span[1], 600)

    # Arrays to store the tree cover and rainfall levels
    T_arr = np.zeros((len(t_eval), len(T0), len(R0)))
    R_arr = np.zeros((len(t_eval), len(T0), len(R0)))

    # Solve the differential equations for each initial condition
    for i, r0 in enumerate(R0):
        for j, t0 in enumerate(T0):
            sol = solve_ivp(fun=lambda t, y: dzdt(t, y), 
                            y0=[t0, r0], t_span=t_span, t_eval=t_eval, 
                            method='RK45')
            T_arr[:, j, i] = sol.y[0]
            R_arr[:, j, i] = sol.y[1]

    return T_arr
