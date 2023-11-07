from scipy.integrate import solve_ivp
import numpy as np




# Defining parameters
R = 2         # [mm/day]
r_m=0.3     # [1/year]
h_R=0.5     # [mm/day]
m_n=0.15    # [1/year]
h_n=0.1
m_f=0.11    # [1/year]
h_f=0.6
p=7
k=0.9
b=2         # [mm/day]
r_R=1       # [1/year]


def run_model1():

    r = lambda R: r_m * R / (h_R + R)


    def dTdt(t, T, R):
        return r(R) * T * (1 - T/k) - m_n * T * h_n / (T + h_n)  - m_f * T * (h_f**p) / (T**p + h_f**p)

    R_arr = np.linspace(0, 5, 100)  # [mm/day]
    T0 = np.array([0.05, 0.25, 0.50, 0.75, 1])


    t_span = (0, 600)
    t_eval = np.linspace(t_span[0], t_span[1], 600)

    T_arr = np.zeros((len(t_eval), len(T0), len(R_arr)))

    for i, R in enumerate(R_arr):
        sol = solve_ivp(fun=lambda t, y: dTdt(t, y, R), y0=T0, t_span=t_span, t_eval=t_eval, method='RK45')
        T_arr[:, :, i] = sol.y.T

    return T_arr, R_arr, t_span




def run_model2():
    R_const = 2

    r = lambda x: r_m * x / (h_R + x)

    def dTdt(t, y1, y2):
        return r(y2) * y1 * (1 - y1/k) - m_n * y1 * h_n / (y1 + h_n)  - m_f * y1 * (h_f**p) / (y1**p + h_f**p)

    def dRdt(t, y1, y2):
        return r_R * ((R_const + b * y1 / k) - y2)

    def dzdt(t, y):
        return [dTdt(t, *y), dRdt(t, *y)]

    R0 = np.linspace(0, 5, 100)  # [mm/day]
    T0 = np.linspace(0, 1, 30)

    t_span = (0, 600)
    t_eval = np.linspace(t_span[0], t_span[1], 600)

    T_arr = np.zeros((len(t_eval), len(T0), len(R0)))
    R_arr = np.zeros((len(t_eval), len(T0), len(R0)))


    for i, r0 in enumerate(R0):
        for j, t0 in enumerate(T0):
            sol = solve_ivp(fun=lambda t, y: dzdt(t, y), y0=[t0, r0] , t_span=t_span, t_eval=t_eval, method='RK45')
            T_arr[:, j, i] = sol.y[0]
            R_arr[:, j, i] = sol.y[1] 




    return T_arr

