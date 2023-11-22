import numpy as np
import matplotlib.pyplot as plt
import sympy as sy
from matplotlib.lines import Line2D

R, T, z = sy.symbols('R T z')

# Constants representing the parameters for the model

r_m = 0.3     # Maximum growth rate of tree cover [1/year]
h_R = 0.5     # Half-saturation constant for rainfall effect on growth [mm/day]
m_n = 0.15    # Natural mortality rate of trees [1/year]
h_n = 0.1     # Half-saturation constant for natural mortality
m_f = 0.11    # Mortality rate of trees due to external factors [1/year]
h_f = 0.6     # Half-saturation constant for external mortality factors
p = 7         # Hill coefficient for external mortality response
k = 0.9       # Carrying capacity of tree cover
b = 2         # Baseline tree reproduction rate [mm/day]
r_R = 1 
    # Rainfall regeneration rate [1/year]



r = lambda x: r_m * x / (h_R + x)

dTdt = r(R) * T * (1 - T/k) - m_n * T * h_n/ (h_n + T) - m_f * T * h_f**p / (T**p + h_f**p)

J = lambda t, x : sy.diff(dTdt, T).subs(R, x).subs(T, t)

F = lambda x: sy.Poly(sy.expand(r(x) * z*(1 - z/k)*(z+h_n) *(z**p + h_f**p) - m_n*z*h_n*(z**p + h_f**p) - m_f*z*(z + h_n)*h_f**p))

a = lambda x: F(x).all_coeffs()



# J = lambda T, R: r(R)*(1 - 2*T/k) - m_n*h_n/(T + h_n)**2 - m_f*h_f**p * (T**p - p*T**(p-1) + h_f**p) / (T**p + h_f**p)**2




# a = lambda x: [-r(x)/k, r(x), r(x)*h_n*(1-1/k)-m_n*h_n, 0, 0, 0, 0, -r(x)*h_f**p / k, r(x)*h_f**p -m_f*h_f**p, -(r(x)*h_n*h_f**p *(1/k - 1) + m_n*h_n*h_f**p+m_f*h_f**p*h_n),0]

R_init = np.linspace(0, 5, 200)

for r_init in R_init:
    T_eq = np.roots(a(r_init))
    for t_eq in T_eq:
        if (np.isreal(t_eq) and t_eq.real >=0):
            if (J(t_eq, r_init) < 0):
                plt.plot(r_init, t_eq.real, 'ro')
            elif (J(t_eq, r_init) > 0):
                plt.plot(r_init, t_eq.real, 'bo')
            else:
                plt.plot(r_init, t_eq.real, 'go')
                
plt.legend(handles = [Line2D([0], [0], marker='o', color='r', label='attracting node', markerfacecolor='r', markersize=10), Line2D([0],[0], marker='o', color='b', label='repelling node', markerfacecolor='b', markersize=10)])
plt.xlabel(r'$\bar{R}$')
plt.ylabel(r'$\bar{T}$')
plt.show()





