import numpy as np
import matplotlib.pyplot as plt
import sympy as sy
from matplotlib.lines import Line2D

### Constants representing the parameters for the model

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


r = lambda x: r_m * x / (h_R + x)

### Defining the Jacobian of the system & the polynomial whose roots give the equilibria

J = lambda T, R: r(R)*(1 - 2*T/k) - m_n*h_n**2/(T + h_n)**2 - m_f*h_f**p * (T**p - p*T**p + h_f**p) / (T**p + h_f**p)**2
a = lambda x: [-r(x)/k, r(x)*(1-h_n/k), h_n*r(x)-m_n*h_n, 0, 0, 0, 0, -r(x)*h_f**p / k, h_f**p *(r(x) - h_n / k *r(x) - m_f), (r(x)*h_n*h_f**p - m_n*h_n*h_f**p-m_f*h_f**p*h_n),0]

### Function to plot the equilibria

def plot_equilibria():

    R_init = np.linspace(0, 5, 200)

    for r_init in R_init:
        T_eq = np.roots(a(r_init))
        for t_eq in T_eq:
            if (np.isreal(t_eq) and t_eq.real >=0):
                if (J(t_eq, r_init) < 0):
                    plt.plot(r_init, t_eq.real, 'go')
                elif (J(t_eq, r_init) > 0):
                    plt.plot(r_init, t_eq.real, 'ro')
                else:
                    plt.plot(r_init, t_eq.real, 'bo')
                    
    plt.legend(handles = [Line2D([0], [0], marker='o', color='g', label='attractive node', markerfacecolor='g', markersize=10), Line2D([0],[0], marker='o', color='r', label='repulsive node', markerfacecolor='r', markersize=10)])
    plt.xlabel(r'$\overline{R}$[mm/day]')
    plt.ylabel(r'$\overline{T}$ [%]')
    plt.title(r"Graph of the equilibria $(\overline{T},\overline{R})$")
    plt.show()

    return 


