import numpy as np
from matplotlib.widgets import Slider
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import scipy.io
import cartopy.crs as ccrs
import cartopy.feature as cfeature


R = 2         # Amount of rainfall [mm/day]
r_m = 0.3    # Maximum growth rate of tree cover [1/year]
h_R = 0.5     # Half-saturation constant for rainfall effect on growth [mm/day]
m_n = 0.15    # Natural mortality rate of trees [1/year]
h_n = 0.1     # Half-saturation constant for natural mortality
m_f = 0.11   # Mortality rate of trees due to external factors [1/year]
h_f = 0.6     # Half-saturation constant for external mortality factors
p = 7         # Hill coefficient for external mortality response
k = 0.9       # Carrying capacity of tree cover
b = 2         # Baseline tree reproduction rate [mm/day]
r_R = 1      # Rainfall regeneration rate [1/year]


# Growth rate of tree cover as a function of rainfall
r = lambda x: r_m * x / (h_R + x)

### New parameters
alpha1= 0.1
r_H = .01  # Human growth rate [1/year]
K = .8  # Carrying capacity of humans 
pi1 = .1
alpha2 = 0.01
lam = .01
lam0 = .1
R_const = 2
# Differential equation for tree cover growth
def dTdt(y1, y3 , y4,K):
    return (r(R_const) * y1 * (1 - y1 / k) - 
            m_n * y1 * h_n / (y1 + h_n) -
            m_f * y1 * (h_f ** p) / (y1 ** p + h_f ** p) - alpha1*y1*y3 - alpha2* y1**2 * y4)

def dHdt(y1, y3, y4,K):
    return r_H * y3*(1 - y3/K) + pi1 * alpha1* y1*y3

def dPdt(y1, y3, y4, K):
    return lam*y3 - lam0*y4

t = np.linspace(0, 200, 201)

def computeVariables(t,T0,H0,P0,K):
    T = [T0]
    H = [H0]
    P = [P0]
    for j in range(1, len(t)):
        next_T = T[j-1] + dTdt(T[j-1], H[j-1], P[-1],K)
        next_H = H[j-1] + dHdt(T[j-1], H[j-1], P[-1],K)
        next_P = P[j-1] + dPdt(T[j-1], H[j-1], P[-1],K)
        T.append(next_T)
        H.append(next_H)
        P.append(next_P)
    T=100*(np.array(T))
    H=100*(np.array(H))
    P=100*(np.array(P))
    return T, H, P

t = np.linspace(0, 200, 201)

# En pleine forêt --> K petit, T grand, H et P petits
K=0.1
T0=0.9
H0=0.05
P0=0.05

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

T, H, P = computeVariables(t,T0, H0, P0,K)
line, = ax.plot(T, H, P)
line2, = ax.plot(T[0], H[0], P[0], 'or', label='t = 0 [day]')
line3, = ax.plot(T[-1], H[-1], P[-1], 'og', label='t = 200 [day]')

ax.set_xlabel('Tree cover [%]')
ax.set_ylabel('Population density [%]')
ax.set_zlabel('Population demand [%]')
ax.legend()
fig.suptitle("Phase portrait of the model in the heart of the forest", fontsize=15)
fig.subplots_adjust(left=0.25, bottom=0.25)
ax.set_ylim(0, 100)  
ax.set_xlim(0, 100) 
ax.set_zlim(0, 100)   
plt.show()

# Loin de la forêt --> K grand, T petit, H et P grands
K=0.8
T0=0.15
H0=0.5
P0=0.3

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

T, H, P = computeVariables(t,T0, H0, P0,K)
line, = ax.plot(T, H, P)
line2, = ax.plot(T[0], H[0], P[0], 'or', label='t = 0 [day]')
line3, = ax.plot(T[-1], H[-1], P[-1], 'og', label='t = 200 [day]')

ax.set_xlabel('Tree cover [%]')
ax.set_ylabel('Population density [%]')
ax.set_zlabel('Population demand [%]')
ax.legend()
fig.suptitle("Phase portrait of the model far from the forest", fontsize=15)
fig.subplots_adjust(left=0.25, bottom=0.25)
ax.set_ylim(0, 100)  
ax.set_xlim(0, 100) 
ax.set_zlim(0, 100)   
plt.show()


# Pas trop loin de la forêt --> K moyen, T, H et P moyen
K=0.5
T0=0.5
H0=0.2
P0=0.15

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

T, H, P = computeVariables(t,T0, H0, P0,K)
line, = ax.plot(T, H, P)
line2, = ax.plot(T[0], H[0], P[0], 'or', label='t = 0 [day]')
line3, = ax.plot(T[-1], H[-1], P[-1], 'og', label='t = 200 [day]')

ax.set_xlabel('Tree cover [%]')
ax.set_ylabel('Population density [%]')
ax.set_zlabel('Population demand [%]')
ax.legend()
fig.suptitle("Phase portrait of the model near the forest", fontsize=15)
fig.subplots_adjust(left=0.25, bottom=0.25)
ax.set_ylim(0, 100)  
ax.set_xlim(0, 100) 
ax.set_zlim(0, 100)   
plt.show()



