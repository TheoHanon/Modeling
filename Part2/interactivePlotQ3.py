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

# Differential equation for tree cover growth
def dTdt(y1, y3 , y4, alpha1,alpha2,K,r_H, pi1,lam,lam0,R_const):
    return (r(R_const) * y1 * (1 - y1 / k) - 
            m_n * y1 * h_n / (y1 + h_n) -
            m_f * y1 * (h_f ** p) / (y1 ** p + h_f ** p) - alpha1*y1*y3 - alpha2* y1**2 * y4)


def dHdt(y1, y3, y4, alpha1,alpha2,K,r_H, pi1,lam,lam0,R_const):
    return r_H * y3*(1 - y3/K) + pi1 * alpha1* y1*y3

def dPdt(y1, y3, y4, alpha1,alpha2,K,r_H, pi1,lam,lam0,R_const):
    return lam*y3*(1-y4) - lam0*y4

t = np.linspace(0, 200, 201)

def computeVariables(t,alpha1,alpha2,K,r_H, pi1,lam,lamb0,R_const):
    T = [0.1]
    H = [.3]
    P = [0]
    for j in range(1, len(t)):
        next_T = T[j-1] + dTdt(T[j-1], H[j-1], P[-1],alpha1,alpha2,K,r_H, pi1,lam,lamb0,R_const)
        next_H = H[j-1] + dHdt(T[j-1], H[j-1], P[-1],alpha1,alpha2,K,r_H, pi1,lam,lamb0,R_const)
        next_P = P[j-1] + dPdt(T[j-1], H[j-1], P[-1],alpha1,alpha2,K,r_H, pi1,lam,lamb0,R_const)
        T.append(next_T)
        H.append(next_H)
        P.append(next_P)
    T=100*(np.array(T))
    H=np.array(H)*100
    P=100*(np.array(P))
    return T, H, P


### New parameters
init_alpha1= 0.3 
init_r_H = .01  # Human growth rate [1/year]
init_K = .8  # Carrying capacity of humans 
init_pi1 = .05
init_alpha2 = 0.1
init_lam = .01
init_lam0 = .01
init_R_const = 2

fig, ax = plt.subplots(figsize=(12, 7))
line1, = ax.plot(t, computeVariables(t, init_alpha1, init_alpha2, init_K,init_r_H,init_pi1,init_lam,init_lam0,init_R_const)[0], lw=2, label='Tree Cover [%]')
line2, = ax.plot(t, computeVariables(t, init_alpha1, init_alpha2, init_K,init_r_H,init_pi1,init_lam,init_lam0,init_R_const)[1], lw=2, label='Human Population [%]')
line3, = ax.plot(t, computeVariables(t, init_alpha1, init_alpha2, init_K,init_r_H,init_pi1,init_lam,init_lam0,init_R_const)[2], lw=2, label='Population pressure [%]')
ax.set_xlabel('Time [year]')
ax.legend()
fig.suptitle("Dynamics of the deforestation caused by human population pressure",fontsize=15)
fig.subplots_adjust(left=0.25, bottom=0.25)

ax.set_ylim(0, 100)  

axAlpha1 = fig.add_axes([0.1, 0.1, 0.15, 0.03])
alpha1_slider = Slider(
    ax=axAlpha1,
    label=' $\\alpha_1$',
    valmin=0,
    valmax=1,
    valinit=0.3,
)
axAlpha2 = fig.add_axes([0.3, 0.1, 0.15, 0.03])
alpha2_slider = Slider(
    ax=axAlpha2,
    label=' $\\alpha_2$',
    valmin=0,
    valmax=1,
    valinit=0.1,
)
axK = fig.add_axes([0.5, 0.1, 0.15, 0.03])
K_slider = Slider(
    ax=axK,
    label=' $K$',
    valmin=0,
    valmax=1,
    valinit=.8,
)
axRh = fig.add_axes([0.7, 0.1, 0.15, 0.03])
Rh_slider = Slider(
    ax=axRh,
    label=' $r_H$',
    valmin=0,
    valmax=1,
    valinit=0.01,
)
axpi1 = fig.add_axes([0.1, 0.05, 0.15, 0.03])
pi1_slider = Slider(
    ax=axpi1,
    label=' $\\pi_1$',
    valmin=0,
    valmax=.2,
    valinit=0.05,
)
axlambd = fig.add_axes([0.3, 0.05, 0.15, 0.03])
lambd_slider = Slider(
    ax=axlambd,
    label=' $\\lambda$',
    valmin=0,
    valmax=1,
    valinit=0.1,
)
axlambd0 = fig.add_axes([0.5, 0.05, 0.15, 0.03])
lambd0_slider = Slider(
    ax=axlambd0,
    label=' $\\lambda_0$',
    valmin=0,
    valmax=1,
    valinit=.01,
)
axR = fig.add_axes([0.7, 0.05, 0.15, 0.03])
R_slider = Slider(
    ax=axR,
    label=' $R$',
    valmin=0,
    valmax=5,
    valinit=2,
)

def update(val):
    alpha1 = alpha1_slider.val
    alpha2 = alpha2_slider.val
    K=K_slider.val
    r_H=Rh_slider.val
    pi1=pi1_slider.val
    lam=lambd_slider.val
    lam0=lambd0_slider.val
    R_const=R_slider.val
    res=computeVariables(t, alpha1,alpha2,K,r_H, pi1,lam,lam0,R_const)
    line1.set_ydata(res[0])
    line2.set_ydata(res[1])
    line3.set_ydata(res[2])

alpha1_slider.on_changed(update)
alpha2_slider.on_changed(update)
K_slider.on_changed(update)
Rh_slider.on_changed(update)
pi1_slider.on_changed(update)
lambd_slider.on_changed(update)
lambd0_slider.on_changed(update)
R_slider.on_changed(update)
plt.show()