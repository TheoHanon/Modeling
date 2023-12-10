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

### New parameters
r_H = .01 #Human growth rate [1/year]
K = .8 #Carrying capacity of humans 
pi1 = .1
alpha1 = .1 # Max cutting tree cover / year 

alpha2 = 0.01
lam = .01
lam0 = .1

def new_model(T0, R0, H0, P0):
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
    def dTdt(t, y1, y2, y3, y4):
        return (r(y2) * y1 * (1 - y1 / k) - 
                m_n * y1 * h_n / (y1 + h_n) -
                m_f * y1 * (h_f ** p) / (y1 ** p + h_f ** p) - alpha1*y1*y3 - alpha2* y1**2 * y4)

    # Differential equation for rainfall dynamics
    def dRdt(t, y1, y2, y3, y4):
        return 0

    def dHdt(t, y1, y2, y3, y4):
        return r_H * y3*(1 - y3/K) + pi1 * alpha1* y1*y3

    def dPdt(t, y1, y2, y3, y4):
        return lam*(y3-y4) - lam0*y4

    # System of differential equations for tree cover and rainfall
    def dzdt(t, y):
        return [dTdt(t, *y), dRdt(t, *y), dHdt(t, *y), dPdt(t, *y)]

    # Array of initial rainfall levels for simulation

    nYear = 600
    # Time span for the simulation
    t_span = (0, nYear)
    # Time points at which to solve the system
    t_eval = np.linspace(t_span[0], t_span[1], nYear)


    sol = solve_ivp(fun=lambda t, y: dzdt(t, y), 
                            y0=[T0, R0, H0, P0], t_span=t_span, t_eval=t_eval, 
                            method='RK45')


    return t_eval, sol.y[0], sol.y[1], sol.y[2], sol.y[3]

# t, T, R, H, P = new_model(.1, 3, .5, 0)


# plt.figure()
# plt.plot(t, T, label='Tree cover')
# plt.plot(t, R, label='Rainfall')
# plt.plot(t, H, label = 'Human population')
# plt.plot(t, P, label = 'Human demand')
# # plt.xticks(np.arange(0, 600*365, 365), np.arange(0, 600, 1))
# # plt.plot(t, U, label='Rainfall')
# plt.xlabel('Time [years]')
# plt.ylabel('Level')
# plt.legend()
# plt.show()



fname_mat = './2023worldforestrainfalldata.mat'
mat_file = scipy.io.loadmat(fname_mat)

forest1999 = mat_file['forest1999']
rainfall1999 = mat_file['rainfall1999']
rainfall2100 = mat_file['rainfall2100']
row,col=rainfall1999.shape
# print(np.isnan(rainfall1999))
lats = np.squeeze(mat_file['lats'], axis=0)
lons = np.squeeze(mat_file['lons'],axis=0)
lons_sa = [-90, -30]  



human_density = np.full((len(lats), len(lons)), np.nan)

human_density[~(np.isnan(forest1999))] = (100 - forest1999[~(np.isnan(forest1999))])*.4
human_density[(np.isnan(forest1999)) ^ (np.isnan(rainfall1999))] = 50


tree_cover_pred = np.full((len(lats), len(lons), 200), np.nan)
human_density_pred = np.full((len(lats), len(lons), 200), np.nan)
for i in range(row):
    for j in range(col):
        if(lons_sa[0]<=lons[j]<=lons_sa[1]):
            if (np.isnan(human_density[i, j])): continue
            H0 = human_density[i, j]*.01 
            T0 = forest1999[i, j]*.01 if ~np.isnan(forest1999[i, j]) else 0.0
            R0 = rainfall2100[i, j]*0.001 if ~np.isnan(rainfall2100[i, j]) else 0.0
            P0 = 0
            # print(H0, T0, R0, P0)
            _ , T , _, H, _ = new_model(T0, R0, H0, P0)
            tree_cover_pred[i, j, :] = T
            human_density_pred[i, j, :] = H

            # print(T0, T[0])



print("Done")

# fig = plt.figure(figsize=(10, 10))
# ax1 = fig.add_subplot(111, projection=ccrs.PlateCarree())

# # Forest cover for 1999
# ax1.set_extent([-90,-30, min(lats), max(lats)], crs=ccrs.PlateCarree())
# ax1.add_feature(cfeature.OCEAN, alpha = .7)
# cs1 = ax1.contourf(lons, lats, tree_cover_pred, cmap='Greens', vmin=0, vmax=1)
# ax1.set_title('Forest Cover in 1999', fontsize = 14)
# ax1.coastlines(resolution='110m', color='black', linewidth=1)
# ax1.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

# plt.show()

# fig = plt.figure(figsize=(10, 10))
# ax1 = fig.add_subplot(111, projection=ccrs.PlateCarree())

# # Forest cover for 1999
# ax1.set_extent([-90,-30, min(lats), max(lats)], crs=ccrs.PlateCarree())
# ax1.add_feature(cfeature.OCEAN, alpha = .7)
# cs1 = ax1.contourf(lons, lats, human_density_pred, cmap='OrRd', vmin=0, vmax=1)
# ax1.set_title('Forest Cover in 1999', fontsize = 14)
# ax1.coastlines(resolution='110m', color='black', linewidth=1)
# ax1.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

# plt.show()



fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(111, projection=ccrs.PlateCarree())
ax1.set_extent([-90, -30, min(lats), max(lats)], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.OCEAN, alpha=0.7)
contour = ax1.contourf(lons, lats, tree_cover_pred[:, :, 0], cmap='OrRd', vmin=0, vmax=1,levels=np.linspace(0,1, 100))
ax1.set_title('Forest Cover in 1999', fontsize=14)
ax1.coastlines(resolution='110m', color='black', linewidth=1)
ax1.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

# Update function for the animation
def update_plot(frame_number, tree_cover_pred, ax):
    # ax.clear()
    ax.set_extent([-90, -30, min(lats), max(lats)], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN, alpha=0.7)
    ax.coastlines(resolution='110m', color='black', linewidth=1)
    ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    ax.set_title(f'Forest Cover in {1999 + frame_number}', fontsize=14)  # Update the year dynamically
    contour = ax.contourf(lons, lats, tree_cover_pred[:, :, frame_number], cmap='OrRd', vmin=0, vmax=.8,levels=np.linspace(0,1, 100))
    return contour,

ani = animation.FuncAnimation(fig, update_plot, frames=200,interval=20, fargs=(tree_cover_pred, ax1))
# plt.show()
# Create the animation
# ani = animation.FuncAnimation(fig, update_plot, frames=100, fargs=(tree_cover_pred, ax1, fig))
ani.save('test2.mp4', fps=30)