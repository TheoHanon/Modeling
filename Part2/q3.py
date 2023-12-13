from scipy.integrate import solve_ivp
import matplotlib.animation as animation
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import scipy.io
import cartopy.crs as ccrs
import cartopy.feature as cfeature




R = 2           # Amount of rainfall [mm/day]
r_m = 0.3       # Maximum growth rate of tree cover [1/year]
h_R = 0.5       # Half-saturation constant for rainfall effect on growth [mm/day]
m_n = 0.15      # Natural mortality rate of trees [1/year]
h_n = 0.1       # Half-saturation constant for natural mortality
m_f = 0.11      # Mortality rate of trees due to external factors [1/year]
h_f = 0.6       # Half-saturation constant for external mortality factors
p = 7           # Hill coefficient for external mortality response
k = 0.9         # Carrying capacity of tree cover
b = 2           # Baseline tree reproduction rate [mm/day]
r_R = 1         # Rainfall regeneration rate [1/year]

### New parameters
r_H = .01       # Human Growth Rate Factor [1/year]
K = .8          # Carrying capacity of humans [-]
pi1 = .05       # Proportional constant of the use of forestry resources [-]
alpha1 = .1     # Max Cutting Tree Rate caused by Human Population [1/year]

alpha2 = 0.01   # Max Cutting Tree Rate caused by Human Demand [1/year]
lam = .1        # Growth Rate of Demand due to the Population [1/year]
lam0 = .05      # Proportion of the Demand fulfill by other alternatives [1/year]




### Define the new model 


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
        return lam*y3*(1-y4) - lam0*y4

    # System of differential equations for tree cover and rainfall
    def dzdt(t, y):
        return [dTdt(t, *y), dRdt(t, *y), dHdt(t, *y), dPdt(t, *y)]

    # Array of initial rainfall levels for simulation

    nYear = 200
    # Time span for the simulation
    t_span = (0, nYear)
    # Time points at which to solve the system
    t_eval = np.linspace(t_span[0], t_span[1], nYear)


    sol = solve_ivp(fun=lambda t, y: dzdt(t, y), 
                            y0=[T0, R0, H0, P0], t_span=t_span, t_eval=t_eval, 
                            method='RK45')


    return t_eval, sol.y[0], sol.y[1], sol.y[2], sol.y[3]


### Import data


def import_data():

    fname_mat = 'data/2023worldforestrainfalldata.mat'
    mat_file = scipy.io.loadmat(fname_mat)

    forest1999 = mat_file['forest1999']
    rainfall1999 = mat_file['rainfall1999']
    rainfall2100 = mat_file['rainfall2100']
    row,col=rainfall1999.shape

    lats = np.squeeze(mat_file['lats'], axis=0)
    lons = np.squeeze(mat_file['lons'],axis=0)
    lons_sa = [-90, -30]  

    return forest1999, rainfall1999, rainfall2100, row, col, lats, lons, lons_sa


### Approximating Inital Human density


def approx_human_density(forest1999, rainfall1999, lats, lons):


    human_density = np.full((len(lats), len(lons)), np.nan)
    human_density[~(np.isnan(forest1999))] = (100 - forest1999[~(np.isnan(forest1999))])*.4
    human_density[(np.isnan(forest1999)) ^ (np.isnan(rainfall1999))] = 30

    return human_density


### Simulating the model across South America


def simulation():

    forest1999, rainfall1999, rainfall2100, row, col, lats, lons, lons_sa = import_data()
    human_density = approx_human_density(forest1999, rainfall1999, lats, lons) 

    tree_cover_pred = np.full((len(lats), len(lons), 200), np.nan)
    human_density_pred = np.full((len(lats), len(lons), 200), np.nan)
    demand_pred = np.full((len(lats), len(lons), 200), np.nan)



    for i in range(row):
        for j in range(col):
            if(lons_sa[0]<=lons[j]<=lons_sa[1]):
                if (np.isnan(human_density[i, j])): continue
                H0 = human_density[i, j]*0.01 
                T0 = forest1999[i, j]*.01 if ~np.isnan(forest1999[i, j]) else 0.0
                R0 = rainfall2100[i, j]*0.001 if ~np.isnan(rainfall2100[i, j]) else 0.0
                P0 = 0

                _ , T , R, H, P = new_model(T0, R0, H0, P0)
                tree_cover_pred[i, j, :] = T
                human_density_pred[i, j, :] = H
                demand_pred[i, j, :] = P


    return tree_cover_pred, human_density_pred, demand_pred, lats, lons

### Define plotting functions

def plot_simulation():

    tree_cover_pred, human_density_pred, demand_pred, lats, lons = simulation()

    ### Final state of the simulation

    fig = plt.figure(figsize=(10, 5))

    ax1 = fig.add_subplot(131, projection=ccrs.PlateCarree())
    ax1.set_extent([-90, -30, min(lats), max(lats)], crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.OCEAN, alpha=0.7)
    contour1 = ax1.contourf(lons, lats, tree_cover_pred[:,:, -1], cmap='Greens', vmin=0, vmax=1, levels=np.arange(0, 1.05, .05))
    ax1.set_title('Predicted Tree Cover in 2150', fontsize=14)
    ax1.coastlines(resolution='110m', color='black', linewidth=1)
    ax1.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    plt.colorbar(contour1, ax=ax1, orientation='horizontal',pad=0.1, boundaries=np.arange(0, 1.1, .1))

    ax2 = fig.add_subplot(132, projection=ccrs.PlateCarree())
    ax2.set_extent([-90, -30, min(lats), max(lats)], crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.OCEAN, alpha=0.7)
    contour2 = ax2.contourf(lons, lats, human_density_pred[:,:, -1], cmap='Reds', vmin=0, vmax=1, levels=np.arange(0, 1.05, .05))
    ax2.set_title('Predicted Human Density in 2150', fontsize=14)
    ax2.coastlines(resolution='110m', color='black', linewidth=1)
    ax2.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    plt.colorbar(contour2, ax=ax2, orientation='horizontal',pad=0.1, boundaries=np.arange(0, 1.1, .1))


    ax3 = fig.add_subplot(133, projection=ccrs.PlateCarree())
    ax3.set_extent([-90, -30, min(lats), max(lats)], crs=ccrs.PlateCarree())
    ax3.add_feature(cfeature.OCEAN, alpha=0.7)
    contour3 = ax3.contourf(lons, lats, demand_pred[:,:, -1], cmap='Blues', vmin=0, vmax=1, levels=np.arange(0, 1.05, .05))
    ax3.set_title('Predicted Demand in 2150', fontsize=14)
    ax3.coastlines(resolution='110m', color='black', linewidth=1)
    ax3.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    plt.colorbar(contour3, ax=ax3, orientation='horizontal',pad=0.1, boundaries=np.arange(0, 1.1, .1))
    
    plt.show()


    ### Full simulation of the model


    fig = plt.figure(figsize=(7, 4))
    plt.suptitle("Simulation of Tree Cover & Human Density from 1999 to 2150", fontsize=14)
    ax1 = fig.add_subplot(121, projection=ccrs.PlateCarree())
    ax1.set_extent([-90, -30, min(lats), max(lats)], crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.OCEAN, alpha=0.7)
    contour1 = ax1.contourf(lons, lats, tree_cover_pred[:,:, 0], cmap='Greens', vmin=0, vmax=1,levels=np.arange(0, 1.05, .05))#,levels=np.linspace(0,30, 100))
    ax1.set_title('Forest Cover in 1999', fontsize=14)
    ax1.coastlines(resolution='110m', color='black', linewidth=1)
    ax1.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    plt.colorbar(contour1, ax=ax1, orientation='horizontal',pad=0.1, boundaries = np.arange(0, 1.1, .1))

    ax2 = fig.add_subplot(122, projection=ccrs.PlateCarree())
    ax2.set_extent([-90, -30, min(lats), max(lats)], crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.OCEAN, alpha=0.7)
    contour2 = ax2.contourf(lons, lats, human_density_pred[:,:, 0], cmap='Reds', vmin=0, vmax=1,levels=np.arange(0, 1.05, .05))
    ax2.set_title('Human Density in 1999', fontsize=14)
    ax2.coastlines(resolution='110m', color='black', linewidth=1)
    ax2.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    plt.colorbar(contour2, ax=ax2, orientation='horizontal',pad=0.1,boundaries=np.arange(0, 1.1, .1))


    # # Update function for the animation
    def update_plot(frame_number, tree_cover_pred, human_density_pred, ax1, ax2):
        ax1.set_title(f'Forest Cover in {1999 + frame_number}', fontsize=14)  # Update the year dynamically
        contour1 = ax1.contourf(lons, lats, tree_cover_pred[:, :, frame_number], cmap='Greens', vmin=0, vmax=1,levels=np.arange(0, 1.05, .05))

        ax2.set_title(f'Human Density in {1999 + frame_number}', fontsize=14)  # Update the year dynamically
        contour2 = ax2.contourf(lons, lats, human_density_pred[:, :, frame_number], cmap='Reds', vmin=0, vmax=1,levels=np.arange(0, 1.05, .05))

        return contour1, contour2,

    ani = animation.FuncAnimation(fig, update_plot, frames=152,interval=100, fargs=(tree_cover_pred, human_density_pred, ax1, ax2))
    # ani.save('plot/ani.mp4', writer='ffmpeg') # Uncomment to save the animation
    plt.show()


    return




