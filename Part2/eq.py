import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import UnivariateSpline
import math 
import matplotlib.pyplot as plt
import numpy as np
import scipy.io


fname = './treecovdata.nc'
fname_mat = './2023worldforestrainfalldata.mat'
mat_file = scipy.io.loadmat(fname_mat)

rainfall1999 = mat_file['rainfall1999']
rainfall2014 = mat_file['rainfall2014']
rainfall2100 = mat_file['rainfall2100']
forest1999 = mat_file['forest1999']
lats = np.squeeze(mat_file['lats'], axis=0)
lons = np.squeeze(mat_file['lons'],axis=0)

row,col=rainfall2100.shape

fname_mat = './worldfixedpointdata.mat'
mat_file = scipy.io.loadmat(fname_mat)
fixedpointdata_southamerica = mat_file['fixedpointdata_southamerica']
fixedpointdata_africa = mat_file['fixedpointdata_africa']
fixedpointdata_aus = mat_file['fixedpointdata_aus']

# longitudes of southamerica, africa and australia [lower boundary value, laster boundary value]
lons_sa = [-90, -30]    # South america coordinates (longitude range)
lons_af = [-30, 60]     # African coordinates (longitude range)
lons_au = [60, 160]     # Australia coordinates (longitude range)

def unpack(x):
    return x[0][0]

vunpack = np.vectorize(unpack)


# we have values of R at equilibirum only for multiples of 25
def get_rainfall_index(r):
    available_rainfall = np.arange(-25,2526,50)  # from -25 to 2525 with step 50
    return np.argmin(np.abs(available_rainfall-r))


def compute_forest2100(fixedpointdata, Lons, worstcase = False):
    forest2100=np.full((row,col),np.nan)
    for i in range (row):
        for j in range(col):
            if(Lons[0]<=lons[j]<=Lons[1]):
                if(np.isnan(rainfall2100[i,j])): continue
                index = get_rainfall_index(rainfall2100[i,j])
                T_arr = vunpack(fixedpointdata[index][1:3])
                T_arr = T_arr[np.logical_not(np.isnan(T_arr))]
                if worstcase:
                    forest2100[i,j] = T_arr[T_arr.argmin()]
                else:
                    forest2100[i,j] = T_arr[T_arr.argmax()]
     
            
    return forest2100


forest2100_america_wc = compute_forest2100(fixedpointdata_southamerica, lons_sa, worstcase = True)
forest2100_america_bc = compute_forest2100(fixedpointdata_southamerica, lons_sa, worstcase = False)
forest2100_africa_wc = compute_forest2100(fixedpointdata_africa, lons_af, worstcase = True)
forest2100_africa_bc = compute_forest2100(fixedpointdata_africa, lons_af, worstcase = False)
forest2100_aus_wc = compute_forest2100(fixedpointdata_aus, lons_au, worstcase = True)
forest2100_aus_bc = compute_forest2100(fixedpointdata_aus, lons_au, worstcase = False)

forest2100_wc = [forest2100_america_wc, forest2100_africa_wc, forest2100_aus_wc]
forest2100_bc = [forest2100_america_bc, forest2100_africa_bc, forest2100_aus_bc]



def plot_rainfall(rainfall1999, rainfall2100):
        fig = plt.figure(figsize=(10, 10))


        ax1 = fig.add_subplot(211, projection=ccrs.PlateCarree())

        ax1.set_extent([min(lons), max(lons), min(lats), max(lats)], crs=ccrs.PlateCarree())
        ax1.add_feature(cfeature.OCEAN, alpha = .7)
        cs1 = ax1.contourf(lons, lats, rainfall1999, cmap='Blues', vmax = 7000, levels=np.linspace(0, 7000, 7))#, vmin=0, vmax=80, levels=np.linspace(0, 80, 80))
        ax1.set_title('Rainfall in 1999', fontsize = 14)
        ax1.coastlines(resolution='110m', color='black', linewidth=1)
        ax1.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        # cbar = plt.colorbar(cs1, ax=ax1, orientation='horizontal', pad = 0.1,boundaries = np.linspace(0, 5000, 100))
        # cbar.set_label('Rainfall (mm/year))')

        ax2 = fig.add_subplot(212, projection=ccrs.PlateCarree())
        
        ax2.set_extent([min(lons), max(lons), min(lats), max(lats)], crs=ccrs.PlateCarree())
        ax2.add_feature(cfeature.OCEAN, alpha = .7)
        cs2 = ax2.contourf(lons, lats, rainfall2100, cmap='Blues', vmax = 7000,levels=np.linspace(0, 7000, 7))#, vmin=0, vmax=80, levels=np.linspace(0, 80, 80))
        ax2.set_title('Rainfall in 2100', fontsize = 14)
        ax2.coastlines(resolution='110m', color='black', linewidth=1)
        ax2.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        cbar = plt.colorbar(cs1, ax=ax2, orientation='horizontal', pad = 0.1, boundaries = np.linspace(0, 7000, 7), aspect = 40)
        cbar.set_label('Rainfall (mm/year)')
        cbar_ticks = np.arange(0, 8000, 1000)  # Creates a list from 0 to 80 in steps of 10
        cbar.set_ticks(cbar_ticks)

        plt.show()







def plot_forest_cover(forest_wc,forest_bc):

    fig = plt.figure(figsize=(10, 10))

    gs = gridspec.GridSpec(4, 1, height_ratios=[1, 1, 1,0.05])

    ax1 = fig.add_subplot(gs[0], projection=ccrs.PlateCarree())

    # Forest cover for 1999
    ax1.set_extent([min(lons), max(lons), min(lats), max(lats)], crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.OCEAN, alpha = .7)
    cs1 = ax1.contourf(lons, lats, forest1999, cmap='Greens', vmin=0, vmax=80, levels=np.linspace(0, 80, 80))
    ax1.set_title('Forest Cover in 1999', fontsize = 14)
    ax1.coastlines(resolution='110m', color='black', linewidth=1)
    ax1.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    # Forest cover for 2100 - best case
    ax2 = fig.add_subplot(gs[1], projection=ccrs.PlateCarree())
    ax2.set_extent([min(lons), max(lons), min(lats), max(lats)], crs=ccrs.PlateCarree())
    for forest in forest_bc:
        cs2 = ax2.contourf(lons, lats, forest, cmap="Greens", vmin=0, vmax=80, levels=np.linspace(0, 80, 80))
    ax2.add_feature(cfeature.OCEAN, alpha = 0.7)
    ax2.set_title(f'Predicted Best Case Forest Cover in 2100', fontsize = 14)
    ax2.coastlines(resolution='110m', color='black', linewidth=1)
    ax2.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.7, linestyle='--')

    # Forest cover for 2100 - worst case
    ax3 = fig.add_subplot(gs[2], projection=ccrs.PlateCarree())
    ax3.set_extent([min(lons), max(lons), min(lats), max(lats)], crs=ccrs.PlateCarree())
    for forest in forest_wc:
        cs3 = ax3.contourf(lons, lats, forest, cmap="Greens", vmin=0, vmax=80, levels=np.linspace(0, 80, 80))
    ax3.add_feature(cfeature.OCEAN, alpha = 0.7)
    ax3.set_title(f'Predicted Worst Case Forest Cover in 2100', fontsize = 14)
    ax3.coastlines(resolution='110m', color='black', linewidth=1)
    ax3.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.7, linestyle='--')

    cax = fig.add_subplot(gs[3])
    cbar = plt.colorbar(cs2, cax=cax, orientation='horizontal', pad = 0.01)
    cbar.set_label('Forest Cover (%)')
    cbar_ticks = np.arange(0, 81, 10)  # Creates a list from 0 to 80 in steps of 10
    cbar.set_ticks(cbar_ticks)  # Sets the positions of the ticks

    # Optionally, set custom tick labels
    cbar_ticklabels = [f'{tick}%' for tick in cbar_ticks]  # List comprehension to add '%' to each tick
    cbar.set_ticklabels(cbar_ticklabels)  # Sets the text of the ticks
    fig.subplots_adjust(top=0.95, bottom=0.1, left=0.05, right=0.95, hspace=0.2)
    plt.savefig("PredictedTreeCover.pdf")
    plt.show()
    return 

plot_rainfall(rainfall1999,rainfall2100)
plot_forest_cover(forest2100_wc,forest2100_bc)



