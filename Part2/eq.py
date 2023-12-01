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
lons_sa = [-90, -30] # South america coordinates (longitude range)
lons_af = [-30, 60] # African coordinates (longitude range)
lons_au = [60, 160] 

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
                T_init = forest1999[i,j]
                if worstcase:
                    forest2100[i,j] = T_arr[T_arr.argmin()]
                else:
                    forest2100[i,j] = T_arr[T_arr.argmax()]
     
            
    return forest2100


forest2100_america_wc = compute_forest2100(fixedpointdata_southamerica, lons_sa, worstcase = True)
forest2100_america_bs = compute_forest2100(fixedpointdata_southamerica, lons_sa, worstcase = False)


forest2100_africa_wc = compute_forest2100(fixedpointdata_africa, lons_af, worstcase = True)
forest2100_africa_bs = compute_forest2100(fixedpointdata_africa, lons_af, worstcase = False)
forest2100_aus_wc = compute_forest2100(fixedpointdata_aus, lons_au, worstcase = True)
forest2100_aus_bs = compute_forest2100(fixedpointdata_aus, lons_au, worstcase = False)

forest2100_wc = [forest2100_america_wc, forest2100_africa_wc, forest2100_aus_wc]
forest2100_bs = [forest2100_america_bs, forest2100_africa_bs, forest2100_aus_bs]



def plot_forest_cover(forest_covers, case):

    fig = plt.figure(figsize=(15, 12))

    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 0.05], hspace=0.05)

    ax1 = fig.add_subplot(gs[0], projection=ccrs.PlateCarree())
    # First subplot for the year 1999
    ax1.set_extent([min(lons), max(lons), min(lats), max(lats)], crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.OCEAN, alpha = .7)
    cs1 = ax1.contourf(lons, lats, forest1999, cmap='Greens', vmin=0, vmax=80, levels=np.linspace(0, 80, 80))
    ax1.set_title('Forest Cover in 1999', fontsize = 14)
    ax1.coastlines(resolution='110m', color='black', linewidth=1)
    ax1.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    # Second subplot for the year 2100
    ax2 = fig.add_subplot(gs[1], projection=ccrs.PlateCarree())
    ax2.set_extent([min(lons), max(lons), min(lats), max(lats)], crs=ccrs.PlateCarree())
    for forest in forest_covers:
        cs2 = ax2.contourf(lons, lats, forest, cmap="Greens", vmin=0, vmax=80, levels=np.linspace(0, 80, 80))
    ax2.add_feature(cfeature.OCEAN, alpha = 0.7)
    ax2.set_title(f'Predicted {case} Case Forest Cover in 2100', fontsize = 14)
    ax2.coastlines(resolution='110m', color='black', linewidth=1)
    ax2.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.7, linestyle='--')

    cax = fig.add_subplot(gs[2])
    cbar = plt.colorbar(cs2, cax=cax, orientation='horizontal', pad = 0.01)
    cbar.set_label('Forest Cover (%)')
    cbar_ticks = np.arange(0, 81, 10)  # Creates a list from 0 to 80 in steps of 10
    cbar.set_ticks(cbar_ticks)  # Sets the positions of the ticks

    # Optionally, set custom tick labels
    cbar_ticklabels = [f'{tick}%' for tick in cbar_ticks]  # List comprehension to add '%' to each tick
    cbar.set_ticklabels(cbar_ticklabels)  # Sets the text of the ticks
    fig.subplots_adjust(top=0.95, bottom=0.35, left=0.05, right=0.95, hspace=0.02)
    plt.show()

    return 


plot_forest_cover(forest2100_wc, case = 'Worst')
plot_forest_cover(forest2100_bs, case = 'Best')


