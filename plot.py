from model import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap


def create_smooth_colormap():
    # Define colors as (R, G, B)
    sandybrown = (244/255, 164/255, 96/255)
    gold = (255/255, 215/255, 0/255)
    forestgreen = (34/255, 139/255, 34/255)

    # Create a colormap with smooth transitions
    cdict = {
        'red':   [(0.0, sandybrown[0], sandybrown[0]),
                    (0.05, sandybrown[0], sandybrown[0]),
                  (0.06, gold[0], gold[0]),
                  (0.7, gold[0], gold[0]),
                  (0.71, forestgreen[0], forestgreen[0]),
                  (1.0, forestgreen[0], forestgreen[0])],  # Ensure it ends with x=1
        'green': [(0.0, sandybrown[1], sandybrown[1]),
                    (0.05, sandybrown[1], sandybrown[1]),
                  (0.06, gold[1], gold[1]),
                  (0.7, gold[1], gold[1]),
                  (0.71, forestgreen[1], forestgreen[1]),
                  (1.0, forestgreen[1], forestgreen[1])],  # Ensure it ends with x=1
        'blue':  [(0.0, sandybrown[2], sandybrown[2]),
                (0.05, sandybrown[2], sandybrown[2]),
                  (0.06, gold[2], gold[2]),
                  (0.7, gold[2], gold[2]),
                  (0.71, forestgreen[2], forestgreen[2]),
                  (1.0, forestgreen[2], forestgreen[2])]   # Ensure it ends with x=1
    }
    
    return LinearSegmentedColormap('CustomMap', cdict)

def apply_custom_colormap(T_results, colormap):
    # Normalize T_results to be in the range 0-1
    T_norm = T_results / T_results.max()

    # Use the colormap to get the corresponding color for each value
    return colormap(T_norm)

def plot_q11(T_arr, t_span):
    R_arr = np.linspace(0, 5, 100)  # Define the range of R values
    T0 = np.array([0.05, 0.25, 0.50, 0.75, 1])  # Define the initial T values

    # Create a custom colormap
    custom_colormap = create_smooth_colormap()

    for i, T_init in enumerate(T0):
        T_results = T_arr[:, i, :]  # Get the results for the current T_init

        plt.figure(figsize=(8, 6))
        extent = [t_span[0], t_span[1], R_arr.min(), R_arr.max()]
        
        # Plot the heatmap using imshow with the smooth custom colormap
        im = plt.imshow(T_results.T, extent=extent, origin='lower', aspect='auto', cmap=custom_colormap, vmax=1, vmin=0)
        
        # Add colorbar for reference
        cbar = plt.colorbar(im)

        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label='Arid (T < 0.05)', markersize=15, markerfacecolor='sandybrown', linestyle='None'),
            Line2D([0], [0], marker='o', color='w', label='Savannah (0.05 ≤ T ≤ 0.7)', markersize=15, markerfacecolor='gold', linestyle='None'),
            Line2D([0], [0], marker='o', color='w', label='Forest (T > 0.7)', markersize=15, markerfacecolor='forestgreen', linestyle='None'),
        ]
        plt.legend(handles=legend_elements, loc='upper right', shadow=True)        
        # Label the axes and the plot
        plt.xlabel('Time [year]')
        plt.ylabel('R [mm/day]')
        plt.title(rf'$T_0={T_init}$', fontsize=25)
        
        # Show the plot
        plt.show()
# Rest of the code remains unchanged...



def plot_q12(T_arr, R_arr):


    R_arr = np.linspace(0, 5, 100)  # [mm/day]
    T0 = np.array([0.05, 0.25, 0.50, 0.75, 1])

    T600 = T_arr[-1, :, :].flatten()
    R_mesh, T0_mesh = np.meshgrid(R_arr, T0)

    # Define the color for each point based on the T600 value
    colors = np.where(T600 < 0.05, 'sandybrown', np.where(T600 <= 0.7, 'gold', 'forestgreen'))

    # Create a 3D scatter plot
    fig = plt.figure(figsize=(8, 6))  # Set figure size
    ax = fig.add_subplot(111, projection='3d', adjustable='box')

    # Increase marker size and make them semi-transparent for better visibility
    sc = ax.scatter(R_mesh.flatten(), T0_mesh.flatten(), T600.flatten(), c=colors, s=35, alpha=0.6)

    # Custom legend with more descriptive labels and larger icons
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Arid (T < 0.05)', markersize=15, markerfacecolor='sandybrown', linestyle='None'),
        Line2D([0], [0], marker='o', color='w', label='Savannah (0.05 ≤ T ≤ 0.7)', markersize=15, markerfacecolor='gold', linestyle='None'),
        Line2D([0], [0], marker='o', color='w', label='Forest (T > 0.7)', markersize=15, markerfacecolor='forestgreen', linestyle='None'),
    ]
    ax.legend(handles=legend_elements, shadow=True)

    # Grid lines and 3D effect
    ax.grid(True)  # Show grid lines for better orientation
    ax.set_xlabel('R')
    ax.set_ylabel('T0')
    ax.set_zlabel('T(600)')
    # ax.set_title(r"Equilibrium states $T(600)$ with respect to R and $T_0$", fontsize=16)

    # Adjust the aspect ratio
    # ax.set_box_aspect([np.ptp(i) for i in [R_mesh.flatten(), T0_mesh.flatten(), T600.flatten()]])

    # Set view angle for better 3D effect
    ax.view_init(elev=25, azim=130)  # Elevation and azimuth for good 3D view

    plt.tight_layout()  # Adjust the padding between and around subplots.
    plt.show()

    return 



def plot_q2_3D(T_arr):


    R0 = np.linspace(0, 5, 100)  # [mm/day]
    T0 = np.linspace(0, 1, 30)

    R_mesh, T0_mesh = np.meshgrid(R0, T0)
    T600 = T_arr[-1, :, :]  # Assuming T_arr is a 3D array with the correct shape


    fig = plt.figure(figsize=(8, 6))  # Set figure size
    ax = fig.add_subplot(111, projection='3d', adjustable='box')

    # Scatter plot
    colors = np.where(T600.flatten() < 0.05, 'sandybrown', np.where(T600.flatten() <= 0.7, 'gold', 'forestgreen'))
    sc = ax.scatter(R_mesh.flatten(), T0_mesh.flatten(), T600.flatten(), c=colors, s=35, alpha=0.6)


    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Arid (T < 0.05)', markersize=15, markerfacecolor='sandybrown', linestyle='None'),
        Line2D([0], [0], marker='o', color='w', label='Savannah (0.05 ≤ T ≤ 0.7)', markersize=15, markerfacecolor='gold', linestyle='None'),
        Line2D([0], [0], marker='o', color='w', label='Forest (T > 0.7)', markersize=15, markerfacecolor='forestgreen', linestyle='None'),
    ]
    ax.legend(handles=legend_elements, shadow=True)

    # Grid lines and 3D effect
    ax.grid(True)  # Show grid lines for better orientation
    ax.set_xlabel(r'R_0')
    ax.set_ylabel(r'T_0')
    ax.set_zlabel('T(600)')
    # ax.set_title(r"Equilibrium states $T(600)$ with respect to $R_0$ and $T_0$", fontsize=16)

    # Adjust the aspect ratio
    # ax.set_box_aspect([np.ptp(i) for i in [R_mesh.flatten(), T0_mesh.flatten(), T600.flatten()]])

    # Set view angle for better 3D effect
    ax.view_init(elev=25, azim=130)  # Elevation and azimuth for good 3D view

    plt.tight_layout()  # Adjust the padding between and around subplots.
    plt.show()


if __name__ == "__main__":
    T_arr, R_arr, t_span = run_model1()

    plot_q11(T_arr,t_span)
    plot_q12(T_arr, R_arr)

    T_arr = run_model2()

    plot_q2_3D(T_arr)


