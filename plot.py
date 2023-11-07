from model import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D





def plot_q11(T_arr, R_arr, t_span):

    R_arr = np.linspace(0, 5, 100)  # [mm/day]
    T0 = np.array([0.05, 0.25, 0.50, 0.75, 1])

    for i in range(len(T0)):

        T_init = T0[i]
        T_results = T_arr[:, i, :]

        plt.figure(figsize=(8, 6))
        
        # Plot the heatmap using imshow with the custom colormap
        extent = [t_span[0], t_span[1], R_arr.min(), R_arr.max()]
        im = plt.imshow(T_results.T, extent=extent, origin='lower', aspect='auto', cmap="viridis", vmax = 1, vmin = 0)
        
        # Add contour lines if desired
        levels = np.linspace(0, 1, num=100)  # Adjust number of levels
        
        # Add colorbar for reference
        plt.colorbar(im)
        
        # Label the axes and the plot
        plt.xlabel('Time [year]')
        plt.ylabel('R [mm/day]')
        plt.title(rf'$T_0 = {T_init}$', fontsize = 16)
        
        # Show the plot
        plt.show()
    return 


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
    ax.set_title(r"Equilibrium states $T(600)$ with respect to R and $T_0$", fontsize=16)

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
    ax.set_xlabel('R')
    ax.set_ylabel('T0')
    ax.set_zlabel('T(600)')
    ax.set_title(r"Equilibrium states of $T(600)$ with respect to $R_0$ and $T_0$", fontsize=16)

    # Adjust the aspect ratio
    # ax.set_box_aspect([np.ptp(i) for i in [R_mesh.flatten(), T0_mesh.flatten(), T600.flatten()]])

    # Set view angle for better 3D effect
    ax.view_init(elev=25, azim=130)  # Elevation and azimuth for good 3D view

    plt.tight_layout()  # Adjust the padding between and around subplots.
    plt.show()


if __name__ == "__main__":
    T_arr, R_arr, t_span = run_model1()

    plot_q11(T_arr, R_arr, t_span)
    plot_q12(T_arr, R_arr)

    T_arr = run_model2()

    plot_q2_3D(T_arr)


