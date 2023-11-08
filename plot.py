from model import *  # Import all functions and variables from model script
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def plot_q11(T_arr, t_span):
    """
    Plot a series of heatmaps showing the change in tree cover over time for different initial conditions and rainfall rates.

    Parameters:
    T_arr : ndarray
        The array containing the tree cover values over time and rainfall rates.
    t_span : tuple
        The time span over which the simulation was run.
    """
    # Define the range of rainfall rates and initial tree cover values
    R_arr = np.linspace(0, 5, 100)
    T0 = np.array([0.05, 0.25, 0.50, 0.75, 1])

    # Loop over initial tree cover values to generate heatmaps
    for i, T_init in enumerate(T0):
        # Extract results for the current initial tree cover value
        T_results = T_arr[:, i, :]  
        
        # Initialize plot
        plt.figure(figsize=(8, 6))
        extent = [t_span[0], t_span[1], R_arr.min(), R_arr.max()]
        
        # Create and display heatmap
        im = plt.imshow(T_results.T, extent=extent, origin='lower', aspect='auto', cmap="viridis", vmax=1, vmin=0)
        cbar = plt.colorbar(im)
        cbar.set_label('Tree Cover (T)', rotation=270, labelpad=15)
        plt.xlabel('Time [year]')
        plt.ylabel('Rainfall (R) [mm/day]')
        plt.title(f'Initial Tree Cover $T_0={T_init}$', fontsize=25)
        plt.show()

def plot_q12(T_arr, R_arr):
    """
    Create a 3D scatter plot representing the final tree cover values against initial tree cover and rainfall rates.

    Parameters:
    T_arr : ndarray
        The array containing the tree cover values over time and rainfall rates.
    R_arr : ndarray
        The array of rainfall rates used in the simulation.
    """
    # Extract final tree cover values for all simulations
    T600 = T_arr[-1, :, :].flatten()
    T0 = np.array([0.05, 0.25, 0.50, 0.75, 1])
    
    R_mesh, T0_mesh = np.meshgrid(R_arr, T0)

    # Determine color for points based on final tree cover
    colors = np.where(T600 < 0.05, 'sandybrown', np.where(T600 <= 0.7, 'gold', 'forestgreen'))

    # Set up 3D scatter plot
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d', adjustable='box')
    sc = ax.scatter(R_mesh.flatten(), T0_mesh.flatten(), T600.flatten(), c=colors, s=35, alpha=0.6)

    # Define custom legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Arid (T < 0.05)', markersize=15, markerfacecolor='sandybrown', linestyle='None'),
        Line2D([0], [0], marker='o', color='w', label='Savannah (0.05 ≤ T ≤ 0.7)', markersize=15, markerfacecolor='gold', linestyle='None'),
        Line2D([0], [0], marker='o', color='w', label='Forest (T > 0.7)', markersize=15, markerfacecolor='forestgreen', linestyle='None'),
    ]
    ax.legend(handles=legend_elements, shadow=True)

    # Add axis labels and grid
    ax.grid(True)
    ax.set_xlabel('Rainfall (R) [mm/day]')
    ax.set_ylabel('Initial Tree Cover (T0)')
    ax.set_zlabel('Tree Cover at T=600')
    ax.view_init(elev=25, azim=130)
    plt.tight_layout()
    plt.show()

def plot_q12_add():
    """
    Plot the rate of change of tree cover (dT/dt) over a range of tree cover values for different rainfall rates.
    """
    # Lambda function for rate of change of tree cover
    r = lambda x: r_m * x / (h_R + x)

    # Define range of rainfall rates and tree cover values for plotting
    R_arr = np.linspace(0, 5, 6)
    T_array = np.linspace(0, 1, 100)

    # Plot dT/dt for each rainfall rate
    for i, R in enumerate(R_arr):
        plt.plot(T_array, r(R) * T_array * (1 - T_array/k) - m_n * T_array * h_n / (T_array + h_n) - m_f * T_array * (h_f**p) / (T_array**p + h_f**p), label=f"R={R}")
    plt.axhline(y=0, color='m', linestyle='--')
    plt.xlabel("Tree Cover (T)")
    plt.ylabel("Rate of Change of Tree Cover (dT/dt)")
    plt.legend()
    plt.show()

def plot_q2_3D(T_arr):
    """
    Create a 3D scatter plot showing the distribution of final tree cover against initial rainfall and initial tree cover.
    """
    # Define range of initial rainfall rates and initial tree cover values
    R0 = np.linspace(0, 5, 100)
    T0 = np.linspace(0, 1, 30)

    # Generate meshgrid for the 3D scatter plot
    R_mesh, T0_mesh = np.meshgrid(R0, T0)
    T600 = T_arr[-1, :, :]  # Extract final tree cover values

    # Set up the figure for 3D plotting
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d', adjustable='box')

    # Determine color for points and create scatter plot
    colors = np.where(T600.flatten() < 0.05, 'sandybrown', np.where(T600.flatten() <= 0.7, 'gold', 'forestgreen'))
    sc = ax.scatter(R_mesh.flatten(), T0_mesh.flatten(), T600.flatten(), c=colors, s=35, alpha=0.6)

    # Define custom legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Arid (T < 0.05)', markersize=15, markerfacecolor='sandybrown', linestyle='None'),
        Line2D([0], [0], marker='o', color='w', label='Savannah (0.05 ≤ T ≤ 0.7)', markersize=15, markerfacecolor='gold', linestyle='None'),
        Line2D([0], [0], marker='o', color='w', label='Forest (T > 0.7)', markersize=15, markerfacecolor='forestgreen', linestyle='None'),
    ]
    ax.legend(handles=legend_elements, shadow=True)

    # Add axis labels, grid, and set view angle for better visualization
    ax.grid(True)
    ax.set_xlabel('Initial Rainfall (R0) [mm/day]')
    ax.set_ylabel('Initial Tree Cover (T0)')
    ax.set_zlabel('Final Tree Cover (T600)')
    ax.view_init(elev=25, azim=130)
    plt.tight_layout()
    plt.show()

# The main block runs the model and calls the plotting functions.
if __name__ == "__main__":
    # Run the models and get the results
    T_arr, R_arr, t_span = run_model1()
    
    # Call the plotting functions with the results
    plot_q11(T_arr, t_span)
    plot_q12(T_arr, R_arr)
    plot_q12_add()
    
    # Run the second model
    T_arr = run_model2()
    
    # Plot the results in 3D
    plot_q2_3D(T_arr)
