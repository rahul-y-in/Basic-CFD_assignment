import numpy as np
import matplotlib.pyplot as plt

file_names = [
    'Temp_dist_ADI_dt_0.2.dat', 
    'Temp_dist_CN_imex_dt_0.2.dat', 
    'Temp_dist_implicit_dt_0.2.dat', 
    'Temp_dist_explicit_dt_0.2.dat'
]
labels = ['ADI', 'CN-IMEX', 'Implicit', 'Explicit']

# Distinct styles for overlapping lines
linestyles = ['-', '--', '-.', ':']
markers = ['o', 's', '^', 'x']
linewidths = [5, 3.5, 2, 1.5]  # Thicker lines first, thinner lines on top

target_x = 0.15
target_y = 0.2
dt = 0.2
tolerance = 1e-5  

plt.figure(figsize=(10, 6))

for i, (file, label) in enumerate(zip(file_names, labels)):
    try:
        data = np.loadtxt(file, skiprows=1)
        
        x_coords = data[:, 0]
        y_coords = data[:, 1]
        
        match_mask = np.isclose(x_coords, target_x, atol=tolerance) & np.isclose(y_coords, target_y, atol=tolerance)
        match_indices = np.where(match_mask)[0]
        
        if len(match_indices) == 0:
            print(f"Warning: Point ({target_x}, {target_y}) not found in {file}.")
            continue
            
        row_idx = match_indices[0]
        temp_history = data[row_idx, 2:]
        time_array = np.arange(len(temp_history)) * dt
        
        # Plot with distinct styles
        plt.plot(
            time_array, temp_history, 
            label=label,
            linestyle=linestyles[i],
            marker=markers[i],
            linewidth=linewidths[i],
            markersize=7,
            markerfacecolor='none'  # Makes markers hollow so lines underneath show through
        )
        
    except FileNotFoundError:
        print(f"Error: Could not find file '{file}'.")
    except Exception as e:
        print(f"An error occurred while processing '{file}': {e}")

# Formatting the plot
plt.xlabel('Time', fontsize=12)
plt.ylabel('Temperature', fontsize=12)
plt.title(f'Temperature Time History at (x={target_x}, y={target_y})', fontsize=14)
plt.grid(True, linestyle='-', alpha=0.3)
plt.legend()
plt.tight_layout()

# Display the plot
plt.show()
