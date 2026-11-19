import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1. Define the column names based on the 10 data columns
columns = [
    'time_t', 'position_z', 
    'fo_upwind_sine', 'fo_upwind_square', 
    'lax_sine', 'lax_square', 
    'MacCormack_sine', 'MacCormack_square', 
    'so_upwind_sine', 'so_upwind_square'
]

# 2. Define the parameter space and target times
n_values = [1, 3]
dz_values = [1.0, 0.5]
c0_values = [1.0, 0.6, 0.3]
target_times = [0, 2, 10, 18]
domain_length = 40.0
base_path = '/home/ar09/Documents/CFD_Assignments/assignment2/Q3/'

# 3. Iterate through all combinations of n, dz, and C0
for n in n_values:
    for dz in dz_values:
        for c0 in c0_values:
            # Dynamically modify filename for each parameter combination
            filename = f'{base_path}sol_dz_{dz}_C0_{c0}_n_{n}.dat'
            
            if not os.path.exists(filename):
                print(f"File not found, skipping: {filename}")
                continue
                
            df = pd.read_csv(filename, sep=r'\s+', names=columns, skiprows=1)
            
            # 4. Iterate over each requested time
            for t_target in target_times:
                # Find the nearest time available in the dataset
                idx_nearest = (df['time_t'] - t_target).abs().idxmin()
                actual_time = df.loc[idx_nearest, 'time_t']
                df_time = df[df['time_t'] == actual_time]
                
                x = df_time['position_z'].values
                
                # 5. Calculate Exact Solutions
                # PDE: u_t + u_x = 0 implies wave speed c = 1.
                # Shift position by (x - c*t). Modulo applies periodic boundary conditions.
                x_shifted = (x - actual_time) % domain_length
                
                # Apply initial condition formulas to the shifted spatial coordinates
                exact_sine = np.sin((2 * n * np.pi * x_shifted) / 40.0)
                exact_square = np.where(x_shifted <= 10, 1.0, 0.0)

                # 6. Set up the figure with 2 subplots (1 row, 2 columns)
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

                # --- Plot 1: Sine Wave ---
                ax1.plot(x, df_time['fo_upwind_sine'], label='FO Upwind', marker='o', markersize=3)
                ax1.plot(x, df_time['lax_sine'], label='Lax', marker='s', markersize=3)
                ax1.plot(x, df_time['MacCormack_sine'], label='MacCormack', marker='^', markersize=3)
                ax1.plot(x, df_time['so_upwind_sine'], label='SO Upwind', marker='d', markersize=3)
                ax1.plot(x, exact_sine, 'k--', linewidth=2, label='Exact Solution')

                ax1.set_title(f'Sine Wave (n={n}, dz={dz}, C0={c0}, t≈{actual_time:.2f})')
                ax1.set_xlabel('Position (x)')
                ax1.set_ylabel('Amplitude')
                ax1.legend()
                ax1.grid(True)

                # --- Plot 2: Square Wave ---
                ax2.plot(x, df_time['fo_upwind_square'], label='FO Upwind', marker='o', markersize=3)
                ax2.plot(x, df_time['lax_square'], label='Lax', marker='s', markersize=3)
                ax2.plot(x, df_time['MacCormack_square'], label='MacCormack', marker='^', markersize=3)
                ax2.plot(x, df_time['so_upwind_square'], label='SO Upwind', marker='d', markersize=3)
                ax2.plot(x, exact_square, 'k--', linewidth=2, label='Exact Solution')

                ax2.set_title(f'Square Wave (dz={dz}, C0={c0}, t≈{actual_time:.2f})')
                ax2.set_xlabel('Position (x)')
                ax2.set_ylabel('Amplitude')
                ax2.legend()
                ax2.grid(True)

                # 7. Save the plots
                plt.tight_layout()
                save_name = f'plot_dz_{dz}_C0_{c0}_n_{n}_t_{t_target}.png'
                plt.savefig(save_name)
                plt.close() # Close figure to prevent high memory usage
                print(f"Saved: {save_name}")
