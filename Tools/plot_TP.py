import glob
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

# Set seaborn style
sns.set_theme()
sns.set_style("ticks")

# Define a set of distinct colors that are easier on the eyes
cyberpunk_colors = [
    '#E64980',  # Muted pink
    '#37B24D',  # Forest green 
    '#845EF7',  # Soft purple
    '#3BC9DB',  # Light blue
    '#F76707',  # Burnt orange
    '#94D82D',  # Olive green
    '#DA77F2',  # Soft fuchsia
    '#FAB005'   # Golden yellow
]

# Get all .txt files in the specified directory
tp_files = glob.glob('./HELIOS_output/bench_comp5/*.txt')

# Create empty lists to store data from all files
all_data = []
file_names = []

# Read each file
for file in tp_files:
    try:
        # Load data while skipping the header lines
        data = np.loadtxt(file, skiprows=2)  # Skip the first two header lines
        all_data.append(data)
        # Get just the filename without path for legend
        file_names.append(os.path.basename(file))
    except Exception as e:
        print(f"Could not read file: {file}. Error: {e}")

# If no files found, print a message and exit
if not tp_files:
    print("No .txt files found in the specified directory.")
else:
    # Plot each temperature profile with distinct colors and increased line thickness
    for i, data in enumerate(all_data):
        # Check if filename contains 100, 101, or 102 and use black color
        if any(str(num) in file_names[i] for num in [103, 104, 105]):
            color = 'black'
            zord = -10
        else:
            color = cyberpunk_colors[i % len(cyberpunk_colors)]
            zord = 100
        plt.plot(data[:, 1], data[:, 0] / 1e6, label=file_names[i], color=color, linewidth=2.5, zorder=zord)

    # List of suffixes to plot - modify this list to choose which profiles to show
    suffixes = ['comp5_f025', 'comp5_f033', 'comp5_f067']

    # Plot each 55cnce profile with the same color set
    for i, suffix in enumerate(suffixes):
        filepath = f'../Results/55cnce/helios_benchmark/55cnce_{suffix}/NewTemperature.dat'
        try:
            # Skip first line (header), use column 2 for pressure and 3 for temperature
            data = np.loadtxt(filepath, skiprows=1, usecols=(1, 2))
            plt.plot(data[:, 1], 10**(data[:, 0] - 5), label=f'55cnce_{suffix}', linestyle='--', color=cyberpunk_colors[i % len(cyberpunk_colors)], linewidth=2.5)
        except Exception as e:
            print(f"Could not read file: {filepath}. Error: {e}")

    # Add labels and title
    plt.ylabel('Pressure (bar)')  # Y-axis label
    plt.xlabel('Temperature (K)')  # X-axis label
    plt.title('Temperature Profiles')
    plt.legend()  # Show legend with filenames
    plt.yscale('log')  # Set y-axis to log scale
    plt.ylim(1e-8, 1e1)  # Set y-axis limits from 1e-8 to 10 bar
    plt.gca().invert_yaxis()  # Invert the pressure axis
    plt.tick_params(axis='both', direction='in')  # Make ticks point inwards

    os.makedirs('output', exist_ok=True)
    # Save the figure in high resolution
    plt.savefig('output/temperature_profiles_bench_comp5.png', dpi=400, bbox_inches='tight')

    plt.show()  # Display the plot

