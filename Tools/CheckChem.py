import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# Create necessary directories
os.makedirs('output', exist_ok=True)

# Read the data file
case_name = 'AlphaAb'
case_number = 100  # Default case number for FastChem data

filename = f'../Results/{case_name}/ConcentrationSTD_T.dat'
temp_profile_file = f'../Results/{case_name}/NewTemperature.dat'
helios_profile_file = f'HELIOS_output/{case_name}/00100.txt'
species_file = '../Condition/SpeciesList/species_helios_comp3.dat'
fastchem_file = f'HELIOS_output/{case_name}/fc_00100/EQ.txt'  # Define FastChem file path here
cloud_file = f'../Results/{case_name}/Diagnostic_condens.dat'  # Add cloud data file
data = pd.read_csv(filename, sep='\t', skiprows=2, header=None)

# Read cloud data
try:
    cloud_data = pd.read_csv(cloud_file, sep='\t', header=0)
    print(f"Successfully loaded cloud data from {cloud_file}")
    print("Cloud data columns:")
    for col in cloud_data.columns:
        print(f"  {col}")
except Exception as e:
    print(f"Warning: Could not load cloud data: {e}")
    cloud_data = None

# Read temperature profile data
try:
    temp_profile = pd.read_csv(temp_profile_file, sep='\t', comment='#', header=None,
                              names=['Height', 'LogP', 'Temperature'])
    print(f"Successfully loaded EPACRIS temperature profile from {temp_profile_file}")
except Exception as e:
    print(f"Warning: Could not load EPACRIS temperature profile: {e}")
    temp_profile = None

# Read HELIOS temperature profile
try:
    # Skip the first row (units) and use the second row as header
    helios_profile = pd.read_csv(helios_profile_file, sep='\t', skiprows=1)
    
    # Convert HELIOS pressure from dyne/cm² to bar
    # 1 bar = 1e6 dyne/cm²
    helios_profile['Pressure_bar'] = helios_profile['Pressure'] / 1e6
    
    print(f"Successfully loaded HELIOS temperature profile from {helios_profile_file}")
    print(f"HELIOS pressure range: {helios_profile['Pressure_bar'].min():.2e} - {helios_profile['Pressure_bar'].max():.2e} bar")
    print(f"HELIOS temperature range: {helios_profile['Temp'].min():.1f} - {helios_profile['Temp'].max():.1f} K")
except Exception as e:
    print(f"Warning: Could not load HELIOS temperature profile: {e}")
    helios_profile = None

# Boltzmann constant in J/K
kb = 1.380658E-23

# Read the species list file to map species numbers to names
species_names = {}
try:    
    # Read the file line by line to properly handle the format
    with open(species_file, 'r') as f:
        # Skip the header line
        header = f.readline()
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                species_name = parts[0]
                try:
                    species_number = int(parts[2])
                    species_names[species_number] = species_name
                except ValueError:
                    continue
    print(f"Successfully loaded species names from {species_file}")
except Exception as e:
    print(f"Warning: Could not load species names: {e}")
    print("Will use species numbers instead of names for plotting")

# Print all species names and their respective numbers
print("\nSpecies Numbers and Names:")
print("--------------------------")
if species_names:
    # Sort by species number for better readability
    for species_num in sorted(species_names.keys()):
        print(f"Species {species_num}: {species_names[species_num]}")
else:
    print("No species names available.")
print("--------------------------\n")


########################################################################################
# Extract data columns
# temperature = data.iloc[:, 3].values  # Temperature (K)
# pressure = data.iloc[:, 4].values  # Pressure (Pa)

# If ConcentrationSTD_T.dat already contains VMR, read them directly
pressure = data.iloc[:, 4].values  # Pressure (Pa)
temperature = data.iloc[:, 3].values  # Temperature (K)

# Create a dictionary to store species VMR data directly
species_vmr = {}

# Process all columns from index 5 onwards (these are the species VMRs)
for col_idx in range(5, data.shape[1]):
    species_idx = col_idx - 5 + 1  # Add 1 to make species index start from 1
    if species_idx <= 111:
        species_vmr[species_idx] = data.iloc[:, col_idx].values
    else:
        print(f"Warning: Skipping column {col_idx} as it would create species index {species_idx} > 111")

# Print the first species column (column index 5, species index 1)
first_species_idx = 7  # This should be 1, not 112
first_species_name = species_names.get(first_species_idx, f"Species {first_species_idx}")
print(f"\nFirst species column ({first_species_name}, ID: {first_species_idx}):")
print("--------------------------")
if first_species_idx in species_vmr:
    print(species_vmr[first_species_idx])
else:
    print(f"Error: Species index {first_species_idx} not found in species_vmr!")
    print(f"Available species indices: {sorted(species_vmr.keys())}")
print("--------------------------\n")

# Print some information about the extracted data
print(f"Temperature range: {temperature.min():.2f} - {temperature.max():.2f} K")
print(f"Pressure range: {pressure.min()/1.0e5:.2e} - {pressure.max()/1.0e5:.2e} bar")
print(f"Number of species: {len(species_vmr)}")

# VMR is already read directly
vmr = species_vmr

# Check if the number of species in VMR matches the number in species_vmr
debug_vmr = vmr
print(f"\nNumber of species in species_vmr: {len(species_vmr)}")
print(f"Number of species in VMR: {len(vmr)}")

if len(vmr) == len(species_vmr):
    print("✓ The number of species in VMR matches the number in species_vmr")
else:
    print("⚠ Warning: The number of species in VMR does not match the number in species_vmr!")
    missing_in_vmr = set(species_vmr.keys()) - set(vmr.keys())
    if missing_in_vmr:
        print(f"Species in species_vmr but not in VMR: {missing_in_vmr}")
    missing_in_species_vmr = set(vmr.keys()) - set(species_vmr.keys())
    if missing_in_species_vmr:
        print(f"Species in VMR but not in species_vmr: {missing_in_species_vmr}")

print("--------------------------\n")

# Extract cloud molecular abundances for species 7 (water)
cloud_molec_abundance = None
cloud_pressure_bar = None

# if cloud_data is not None:
#     try:
#         # Find the column with nCloudsMolec[7]
#         cloud_molec_col = [col for col in cloud_data.columns if 'Molec' in col and '7' in col]
#         if cloud_molec_col:
#             cloud_molec_col = cloud_molec_col[0]
#             print(f"Found cloud molecular abundance column: {cloud_molec_col}")
#             cloud_molec_abundance = cloud_data[cloud_molec_col].values
            
#             # Extract pressure data
#             if 'log(P)' in cloud_data.columns:
#                 # Convert log10(Pa) to bar
#                 cloud_pressure_bar = 10**cloud_data['log(P)'].values / 1.0e5
#                 print(f"Successfully extracted cloud data for plotting")
#                 print(f"Cloud pressure range: {cloud_pressure_bar.min():.2e} - {cloud_pressure_bar.max():.2e} bar")
#                 print(f"Non-zero cloud molecular abundances: {np.sum(cloud_molec_abundance > 0)}")
#         else:
#             print("Warning: Could not find cloud molecular abundance column for species 7")
#     except Exception as e:
#         print(f"Error processing cloud data: {e}")

# Define colors for cycling through - using a larger set of distinct colors
# Combine basic colors with more sophisticated colormaps for better distinction
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# Create an extended color palette with more distinct colors
basic_colors = ['k', 'r', 'b', 'g', 'm', 'c', 'y']
extended_colors = ['navy', 'darkgreen', 'crimson', 'purple', 'orange', 'teal', 'brown',
                  'olive', 'darkred', 'steelblue', 'darkviolet', 'goldenrod', 'indigo', 
                  'lightseagreen', 'saddlebrown', 'darkorange', 'forestgreen', 'midnightblue',
                  'darkslategray', 'maroon', 'darkgoldenrod', 'slateblue', 'sienna', 'mediumvioletred']

# Add tableau colors for even more variety
tableau_colors = list(mcolors.TABLEAU_COLORS.values())

# Add some colors from the viridis colormap for good measure (perceptually uniform)
viridis_colors = [cm.viridis(i/10) for i in range(0, 10, 2)]
viridis_hex = [mcolors.rgb2hex(color) for color in viridis_colors]

# Combine all color sets, prioritizing the most distinct ones
COLORS = basic_colors + extended_colors + tableau_colors + viridis_hex

# Create a function to get a unique color for a species index or name
def get_unique_color(species_name_or_idx, used_colors=None):
    if used_colors is None:
        used_colors = []
    
    # For string inputs, hash the name to get a stable index
    if isinstance(species_name_or_idx, str):
        # Use hash to create a stable index based on the species name
        hash_val = hash(species_name_or_idx.lower())
        idx = abs(hash_val) % len(COLORS)
    else:
        idx = species_name_or_idx % len(COLORS)
    
    # If this color is already used, find the next available one
    base_color = COLORS[idx]
    if base_color in used_colors:
        # Find the next color that hasn't been used
        for i in range(len(COLORS)):
            alt_idx = (idx + i) % len(COLORS)
            if COLORS[alt_idx] not in used_colors:
                return COLORS[alt_idx]
        # If all colors are used, just return the original (this should not happen with our large palette)
        return base_color
    return base_color

# Find the top 15 species by average VMR
top_species = {}
for species_idx in vmr:
    if not np.all(vmr[species_idx] == 0):
        avg_vmr = np.mean(vmr[species_idx][vmr[species_idx] > 0])
        top_species[species_idx] = avg_vmr
        species_name = species_names.get(species_idx, f"Species {species_idx}")
        print(f"{species_name} (ID: {species_idx}, column {species_idx + 4}): Average VMR = {avg_vmr:.2e}")

# Sort species by average VMR and get top 15
top_15_species = sorted(top_species.items(), key=lambda x: x[1], reverse=True)[:15]

# Convert pressure from Pa to bar (1 bar = 100,000 Pa)
pressure_bar = pressure / 1.0e5

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))

# Left panel - Volume mixing ratios from EPACRIS
# Create a dictionary to store colors by species name for consistent coloring
species_colors = {}
used_colors = []  # Track which colors we've used to avoid repetition

# Plot only the top 15 species from EPACRIS
for i, (species_idx, avg_vmr) in enumerate(top_15_species):
    species_name = species_names.get(species_idx, f"Species {species_idx}")
    # Get a unique color that hasn't been used yet
    color = get_unique_color(i, used_colors)
    used_colors.append(color)
    species_colors[species_idx] = color  # Store color by species index for matching with cloud data
    species_colors[species_name.lower()] = color  # Also store by species name for FastChem matching
    label = f"{species_name} ({avg_vmr:.2e})"
    ax1.loglog(vmr[species_idx], pressure_bar, linestyle='solid', color=color, linewidth=2, label=label)

# Plot cloud molecular abundance data if available
if cloud_molec_abundance is not None and cloud_pressure_bar is not None and np.any(cloud_molec_abundance > 0):
    # Use the same color as water vapor (species 7) if available
    cloud_color = species_colors.get(7, 'r')  # Default to red if water color not found
    
    # Plot condensed molecules with dotted line
    ax1.loglog(cloud_molec_abundance, cloud_pressure_bar, linestyle='dotted', color=cloud_color, 
              linewidth=3, label=f"H2O Condensed Molecules")
    print("Added condensed H2O molecular abundance to plot")

# Try to add FastChem species to the same plot
try:
    # Read the FastChem equilibrium chemistry file - use the path defined earlier
    # fastchem_file = f'../output/fastchem/{case_number:05d}/EQ.txt'  # Commenting out incorrect path
    # Use the path defined at the top of the file
    print(f"Attempting to read FastChem data from: {fastchem_file}")
    
    # Read the data
    with open(fastchem_file, 'r') as f:
        lines = f.readlines()
    
    # Find the header line
    header = None
    for i, line in enumerate(lines):
        if 'Pbar' in line and not line.startswith('#'):
            header = line.strip().split()
            break
    
    if header is not None:
        # Get pressure column index
        pbar_idx = header.index('Pbar')
        
        # Skip non-species columns
        skip_columns = ['Pbar', 'Tk', 'n_<tot>', 'n_g', 'mu']
        species_indices = {header.index(col): col for col in header if col not in skip_columns}
        
        # Read data
        fastchem_pressure = []
        species_data = {species: [] for species in species_indices.values()}
        
        for line in lines[i+1:]:  # Start after header
            if line.strip() and not line.startswith('#'):
                values = line.strip().split()
                if len(values) > max(species_indices.keys()):
                    fastchem_pressure.append(float(values[pbar_idx]))
                    for idx, species in species_indices.items():
                        species_data[species].append(float(values[idx]))
        
        # Convert to numpy arrays
        fastchem_pressure = np.array(fastchem_pressure)
        for species in species_data:
            species_data[species] = np.array(species_data[species])
        
        # Find top 15 species by maximum abundance
        species_max_abundance = {species: np.max(species_data[species]) for species in species_data}
        top_fastchem_species = sorted(species_max_abundance.keys(), key=lambda x: species_max_abundance[x], reverse=True)[:15]
        
        print(f"Found {len(top_fastchem_species)} FastChem species to plot")
        
        # Plot top FastChem species
        for i, species in enumerate(top_fastchem_species):
            # Try to match color with EPACRIS species
            species_lower = species.lower()
            if species_lower in species_colors:
                color = species_colors[species_lower]
            else:
                # Get a unique color that hasn't been used yet
                color = get_unique_color(species_lower, used_colors)
                used_colors.append(color)
            
            avg_vmr = species_max_abundance[species]
            ax1.loglog(species_data[species], fastchem_pressure, linestyle='dashed', color=color, 
                     linewidth=3, label=f"{species} (FastChem, {avg_vmr:.2e})")
        
except Exception as e:
    print(f"Warning: Could not plot FastChem species: {e}")
    print(f"Attempted path was: {fastchem_file}")

# Set up the VMR plot
ax1.invert_yaxis()  # Invert y-axis to have pressure increasing downward
ax1.set_xlabel('Volume Mixing Ratio', fontsize=14)
ax1.set_ylabel('Pressure [bar]', fontsize=14)
ax1.set_title('Volume Mixing Ratios (solid: gas, dotted: condensed molecules, dashed: FastChem)', fontsize=16)

# Set axis limits for VMR plot
ax1.set_xlim(1e-14, 1)
ax1.set_ylim(max(pressure_bar), min(pressure_bar))

# Add legend with reasonable size and position
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)

# Right panel - Temperature profile
# Plot EPACRIS profile
if temp_profile is not None:
    temp_profile_p = 10**temp_profile['LogP'] / 1e5  # Convert from Pa to bar
    ax2.plot(temp_profile['Temperature'], temp_profile_p, 'k-', linewidth=2, label='EPACRIS')

# Plot HELIOS profile
if helios_profile is not None:
    ax2.plot(helios_profile['Temp'], helios_profile['Pressure_bar'], 'r--', 
             linewidth=2, label='HELIOS')

if temp_profile is not None or helios_profile is not None:
    ax2.set_yscale('log')  # Set y-axis to logarithmic scale
    ax2.set_xlabel('Temperature [K]', fontsize=14)
    ax2.set_ylabel('Pressure [bar]', fontsize=14)
    ax2.set_title('Temperature Profile', fontsize=16)
    ax2.legend(fontsize=12)
    
    if helios_profile is not None and temp_profile is not None:
        max_temp = max(helios_profile['Temp'].max(), temp_profile['Temperature'].max())
        min_temp = min(helios_profile['Temp'].min(), temp_profile['Temperature'].min())
    elif helios_profile is not None:
        max_temp = helios_profile['Temp'].max()
        min_temp = helios_profile['Temp'].min()
    else:
        max_temp = temp_profile['Temperature'].max()
        min_temp = temp_profile['Temperature'].min()
    
    ax2.set_xlim(min_temp * 0.9, max_temp * 1.1)
    ax2.set_ylim(max(pressure_bar), min(pressure_bar))
    ax2.grid(True, which="both", ls="-", alpha=0.2)
else:
    ax2.text(0.5, 0.5, 'Temperature profiles not available', 
             ha='center', va='center', transform=ax2.transAxes)

# Ensure y-axis range is consistent if FastChem data exists
if 'fastchem_pressure' in locals():
    y_min = min(min(pressure_bar), min(fastchem_pressure))
    y_max = max(max(pressure_bar), max(fastchem_pressure))
    ax1.set_ylim(y_max, y_min)
    ax2.set_ylim(y_max, y_min)

plt.tight_layout()
# Create output directory if it doesn't exist
os.makedirs('output', exist_ok=True)
plt.savefig(f'output/{case_name}_chemistry_and_temp.png', dpi=300, bbox_inches='tight')
plt.close()
