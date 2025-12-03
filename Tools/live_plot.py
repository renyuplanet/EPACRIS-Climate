import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import pandas as pd

# Live plotting script for EPACRIS atmospheric modeling results
# 
# Features:
# - Plots top 10 species by average VMR
# - Can force specific species to be plotted (see force_plot list below)
# - Plots temperature profiles, gas abundances, cloud data, and particle sizes
# - Handles saturation ratios and heat capacity data
#
# Usage: python live_plot.py <step_count> <live_plot_dir> <nmax_iteration> <species_file>

# Get step count, output directory, NMAX iteration, and species file from command line arguments
if len(sys.argv) != 5:
    print("Usage: python live_plot.py <step_count> <live_plot_dir> <nmax_iteration> <species_file>")
    sys.exit(1)
    
total_step_count = int(sys.argv[1])
live_plot_dir = sys.argv[2]
nmax_iteration = int(sys.argv[3])
species_file = sys.argv[4]  # Path to species mapping file from config
species_names = {}
if os.path.exists(species_file):
    with open(species_file, 'r') as f:
        f.readline()  # skip header
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    std_num = int(parts[2])
                    name = parts[0]
                    species_names[std_num] = name
                except Exception:
                    continue

# Define colors for cycling through - using 10 different colors
COLORS = ['k', 'r', 'g', 'b', 'm', 'c', 'y', 'orange', 'purple', 'brown']

# Load data
data = np.loadtxt('Tools/temp_data.txt', skiprows=1)  # Skip header row
t_new = data[:, 0]
pressure = data[:, 1]  # Pressure in bar

# Read radiative diagnostics from file comments
radiative_diagnostics = {}
convergence_tolerances = {}
with open('Tools/temp_data.txt', 'r') as f:
    for line in f:
        if line.startswith('# RADIATIVE_DIAGNOSTICS:'):
            # Updated to parse the new radiation flux values
            parts = line.split()
            for part in parts:
                if '=' in part:
                    key, value = part.split('=')
                    try:
                        radiative_diagnostics[key] = float(value)
                    except ValueError:
                        pass
        elif line.startswith('# CONVERGENCE_TOLERANCES:'):
            parts = line.split()
            for part in parts:
                if '=' in part:
                    key, value = part.split('=')
                    try:
                        convergence_tolerances[key] = float(value)
                    except ValueError:
                        pass

# Get species numbers from header
with open('Tools/temp_data.txt', 'r') as f:
    header = f.readline().strip('# \n').split()

# Parse column headers for different data types
gas_indices = []
cloud_indices = []
gas_species_numbers = []
cloud_species_numbers = []
saturation_ratio_indices = []
saturation_ratio_species = []
particle_size_indices = []
particle_size_species = []
cp_index = None
lapse_index = None
isconv_index = None
pot_temp_index = None

for i, col_name in enumerate(header[2:], 2):  # Skip Temperature and Pressure
    if col_name == 'cp':
        # This is the heat capacity column
        cp_index = i
    elif col_name == 'lapse':
        # This is the lapse rate column
        lapse_index = i
    elif col_name == 'isconv':
        # This is the convective layer flag column
        isconv_index = i
    elif col_name == 'pot_temp':
        # This is the potential temperature column
        pot_temp_index = i
    elif col_name.startswith('sat_ratio_'):
        # This is a saturation ratio column (e.g., sat_ratio_7, sat_ratio_9)
        saturation_ratio_indices.append(i)
        species_num = int(col_name.split('_')[2])  # Extract species number from sat_ratio_X
        saturation_ratio_species.append(species_num)
    elif col_name.startswith('particle_number_density_'):
        # This is a particle number density column (e.g., particle_number_density_7, particle_number_density_9)
        particle_size_indices.append(i)
        species_num = int(col_name.split('_')[3])  # Extract species number from particle_number_density_X
        particle_size_species.append(species_num)
    elif col_name.startswith('c'):
        # This is a cloud species (with 'c' prefix)
        cloud_indices.append(i)
        cloud_species_numbers.append(int(col_name[1:]))  # Remove 'c' prefix to get species number
    else:
        # This is a gas species
        gas_indices.append(i)
        gas_species_numbers.append(int(col_name))

# Extract gas, cloud, heat capacity, lapse rate, potential temperature, saturation ratio, particle size, and convective data
gas_data = data[:, gas_indices]
cloud_data = data[:, cloud_indices]
if saturation_ratio_indices:
    saturation_ratio_data = data[:, saturation_ratio_indices]
else:
    saturation_ratio_data = None

if particle_size_indices:
    particle_size_data = data[:, particle_size_indices]
else:
    particle_size_data = None

if cp_index is not None:
    cp_data = data[:, cp_index]

if lapse_index is not None:
    lapse_data = data[:, lapse_index]

if pot_temp_index is not None:
    pot_temp_data = data[:, pot_temp_index]

if isconv_index is not None:
    isconv_data = data[:, isconv_index].astype(int)
    num_convective_layers = np.sum(isconv_data)
else:
    isconv_data = None
    num_convective_layers = 0

# Calculate average VMR for each gas species (excluding zeros and NaNs)
gas_avg_vmr = {}
for i, col_idx in enumerate(gas_indices):
    vmr = data[:, col_idx]
    vmr = vmr[~np.isnan(vmr)]  # Remove NaNs
    vmr = vmr[vmr > 0]         # Remove zeros
    if len(vmr) > 0:
        gas_avg_vmr[i] = np.mean(vmr)

# Get top 10 gas species by average VMR
top_10_gas = sorted(gas_avg_vmr.items(), key=lambda x: x[1], reverse=True)[:10]

# Create figure with four subplots - add saturation ratios
fig = plt.figure(figsize=(26, 8))  # Wider to accommodate fourth panel
gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1], wspace=0.2)  
ax1 = fig.add_subplot(gs[0])  # Temperature plot
ax2 = fig.add_subplot(gs[1])  # Abundance plot
ax3 = fig.add_subplot(gs[2])  # Heat capacity plot
ax4 = fig.add_subplot(gs[3])  # Saturation ratio plot

# Temperature plot
# First, plot the complete continuous temperature profile
ax1.plot(t_new, pressure, 'k-', linewidth=3, label='Temperature Profile')

# Add potential temperature if available
if pot_temp_index is not None:
    # Filter out zero/very small potential temperature values (which occur in radiative layers)
    pot_temp_plot = np.copy(pot_temp_data)
    pot_temp_plot[pot_temp_plot < 1.0] = np.nan  # Set very small values to NaN for cleaner plotting
    ax1.plot(pot_temp_plot, pressure, 'g--', linewidth=2, alpha=0.8, label='Potential Temperature')

# Add H2O condensation curve for comparison
def h2o_saturation_pressure(temp):
    """Calculate H2O saturation pressure in Pa using same formulation as EPACRIS"""
    if temp < 273.16:
        # Over ice (Murphy & Koop 2005)
        psat = np.exp(9.550426 - 5723.265/temp + 3.53068*np.log(temp) - 0.00728332*temp)
    else:
        # Over water (Seinfeld & Pandis 2006)
        a = 1 - 373.15/temp
        psat = 101325 * np.exp(13.3185*a - 1.97*a*a - 0.6445*a*a*a - 0.1229*a*a*a*a)
    return psat

# Find H2O data in the gas species (species 7 is H2O)
h2o_vmr_data = None
h2o_gas_index = None
for i, species_num in enumerate(gas_species_numbers):
    if species_num == 7:  # H2O is species 7
        h2o_vmr_data = gas_data[:, i]
        h2o_gas_index = i
        break

# Calculate H2O condensation curve using ACTUAL H2O abundances from simulation
if h2o_vmr_data is not None:
    condensation_temps = []
    valid_pressures = []
    
    for p_bar, vmr in zip(pressure, h2o_vmr_data):
        if vmr > 1e-15:  # Only calculate for layers with significant H2O
            p_total = p_bar * 1e5  # Convert bar to Pa
            p_partial = vmr * p_total  # Actual H2O partial pressure at this layer
            
            # Simple approach: solve p_partial = p_sat(T) for T
            temp_guess = 200.0
            for _ in range(1000):  # Simple iteration
                p_sat = h2o_saturation_pressure(temp_guess)
                if abs(p_partial - p_sat) < 1e-6:
                    break
                # Simple adjustment
                if p_partial > p_sat:
                    temp_guess += 1.0
                else:
                    temp_guess -= 1.0
                temp_guess = max(50.0, min(2000.0, temp_guess))
            
            condensation_temps.append(temp_guess)
            valid_pressures.append(p_bar)
    
    # Plot the condensation curve
    if len(condensation_temps) > 0:
        ax1.plot(condensation_temps, valid_pressures, 'cyan', linewidth=3, linestyle='-', 
                 label=f'H2O Condensation', alpha=0.9,zorder=100)
    

# Mark convective layers if data is available
if isconv_data is not None:
    # Create mask for convective layers
    convective_mask = isconv_data == 1
    
    # Overlay convective layers with red markers and thicker segments
    if np.any(convective_mask):
        # Plot red circles on convective layers
        ax1.scatter(t_new[convective_mask], pressure[convective_mask], 
                   c='red', s=40, alpha=0.4, zorder=1, marker='o',
                   label=f'Convective layers ({num_convective_layers})')
        
        # Optionally add a thicker red line segment over convective parts
        # This creates a visual emphasis without breaking continuity
        # for i in range(len(t_new)):
        #     if convective_mask[i]:
        #         # Draw short red line segments for convective layers
        #         if i > 0:
        #             ax1.plot([t_new[i-1], t_new[i]], [pressure[i-1], pressure[i]], 
        #                     'r-', linewidth=6, alpha=0.6, zorder=3)

# Add text showing number of convective layers
# ax1.text(0.02, 0.98, f'Convective layers: {num_convective_layers}', 
#          transform=ax1.transAxes, fontsize=12, fontweight='bold',
#          verticalalignment='top', horizontalalignment='left',
#          bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Add radiative diagnostics box
if radiative_diagnostics and convergence_tolerances:
    # Calculate convergence status
    rflux_max = radiative_diagnostics.get('Rfluxmax', 0.0)
    drflux_max = radiative_diagnostics.get('dRfluxmax', 0.0)
    rt_step = radiative_diagnostics.get('RTstep', 0)
    equilib_layers = radiative_diagnostics.get('equilib_layers', 0)
    ncl = radiative_diagnostics.get('ncl', 0)
    nrl = radiative_diagnostics.get('nrl', 0)
    
    tol_rc = convergence_tolerances.get('Tol_RC', 1e-3)
    tol_rc_r = convergence_tolerances.get('Tol_RC_R', 1e-3)
    bb_flux = convergence_tolerances.get('BBflux', 1.0)
    
    # Calculate convergence ratios
    flux_ratio = rflux_max / tol_rc if tol_rc > 0 else float('inf')
    bb_flux_ratio = rflux_max / (tol_rc_r * bb_flux) if (tol_rc_r * bb_flux) > 0 else float('inf')
    gradient_ratio = drflux_max / (0.2 * tol_rc) if tol_rc > 0 else float('inf')
    
    # Determine convergence status
    converged = (flux_ratio < 1.0) or (bb_flux_ratio < 1.0) or (gradient_ratio < 1.0)
    status_color = 'lightgreen' if converged else 'lightcoral'
    status_text = 'CONVERGED' if converged else 'ITERATING'
    
    # Enhanced climate balance calculations
    # Calculate global energy balance indicators
    stefan_boltzmann = 5.67e-8  # W m^-2 K^-4
    
    # Surface temperature (bottom of atmosphere)
    t_surface = t_new[-1] if len(t_new) > 0 else 0.0
    
    # TOA temperature (top of atmosphere) 
    t_toa = t_new[0] if len(t_new) > 0 else 0.0
    
    # Estimate outgoing longwave radiation from surface temperature
    olr_surface = stefan_boltzmann * t_surface**4 if t_surface > 0 else 0.0
    
    # Estimate effective TOA radiation from TOA temperature  
    olr_toa = stefan_boltzmann * t_toa**4 if t_toa > 0 else 0.0
    
    # Get actual stellar radiation from EPACRIS
    incoming_flux = radiative_diagnostics.get('radiationI0', 0.0) # TOA incoming flux [W/m²]
    boa_incoming_flux = radiative_diagnostics.get('radiationI1', 0.0) # BOA incoming flux [W/m²]  
    toa_net_flux = radiative_diagnostics.get('radiationO', 0.0) # TOA NET flux (upward-downward) [W/m²]
    # IMPORTANT: radiationO = NetFlux[0][0] = upward_flux - downward_flux
    # Negative values = planet gaining energy, Positive values = planet losing energy
    
    # Energy balance percentage (how close are we to equilibrium?)
    # Perfect equilibrium means net flux = 0
    # Use actual TOA net flux if available
    if 'radiationO' in radiative_diagnostics and incoming_flux > 0.0:
        # Energy balance based on how close net flux is to zero relative to incoming flux
        energy_balance_pct = (1.0 - abs(toa_net_flux) / incoming_flux) * 100
    else:
        # Fall back to using flux residual vs blackbody flux
        energy_balance_pct = (1.0 - abs(rflux_max) / bb_flux) * 100 if bb_flux > 0 else 0.0
    energy_balance_pct = max(0.0, min(100.0, energy_balance_pct))  # Clamp to 0-100%
    
    # Temperature stability (difference between current surface/TOA temps)
    temp_gradient = abs(t_surface - t_toa) if (t_surface > 0 and t_toa > 0) else 0.0
    
    # Format diagnostic information 
    diag_text = f"Energy Balance: {energy_balance_pct:.1f}%\n"
    
    # Status based on energy balance
    if energy_balance_pct >= 99:
        status_text = "EXCELLENT"
        status_color = 'green'
    elif energy_balance_pct >= 85:
        status_text = "GOOD"
        status_color = 'darkgreen'
    elif energy_balance_pct >= 70:
        status_text = "FAIR"
        status_color = 'orange'
    else:
        status_text = "POOR"
        status_color = 'red'
    
    diag_text += f"Status: {status_text}\n\n"
    
    # Enhanced radiation budget section using actual values
    diag_text += "Radiation Budget (W/m²):\n"
    diag_text += f"TOA Incoming: {incoming_flux:.1e}\n"
    if boa_incoming_flux != 0.0:
        diag_text += f"BOA Incoming: {boa_incoming_flux:.1e}\n"
    if 'radiationO' in radiative_diagnostics:
        # Calculate outgoing flux from net flux
        toa_outgoing_flux = incoming_flux - toa_net_flux
        diag_text += f"TOA Outgoing: {toa_outgoing_flux:.1e}\n"
        diag_text += f"TOA Net (In-Out): {toa_net_flux:.1e}\n"
        # Physical interpretation
        if toa_net_flux > 1e-6:  # Small threshold for floating point precision
            diag_text += f"Status: Planet losing energy\n"
        elif toa_net_flux < -1e-6:
            diag_text += f"Status: Planet gaining energy\n"
        else:
            diag_text += f"Status: Near equilibrium\n"
    else:
        # Fallback to surface OLR estimates
        diag_text += f"Surface OLR: {olr_surface:.1e}\n"
        diag_text += f"TOA OLR: {olr_toa:.1e}\n"
        net_balance = incoming_flux - olr_toa
        diag_text += f"Net Balance: {net_balance:.1e}\n"
    
    diag_text += "\n"
    
    # Enhanced error metrics section 
    diag_text += f"Error Metrics:\n"
    diag_text += f"Max flux error: {rflux_max:.2e} W/m²\n"
    diag_text += f"Max flux gradient: {drflux_max:.2e}\n"
    
    # Convergence info
    diag_text += f"Convergence Info:\n"
    diag_text += f"RT Step: {int(rt_step):d}\n"
    diag_text += f"Equilibrium layers: {int(equilib_layers):d}/{int(equilib_layers + ncl):d}\n"
    
    # Thermal structure 
    diag_text += f"Thermal Structure:\n"
    diag_text += f"Surface temp: {t_surface:.1f} K\n"
    diag_text += f"TOA temp: {t_toa:.1f} K\n"
    temp_gradient = abs(t_surface - t_toa) if (t_surface > 0 and t_toa > 0) else 0.0
    diag_text += f"Temp gradient: {temp_gradient:.1f} K\n\n"
    
    # Convection summary
    diag_text += f"Convection:\n"
    diag_text += f"Convective: {int(ncl):d} layers\n"
    diag_text += f"Radiative: {int(nrl):d} layers"
    
    ax1.text(0.98, 0.98, diag_text, 
             transform=ax1.transAxes, fontsize=9, fontweight='normal',
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8),
             color=status_color)

ax1.set_yscale('log')
ax1.set_ylim(max(pressure), min(pressure))
ax1.set_xlim(0,1500)
ax1.set_xlabel('Temperature (K)', fontsize=12)
ax1.set_ylabel('Pressure (Bar)', fontsize=12)
ax1.set_title(f'Temperature Profile (Step {total_step_count})', fontsize=14)
ax1.tick_params(axis='both', which='major', labelsize=10)
ax1.legend(fontsize=9, loc='upper left')
ax1.grid(True, alpha=0.3)

# Create a dictionary to map species to colors for consistent coloring
species_colors = {}

# Abundance plot - top 10 gas species
for i, (idx, avg) in enumerate(top_10_gas):
    color = COLORS[i % len(COLORS)]
    vmr = gas_data[:, idx]
    std_num = gas_species_numbers[idx]
    name = species_names.get(std_num, f"Species {std_num}")
    label = f"{name} ({avg:.2e})"
    ax2.plot(vmr, pressure, color=color, linewidth=2, label=label)
    
    # Store color for this species
    species_colors[std_num] = color

# Get list of species already plotted in top 10 to avoid duplicates
already_plotted = set(gas_species_numbers[idx] for idx, _ in top_10_gas)

# Force specific species to be plotted (in addition to top 10)
force_plot = []
for std_num_str in force_plot:
    try:
        std_num = int(std_num_str)
        # Only plot if not already plotted in top 10
        if std_num not in already_plotted:
            # Find the index of the species in gas_species_numbers
            try:
                species_idx = gas_species_numbers.index(std_num)
                vmr = gas_data[:, species_idx]
                
                # Calculate average VMR for the label (excluding NaNs and very small values)
                vmr_clean = vmr[~np.isnan(vmr)]  # Remove NaNs
                vmr_clean = vmr_clean[vmr_clean > 1e-30]  # Remove very small values (more appropriate threshold)
                if len(vmr_clean) > 0:
                    avg_vmr = np.mean(vmr_clean)
                    label = f"{name} (forced, {avg_vmr:.2e})"
                else:
                    label = f"{name} (forced, <1e-30)"
                
                # Use a distinct color for forced species (avoiding index errors)
                color_idx = (len(COLORS) + force_plot.index(std_num_str)) % len(COLORS)
                color = COLORS[color_idx]
                ax2.plot(vmr, pressure, color=color, linewidth=2, label=label)
                # Store color for this species
                species_colors[std_num] = color
            except ValueError:
                print(f"Warning: Species number {std_num} not found in gas_species_numbers for forced plot.")
        else:
            print(f"Info: Species {std_num} already plotted in top 10, skipping duplicate.")
    except ValueError:
        print(f"Warning: Skipping invalid species number '{std_num_str}' in force_plot list.")


# Check for and plot non-zero cloud data
has_clouds = False
for i, std_num in enumerate(cloud_species_numbers):
    cloud_vmr = cloud_data[:, i]
    
    # Check if cloud has any significant values (above threshold)
    if np.any(cloud_vmr > 1e-10):
        has_clouds = True
        # Use same color as gas species if available
        color = species_colors.get(std_num, COLORS[i % len(COLORS)])
        name = species_names.get(std_num, f"Species {std_num}")
        
        # Get non-zero values for proper averaging
        significant_values = cloud_vmr[cloud_vmr > 1e-10]
        avg_vmr = np.mean(significant_values)
        
        # Plot clouds with dashed lines and different marker
        ax2.plot(
            cloud_vmr, 
            pressure, 
            color=color, 
            linestyle='--',
            marker='o', 
            markersize=0,
            linewidth=2,
            label=f"{name} (cloud, {avg_vmr:.2e})"
        )


ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylim(max(pressure), min(pressure))
ax2.set_xlim(1e-14, 1)
ax2.set_xlabel('Volume Mixing Ratio', fontsize=12)
ax2.set_ylabel('Pressure (Bar)', fontsize=12)
ax2.set_title('Gas & Cloud Abundances + Cloud Number Density', fontsize=14)
ax2.tick_params(axis='both', which='major', labelsize=10)
ax2.grid(True, alpha=0.3)

# Create twin axis for particle number density
if particle_size_data is not None:
    ax2_twin = ax2.twiny()
    
    # Plot particle number density for condensible species
    for i, species_num in enumerate(particle_size_species):
        particle_number_density = particle_size_data[:, i]
        
        # Only plot if there are significant number densities
        if np.any(particle_number_density > 1e0):  # Threshold: 1 particle/m³
            # Use same color as corresponding species if available
            color = species_colors.get(species_num, COLORS[i % len(COLORS)])
            name = species_names.get(species_num, f"Species {species_num}")
            
            # Get non-zero values for proper averaging
            significant_density = particle_number_density[particle_number_density > 1e0]
            avg_density = np.mean(significant_density) if len(significant_density) > 0 else 0
            
            # Plot particle number density with dotted lines and triangle markers
            ax2_twin.plot(
                particle_number_density, 
                pressure, 
                color='red', 
                linestyle='-',
                linewidth=2,
                alpha=1,
                label=f"{name} (n={avg_density:.1e} m⁻³)"
            )
    
    ax2_twin.set_xscale('log')
    ax2_twin.set_yscale('log')
    ax2_twin.set_ylim(max(pressure), min(pressure))
    ax2_twin.set_xlim(1e0, 1e15)  # Particle number density range: 1 to 1e15 particles/m³
    ax2_twin.set_xlabel('Particle Number Density (m⁻³)', fontsize=12, color='darkred')
    ax2_twin.tick_params(axis='x', colors='darkred')
    
    # Combined legend for both axes
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_twin.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, bbox_to_anchor=(0.3, 0.91), loc='center', fontsize=7, ncol=2, framealpha=0.2)
else:
    ax2.legend(bbox_to_anchor=(0.3, 0.91), loc='center', fontsize=8, ncol=2, framealpha=0.2)

# Heat capacity and lapse rate plot with twin axes
ax3.plot(cp_data, pressure, 'r-', linewidth=3, label='Heat Capacity')
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_ylim(max(pressure), min(pressure))
# Set reasonable x-limits for heat capacity (typical range for atmospheric gases)
ax3.set_xlim(1e-1, 500)
ax3.set_xlabel('Heat Capacity (J mol⁻¹ K⁻¹)', fontsize=12, color='red')
ax3.set_ylabel('Pressure (Bar)', fontsize=12)
ax3.set_title('Heat Capacity & Lapse Rate', fontsize=14)
ax3.tick_params(axis='both', which='major', labelsize=10)
ax3.tick_params(axis='x', colors='red')
ax3.grid(True, alpha=0.3)

# Create twin axis for lapse rate
if lapse_index is not None:
    ax3_twin = ax3.twiny()
    ax3_twin.plot(lapse_data, pressure, 'b-', linewidth=3, label='Lapse Rate')
    ax3_twin.set_yscale('log')
    ax3_twin.set_ylim(max(pressure), min(pressure))
    # Set reasonable x-limits for lapse rate (typical values 0.1 to 1.0)
    ax3_twin.set_xlim(0.1, 0.4)
    ax3_twin.set_xlabel('Lapse Rate (d ln T / d ln P)', fontsize=12, color='blue')
    ax3_twin.tick_params(axis='x', colors='blue')
    
    # Combined legend
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3_twin.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2, fontsize=10, loc='upper right')
else:
    ax3.legend(fontsize=10)

# Saturation ratio plot
if saturation_ratio_data is not None and len(saturation_ratio_indices) > 0:
    # Species names mapping for condensibles
    condensible_names = {7: 'H₂O', 9: 'NH₃', 20: 'CO', 21: 'CH₄', 52: 'CO₂', 
                        53: 'H₂', 54: 'O₂', 55: 'N₂'}
    
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
    
    for i, species_num in enumerate(saturation_ratio_species):
        sat_ratio = saturation_ratio_data[:, i]
        color = colors[i % len(colors)]
        name = condensible_names.get(species_num, f"Species {species_num}")
        
        ax4.plot(sat_ratio, pressure, color=color, linewidth=2, label=f"{name}")
    
    # Add vertical line at saturation ratio = 1.0 (equilibrium)
    ax4.axvline(x=1.0, color='black', linestyle='--', linewidth=1, alpha=0.7, label='Saturation (S=1)')
    
    # Add shaded regions to highlight supersaturation (S>1) and subsaturation (S<1)
    ax4.axvspan(1.0, 10, alpha=0.1, color='red', label='Supersaturated (S>1)')
    ax4.axvspan(0.01, 1.0, alpha=0.1, color='blue', label='Subsaturated (S<1)')
    
    #ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_ylim(max(pressure), min(pressure))
    ax4.set_xlim(0.6, 1.4)  # Reasonable range for saturation ratios
    ax4.set_xlabel('Saturation Ratio (P/Psat)', fontsize=12)
    ax4.set_ylabel('Pressure (Bar)', fontsize=12)
    ax4.set_title('Saturation Ratios', fontsize=14)
    ax4.tick_params(axis='both', which='major', labelsize=10)
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3)
else:
    # If no saturation ratio data, show placeholder
    ax4.text(0.5, 0.5, 'No Saturation\nRatio Data\nAvailable', 
             transform=ax4.transAxes, fontsize=14, ha='center', va='center',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))
    ax4.set_title('Saturation Ratios', fontsize=14)
    ax4.set_xlabel('Saturation Ratio (P/Psat)', fontsize=12)
    ax4.set_ylabel('Pressure (Bar)', fontsize=12)

plt.savefig(f'{live_plot_dir}TP_profile_{nmax_iteration}_step_{total_step_count}.png', dpi=100)
plt.close()
