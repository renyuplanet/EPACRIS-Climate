import numpy as np
import matplotlib.pyplot as plt

star_name = 'TOI-199'

# Read the converted Alpha A spectrum
alpha_a_data = np.loadtxt(f'{star_name}_spectra_converted.dat')
alpha_a_wl = alpha_a_data[:, 0]  # nm
alpha_a_flux = alpha_a_data[:, 1]  # W/m²/nm

# Read the solar spectrum
solar_data = np.loadtxt('../solar.txt')
solar_wl = solar_data[:, 0]  # nm
solar_flux = solar_data[:, 1]  # W/m²/nm

# Create the plot
plt.figure(figsize=(12, 6))

# Plot both spectra
plt.semilogy(solar_wl, solar_flux, 'b-', label='Solar Spectrum', alpha=0.7)
plt.semilogy(alpha_a_wl, alpha_a_flux, 'r-', label=f'{star_name} Spectrum', alpha=0.7)

# Add labels and title
plt.xlabel('Wavelength (nm)')
plt.ylabel('Flux (W/m²/nm)')
plt.title(f'Comparison of Solar and {star_name} Spectra at 1 AU')
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.legend()

# Set reasonable axis limits
plt.xlim(100, 2500)  # Focus on UV through near-IR
plt.ylim(1e-7, 2.5)  # Cover the main flux range

# Save the plot
plt.savefig(f'{star_name}_spectrum_comparison.png', dpi=300, bbox_inches='tight')
plt.close()

# Also create a zoomed plot of the visible range
plt.figure(figsize=(12, 6))

# Plot both spectra
plt.plot(solar_wl, solar_flux, 'b-', label='Solar Spectrum', alpha=0.7)
plt.plot(alpha_a_wl, alpha_a_flux, 'r-', label=f'{star_name} Spectrum', alpha=0.7)

# Add labels and title
plt.xlabel('Wavelength (nm)')
plt.ylabel('Flux (W/m²/nm)')
plt.title(f'Comparison of Solar and {star_name} Spectra - Visible Range')
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.legend()

# Set visible range limits
plt.xlim(380, 750)  # Visible light range
plt.ylim(0, 2.5)

# Save the zoomed plot
plt.savefig('spectrum_comparison_visible.png', dpi=300, bbox_inches='tight')
plt.close()

# Print some statistics
print("\nSpectrum Statistics:")
print("-" * 50)
print(f"Total irradiance ({star_name}): {np.trapz(alpha_a_flux, alpha_a_wl):.1f} W/m²")
print(f"Total irradiance (Solar): {np.trapz(solar_flux, solar_wl):.1f} W/m²")
print(f"Peak flux ({star_name}): {np.max(alpha_a_flux):.2f} W/m²/nm")
print(f"Peak flux (Solar): {np.max(solar_flux):.2f} W/m²/nm")
print(f"Peak wavelength ({star_name}): {alpha_a_wl[np.argmax(alpha_a_flux)]:.0f} nm")
print(f"Peak wavelength (Solar): {solar_wl[np.argmax(solar_flux)]:.0f} nm") 