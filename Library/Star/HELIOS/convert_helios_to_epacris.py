import numpy as np
import os
from astropy import constants as const
from pathlib import Path

def validate_spectrum(wavelength_nm, flux_wm2_nm, reference_file=None):
    """Validate the converted spectrum by checking:
    1. Wavelength range should start around 100 nm (or lower for EUV)
    2. Flux values should be in reasonable range for stellar spectra at 1 AU
    3. Compare to reference spectrum if provided
    4. Check total irradiance against solar constant (~1366.1 W/m²)
    """
    # Basic checks
    print("\nValidation Results:")
    print("-" * 50)
    
    # Check wavelength range and spacing
    wl_min, wl_max = wavelength_nm.min(), wavelength_nm.max()
    print(f"Wavelength range: {wl_min:.1f} - {wl_max:.1f} nm")
    if wl_min > 100:
        print("WARNING: Spectrum starts above 100 nm - missing UV/EUV range!")
    if wl_max < 2400:
        print("WARNING: Spectrum ends before IR range!")
        
    # Check wavelength spacing is monotonic and reasonable
    wl_spacing = np.diff(wavelength_nm)
    if not np.all(wl_spacing > 0):
        print("ERROR: Wavelength grid is not strictly increasing!")
    
    # Check flux values
    flux_min, flux_max = flux_wm2_nm.min(), flux_wm2_nm.max()
    print(f"Flux range: {flux_min:.2e} - {flux_max:.2e} W/m²/nm")
    if np.any(flux_wm2_nm < 0):
        print("ERROR: Negative flux values detected!")
    if flux_max > 2.5:  # Based on typical solar spectrum peak
        print("WARNING: Maximum flux seems unusually high for 1 AU!")
    
    # Calculate total irradiance (integrate over wavelength)
    total_irr = np.trapz(flux_wm2_nm, wavelength_nm)
    print(f"Total irradiance at 1 AU: {total_irr:.2e} W/m²")
    # Compare to solar constant
    if total_irr > 2000:
        print("WARNING: Total irradiance exceeds expected range (~1366.1 W/m² for Sun)")
    
    # Compare to reference if provided
    if reference_file and os.path.exists(reference_file):
        ref_data = np.loadtxt(reference_file)
        ref_wl, ref_flux = ref_data[:, 0], ref_data[:, 1]
        # Check if wavelength ranges overlap
        wl_overlap = (wl_min <= ref_wl.max()) and (wl_max >= ref_wl.min())
        if not wl_overlap:
            print("WARNING: No wavelength overlap with reference spectrum!")
        else:
            # Interpolate reference to same wavelength grid
            ref_flux_interp = np.interp(wavelength_nm, ref_wl, ref_flux)
            rel_diff = np.abs(flux_wm2_nm - ref_flux_interp) / ref_flux_interp
            print(f"Mean relative difference from reference: {np.mean(rel_diff):.2%}")
            
            # Compare total irradiance
            ref_total = np.trapz(ref_flux, ref_wl)
            print(f"Reference total irradiance: {ref_total:.2e} W/m²")
            print(f"Relative difference in total irradiance: {abs(total_irr-ref_total)/ref_total:.2%}")
    
    print("-" * 50)
    return True

# Define constants

AU_TO_M = 1.495978707e8  # 1 AU in km (matching MATLAB value)
RSUN_TO_KM = 695700  # Solar radius in km (matching MATLAB value)

### IMPORTNT SET STAR RADIUS CORRECTLY
STAR_RADIUS = 0.82 * RSUN_TO_KM  # star radius in km

# Get the directory where this script is located
script_dir = Path(__file__).parent.absolute()

# Path to the input file
#input_file = script_dir / "Alpha_A_spectra.dat"
#output_file = script_dir / "Alpha_A_spectra_converted.dat"

input_file = script_dir / "TOI-199_spectra.dat"
output_file = script_dir / "TOI-199_spectra_converted.dat"

# Check if the input file exists
if not input_file.exists():
    print(f"Error: Input file {input_file} not found.")
    exit(1)

# Read the stellar spectrum data
data = np.loadtxt(input_file, comments='#', skiprows=1)

# Extract wavelength and flux
wavelength_um = data[:, 0]  # wavelength already in microns
flux = data[:, 1]        # flux in erg s⁻¹ cm⁻³

# Convert flux to W/m²/μm at stellar surface
# erg/s/cm³ → W/m²/μm:
# - erg → W: ×1e-7
# - path length conversion: ×1e4
# - area conversion: ×1e-4
flux_wm2_um = flux * 1e-7 * 1e4 * 1e-4

# Rebin to R=1000 resolution (same as MATLAB script)
w4r = [0.1]  # Start at 0.1 μm
while w4r[-1] < 200:
    w4r.append(w4r[-1] * (1 + 1/1000))
w4r = np.array(w4r)

# Interpolate flux to new wavelength grid
f4r = np.interp(w4r[:-1], wavelength_um, flux_wm2_um)
w4rr = 0.5 * (w4r[:-1] + w4r[1:])  # Center wavelengths

# Scale to 1 AU
f4rr = f4r * (STAR_RADIUS / AU_TO_M)**2

# Convert to final units:
# - wavelength: μm → nm (×1000)
# - flux: W/m²/μm → W/m²/nm (÷1000)
wavelength_nm = w4rr * 1000
flux_wm2_nm = f4rr / 1000

# Calculate effective temperature from the spectrum for validation
total_flux = np.trapz(f4r, w4rr)  # Total flux at stellar surface
teff_derived = (total_flux / 5.67e-8)**(1/4)  # Derive Teff from total flux
print(f"\nDerived effective temperature: {teff_derived:.0f} K")
print(f"Expected effective temperature: 5255 K")

# Validate the conversion
reference_file = script_dir.parent / "solar.txt"
print(f"Validating against {reference_file}")
validate_spectrum(wavelength_nm, flux_wm2_nm, str(reference_file))

# Save the converted data
np.savetxt(output_file, np.column_stack((wavelength_nm, flux_wm2_nm)), 
           fmt='%.6f\t%.6e', header='Wavelength (nm)    Flux (W/m²/nm at 1 AU)')

print(f"\nConversion complete. Output saved to {output_file}")
