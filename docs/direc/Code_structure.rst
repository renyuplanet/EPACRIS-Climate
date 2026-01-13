Code Structure
==============

This section describes the structure of EPACRIS, the main program flow, and the organization of key routines. The code is written in C and follows a single-file compilation approach where the main file ``epacris_main.c`` includes all auxiliary C files. The program iteratively solves for the atmospheric temperature-pressure profile and chemical composition by coupling radiative transfer, convective adjustment, and chemistry modules until convergence is achieved.

Program Architecture
--------------------

EPACRIS consists of three main computational levels:

1. **Main driver** ``epacris_main.c``: Initializes atmospheric grids, loads opacity and chemistry data, manages climate-chemistry coupling iterations, and coordinates the overall simulation flow.

2. **Climate solver** ``climate.c``: Computes radiative-convective equilibrium by iteratively solving radiative transfer and applying convective adjustments. This module calls cloud physics and condensation routines.

3. **Supporting modules**: Provide specialized functionality including cloud physics ``conv_cond_funcs.c``, cloud optics ``cloud_optics.c``, opacity reading ``readcross.c``, ``readcia.c``, chemistry ``chemequil.c``, and utilities ``Interpolation.c``, ``Convert.c``.

Main Program Flow ``epacris_main.c``
--------------------------------------

The following sections describe each major step in the execution sequence.

Initialization and Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Line 0: Parameter setup and information display**

* Initiliaze global variables and import headers/C files
* Print run configuration information including planet properties (mass, radius, orbital distance), solver settings, and cloud physics mode.
* Initialize loop counters and local variables.
* Calculate global values like planet surface gravity and convert irradiation angle to radians.

Wavelength Grid Setup
~~~~~~~~~~~~~~~~~~~~~~

**Line ~200: Construct wavelength grid**

* Generate log-spaced wavelength array used for all radiative transfer calculations from ``LAMBDALOW`` to ``LAMBDAHIGH`` with ``NLAMBDA`` points.
* Wavelength stored in nanometers for consistency with opacity databases.

Rayleigh Scattering
~~~~~~~~~~~~~~~~~~~

**Line ~205: Calculate Rayleigh cross-sections**

* Compute refractive index for atmospheric gas (H2, N2, CO2, etc.) at each wavelength using functions from ``RefIdx.c``.
* Calculate wavelength-dependent Rayleigh scattering cross-section
* Store in ``crossr[]`` array for use in radiative transfer calculations.

Stellar Spectrum
~~~~~~~~~~~~~~~~

**Line ~220: Load and process stellar spectrum**

* Read stellar spectrum file specified by ``STAR_SPEC``.
* Interpolate stellar flux onto model wavelength grid.
* Scale flux from 1 AU to planet's orbital distance: ``F_planet = F_1AU / ORBIT²``.
* Apply ``FaintSun`` factor to simulate planetary albedo or stellar variability.
* Extrapolate to longer wavelengths using Rayleigh-Jeans tail if needed.

Temperature-Pressure Profile Initialization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~250: Set up initial T-P profile**

Two modes are supported based on ``TPMODE``:

* **TPMODE = 0**: Generate parametric T-P profile using ``TPPara()``
  
  * Defined using pressure boundaries (``PTOP``, ``PBOTTOM``, ``PMIDDLE``, etc.).
  * Calculate equilibrium temperature at top of atmosphere from stellar irradiation.
  * Initializes an isothermal profile with automatic equilibrium temperature if ``TTOP = 0``


* **TPMODE = 1**: Import T-P profile from file specified by ``TPLIST``
  
  * Read altitude (km), pressure (log₁₀(Pa)), and temperature (K) data from file.
  * Interpolate temperature onto model pressure grid.

 
**Altitude Grid Calculation:**

For each layer from bottom to top:

* Calculate scale height
* Set altitude grid: ``z[j] = z[j-1] - H * ln(P[j]/P[j-1])``.
* Layer centers at midpoints: ``zl[j] = (z[j] + z[j-1]) / 2``.

**Double Grid Construction:**

EPACRIS uses a double-resolution grid ``Tdoub``, ``Pdoub``, ``MMdoub``, ``zdoub`` for accurate radiative transfer in non-isothermal layers:

* Even indices ``2*j``: Layer boundaries.
* Odd indices ``2*j-1``: Layer centers.

This allows the radiative transfer solver to account for temperature variations within layers.

Chemistry Setup
~~~~~~~~~~~~~~~

**Line ~400: Load chemical species and reactions**

**Species List Import:**

* Read species file specified by ``SPECIES_LIST`` e.g., ``Library/SpeciesList/species_HNCSO.dat``.
* Parse species properties: name, type (X/F/C/A), standard number, molecular mass, initial mixing ratio, boundary conditions.
* Classify species:
  
  * **X (Variable)**: Species solved with full continuity equation (``numx`` species).
  * **F (Fast)**: Species in photochemical equilibrium (``numf`` species).
  * **C (Constant)**: Fixed mixing ratio species (``numc`` species).
  * **A (Aerosol)**: Aerosol/haze species treated as variable (subset of X).

**Reaction List Import (Not used in this version):**

* Read reaction file specified by ``REACTION_LIST``.

**Photolysis Cross-Sections:**

* For each photolysis reaction, read wavelength-dependent absorption cross-sections and quantum yields.
* Load from ``Library/PhotochemOpa/[species_name]`` files.
* Interpolate onto model wavelength grid.
* Store in ``cross[][]`` and ``qy[][]`` arrays.

**Initial Atmospheric Composition**

**Line ~790: Load chemical species and reactions**

Three initialization modes based on ``IMODE``:

* **IMODE = 0**: Chemical equilibrium
  
  * Call ``chemquil()`` from ``chemequil.c`` to minimize Gibbs free energy.
  * Read elemental abundances from ``ELE_ABUN`` file.
  * Compute equilibrium mixing ratios at each layer given T, P, and elemental budget.
  * Convert mixing ratios to number densities: ``xx[j][i] = mixing_ratio[j][i] * MM[j]``.

* **IMODE = 1**: Use initial values from species list
  
  * Apply mixing ratios specified in ``SPECIES_LIST`` file at all layers.

* **IMODE = 2**: Import from previous run
  
  * Read ``IMODE_CHEM_FILE`` from previous simulation.
  * File contains mixing ratios (not number densities).
  * Interpolate composition onto current pressure grid using nearest-value extrapolation for out-of-range pressures.
  * Convert mixing ratios back to number densities: ``xx[j][i] = mixing_ratio * MM[j]``.

**Mean Molecular Mass Calculation:**

For each layer:

* Sum total number density and total mass from all species.
* Calculate Helium abundance to fill remainder: ``n_He = n_total - Σ n_i``.
* Compute mean molecular mass: ``μ = Σ(n_i * m_i) / n_total``.

Opacity Loading
~~~~~~~~~~~~~~~

**Line ~990: Initialize opacity arrays and load data**

**Memory Allocation:**

* Allocate 2D arrays for each molecular opacity: ``opacH2O[1..zbin][0..NLAMBDA-1]``.
* Allocate cloud optical property arrays for each defined cloud species ``CLOUD_SPECIES_LIST``: ``cH2O``, ``aH2O``, ``gH2O``... 

**Opacity File Reading:**

* Call ``read_all_opacities()`` from ``readcross.c`` to load gas opacities.
* For each species in ``OPACITY_SPECIES_LIST``:
  
  * Open file ``[CROSSHEADING]/opac[species].dat``.
  * Read 2D opacity table (wavelength × temperature-pressure points).
  * Cache opacity data in memory for efficient reinterpolation.
  * Interpolate onto model's wavelength, temperature, and pressure grid.
  * Store in corresponding ``opacXXX[][]`` array.

**CIA Opacity Loading:**

* Call ``readcia()`` from ``readcia.c`` to load collision-induced absorption.
* Read tables for H₂-H₂, H₂-He, N₂-N₂, N₂-H₂, CO₂-CO₂ (**species currently hardcoded**).
* Interpolate onto model grid and store in ``XXXCIA[][]`` arrays.

**Cloud Optical Property Loading:**

* Call ``read_cloud_optical_tables_mie()`` from ``cloud_optics.c`` (see :ref: `cloud_optics.c`).
* Load Mie scattering lookup tables for condensible species (H₂O, NH₃, etc.) specified in ``CLOUD_SPECIES_LIST``.
* Tables contain extinction cross-section, single-scattering albedo, and asymmetry parameter as function of particle size and wavelength.
* Support two formats: EPACRIS format ``USE_EPACRIS_FORMAT=1`` or LX-Mie format ``USE_EPACRIS_FORMAT=0``.

Output File Preparation
~~~~~~~~~~~~~~~~~~~~~~~~

**Lines ~990: Set up output files and copy configuration**

* Create output directory if it doesn't exist (``OUT_DIR``).
* Define output file paths:
  
  * ``ConcentrationSTD_T.dat``: Final atmospheric composition (mixing ratios).
  * ``NewTemperature.dat``: Final temperature-pressure profile.
  * ``Diagnostic_RT.dat``: Radiative transfer diagnostics.
  * ``Diagnostic_RC.dat``: Radiative-convective diagnostics.
  * ``Diagnostic_condens.dat``: Cloud condensation diagnostics.

* Copy ``config.h`` to output directory as ``config_[IN_FILE_NAME].txt`` for reproducibility.

Initial Climate Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Lines ~1040: First radiative-convective solve**

* Call ``ms_Climate()`` or ``GreyTemp()`` depending on ``RadConv_Solver`` flag.
* ``GreyTemp()`` ``GreyTemp.c``: Simple grey atmosphere radiative equilibrium (not recommended for detailed studies).
* ``ms_Climate()`` ``climate.c``: Full non-grey radiative-convective solver with clouds (see :ref:`climate_solver` section).
* Returns converged temperature profile in ``tempeq[]``.
* Copy to ``Tnew[]`` for use in subsequent iterations.

Climate-Chemistry Coupling Loop (NMAX iterations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~1050: Iterative convergence loop**

After the initial radiative-convective corvengence is achieved, for ``i = 0`` to ``NMAX-1`` iterate climate-chemistry:

a. **Check Climate-Chemistry Convergence:**
   
   * Calculate mean absolute temperature change: ``TVARTOTAL``.
   * If ``TVARTOTAL < TVARTOTAL_TOL`` (typically 1 K), convergence achieved → exit loop.

b. **Reset Arrays:**
   
   * Zero out ``clouds[][]`` array to prevent accumulation between iterations.
   * Reset cloud optical properties to transparent state ``c = 0``, ``a = 0``, ``g = 0``.

c. **Update Temperature Grid:**
   
   * ``T[j] = Tnew[j]`` for all layers.
   * Recalculate double-resolution grid ``Tdoub[]``.
   * Compute layer-center temperatures ``tl[]``.

d. **Recalculate Altitude Grid:**
   
   * Update scale heights with new temperatures.
   * Recalculate altitudes from hydrostatic equilibrium.
   * Update ``z[]``, ``zl[]``, ``zdoub[]`` arrays.

e. **Update Number Densities:**
   
   * Recalculate ``MMZ[]`` and ``MM[]`` from ideal gas law: ``MM[j] = P[j] / (k_B * T[j])``.
   * Update double-resolution grid ``MMdoub[]``.

f. **Re-compute Chemistry**:
   
   * Call ``chemquil()`` again with updated temperature profile.
   * Obtain new equilibrium mixing ratios.
   * Convert to number densities: ``xx[j][i] = mixing_ratio[j][i] * MM[j]``.

g. **Update Mean Molecular Mass:**
   
   * Recalculate ``meanmolecular[]`` with new composition.
   * Adjust Helium abundance.

h. **Save Intermediate Composition:**
   
   * Save current composition to ``ConcentrationSTD_NMAX_[i+1].dat`` for diagnostics.

i. **Run Radiative-Convective Solver:**
   
   * Call ``ms_Climate()`` with iteration counter ``nmax_iteration = i+1``.
   * Cloud physics computed inside (see :ref:`climate_solver` section).
   * Returns updated temperature profile in ``tempeq[]``.

j. **Update Temperature:**
   
   * ``Tnew[j] = tempeq[j]`` for all layers.

Finalization and Cleanup
~~~~~~~~~~~~~~~~~~~~~~~~

**Lines ~1200: Output final results and free memory**

* Write final atmospheric composition to ``ConcentrationSTD_T.dat``.
* Write final T-P profile to ``NewTemperature.dat``.
* Call cleanup routines:
  
  * ``cleanup_opacity_cache()``: Free gas opacity memory.
  * ``cleanup_cia_cache()``: Free CIA opacity memory.
  * ``free_dmatrix()``: Free all allocated 2D arrays.


.. _climate_solver:

Climate Solver ``climate.c``
----------------------------

The ``ms_Climate()`` function is the core of EPACRIS, which calls on the radiative-convective equilibrium calculations. This section describes its internal structure and operation.

Function Overview
~~~~~~~~~~~~~~~~~

**Inputs:**

* ``tempeq[]``: Initial temperature profile (overwritten with solution)
* ``P[]``: Pressure grid
* ``T[]``: Reference temperature profile
* ``Tint``: Internal heat flux temperature
* Output file paths for diagnostics
* ``nmax_iteration``: Current NMAX iteration number

Initialization
~~~~~~~~~~~~~~

**Line ~52: Variable initialization**

* Copy temperature profile to working array ``tempb[]``
* Recalculate Helium abundance for each layer: ``heliumnumber[j] = MM[j] - Σ xx[j][i]``
* Initialize condensibles mode based on ``CONDENSATION_MODE`` setting
* Detect condensible species if ``CONDENSATION_MODE = 1`` or ``2``
* Allocate memory for saturation ratios array
* Initialize convective layer arrays and convergence variables
* Set up RT stop file path and live plotting directory if enabled
* Calculate altitude grid ``znew[]`` from hydrostatic equilibrium using current temperature profile
* Compute scale height and integrate altitude from bottom to top

Radiative-Convective Iteration Loop
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Lines ~50: Main RC iteration loop ``ms_Climate``**

The function performs ``NMAX_RC`` radiative-convective iterations. Each iteration consists of:

1. **Radiative Transfer Phase**
2. **Convective Adjustment Phase** 
3. **Convergence Check**

**Radiative Transfer Phase:**

For each RC iteration ``i``:

* Set RT step limit: ``NMAX_RT`` for first iteration, ``NRT_RC`` for subsequent iterations
* Reset radiative flux array and initialize convergence status

**Line ~245: RT iteration loop**

* While not converged and ``RTstepcount < RTsteplimit``:
  
  * Reinterpolate opacities every 10 steps to update gas opacities for current T/P profile
  * Detect condensibles if ``CONDENSATION_MODE = 1`` or ``2``
  * Call ``ms_RadTrans()`` to compute upward and downward fluxes at all wavelengths and layers
  * Update temperature from flux divergence
  * Check convergence: compute ``Rfluxmax``, ``dRfluxmax``, and ``radiationO``

* Proceed to convective adjustment if ``NMAX > 0``, if ``NMAX = 0`` proceed to ``NMAX = 1``.

**Convective Adjustment Phase:**

**Line ~420: Condensation and lapse rate calculation**

* For each layer, call ``condensation_and_lapse_rate()`` to:
  
  * Calculate saturation vapor pressures for all condensible species
  * Determine condensed mass and vapor abundance
  * Compute adiabatic lapse rate using Graham et al. (2021) formulation
  * Calculate heat capacity accounting for latent heat release

**Line ~450: Cloud freezing**

* If ``FREEZE_CLOUD`` enabled and freeze condition met, store current cloud and gas abundances
* Subsequent iterations use frozen values instead of recalculating

* If ``INCLUDE_CLOUD_PHYSICS > 0``:
  
  * Call ``cloud_redistribution_none()`` or ``exponential_cloud()`` to compute particle sizes and reshape the cloud if needed
  * Call ``calculate_cloud_opacity_arrays()`` to interpolate cloud optical properties from Mie tables (see :ref:`cloud_optics`)

**Line ~515: Convective instability check and adjustment**

* Call ``ms_conv_check()`` to identify convectively unstable layers
* Mark unstable layers and group adjacent layers into convective regions
* Call ``ms_temp_adj()`` to adjust temperature in unstable regions to follow adiabatic profile
* Dry adiabat in cloud-free regions, moist adiabat in cloudy regions (see :ref:`conv_cond`)

* Repeat convective adjustment until no new unstable layers are found

**Rainout code (not fully implemented or tested):**

**Line ~620:**
* Some rainout code exists, which can selectively remove material from the atmosphere and have the alpha parameter from the adiabatic calculations readjust itself self-consistently. This code is not finished or tested.

**Line ~715: Calculation of atmospheric properties for diagnostics and plotting**

**Line ~790: Temperature smoothing**

**Line ~810: Temperature change**

**Line ~850: Final output**


Radiative Transfer Solver ``ms_RadTrans``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~30: ``ms_RadTrans()`` in ``ms_radtrans_test.c``**

* Uses either Toon et al. (1989) delta-2-stream method ``TWO_STR_SOLVER = 0`` or Heng et al. (2018) non-isothermal 2-stream method ``TWO_STR_SOLVER = 1``
* For each wavelength:
  
  * Compute optical depth from absorption and scattering coefficients
  * Solve two-stream equations for upward and downward fluxes
  * Account for gas absorption, Rayleigh scattering, collision-induced absorption, and cloud absorption/scattering

* Compute net flux and calculate heating rate from flux divergence

.. _conv_cond:

Convection and Condensation ``conv_cond_funcs.c``
---------------------------------------------------

The ``conv_cond_funcs.c`` module implements condensation thermodynamics, cloud microphysics, and adiabatic lapse rate calculations.

condensation_and_lapse_rate()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~138: Main function to calculate condensation equilibrium and adiabatic lapse rate**

**Inputs:**

* ``lay``: Layer index
* ``lapse[]``: Output array for adiabatic lapse rate (dT/dP)
* ``xxHe``: Helium number density
* ``cp``: Output pointer for heat capacity
* ``saturation_ratios[]``: Output array for saturation ratios (P_vapor / P_sat)

**Step-by-Step Algorithm for calculating condensibles and adiabatic lapse rate:**

**Step 1: Initialize variables and check cloud freezing state**

* Initialize arrays for mole fractions (``Xv[]`` condensable vapor, ``Xc[]`` condensed vapor), heat capacities, latent heat, and beta parameters
* Set ``beta[] = 0.0`` for all species initially to prevent false latent heat effects
* Check if clouds are frozen: if ``clouds_frozen == 1``, proceed to Step 2a; otherwise proceed to Step 2b

**Step 2a: Frozen cloud mode** if ``clouds_frozen == 1``

* Read frozen cloud and gas abundances directly from ``frozen_clouds[][]`` and ``frozen_xx[][]`` arrays
* Convert to mole fractions: ``Xc[i] = frozen_clouds[lay][species] / MM[lay]``, ``Xv[i] = frozen_xx[lay][species] / MM[lay]``
* Calculate dry gas mole fraction: ``Xd = 1.0 - Σ(Xv[i] + Xc[i])``
* Skip condensation calculation and proceed to Step 3

**Step 2b: Normal condensation mode** if ``clouds_frozen == 0``

* **For each condensible species ``i``:**
  
  * **Saturation vapor pressure:** Compute ``psat[i]`` using species-specific functions (``ms_psat_h2o()``, ``ms_psat_nh3()``, etc.) based on layer temperature ``tl[lay]``
  
  * **Current state:** Read existing cloud and vapor abundances: ``Xc[i] = clouds[lay][species] / MM[lay]``, ``Xv[i] = xx[lay][species] / MM[lay]``
  
  * **Total condensible:** Calculate ``Xtotal = Xv[i] + Xc[i]`` (total material available)
  
  * **Cold trapping** (if ``ENABLE_COLD_TRAP``):
    
    * Check all layers below for condensation of this species
    * If condensation found below, find minimum gas-phase abundance in condensing layers
    * Limit ``Xtotal`` to this minimum value to simulate efficient removal by settling
    * Apply "ghost cold trap fix": if current layer has no condensation but VMR is lower than layer below, set ``Xtotal = gas_below`` to fix artificial reductions
  
  * **Equilibrium partitioning:**
    
    * Calculate saturation mole fraction: ``Xv_sat = psat[i] / pl[lay]``
    * If ``Xtotal > Xv_sat`` (supersaturated):
      
      * Set vapor to saturation: ``Xv[i] = Xv_sat``
      * Condense excess: ``Xc[i] = Xtotal - Xv_sat``
    
    * If ``Xtotal ≤ Xv_sat`` (undersaturated):
      
      * All material in vapor phase: ``Xv[i] = Xtotal``, ``Xc[i] = 0.0``
  
  * **Dry gas fraction:** Update ``Xd -= Xv[i] + Xc[i]``

**Step 3: Calculate heat capacities**

* **Dry species heat capacity:**
  
  * Initialize ``cpxx_dry = 0.0`` and ``MM_dry = 0.0``
  * Add Helium contribution: ``cpxx_dry += xxHe * HeHeat(tl[lay])``, ``MM_dry += xxHe``
  * For each non-condensible species with heat capacity functions (H2O, NH3, CO, CH4, etc. if not in condensibles list), add their contributions
  * Calculate dry heat capacity: ``cp_d = cpxx_dry / MM_dry``

* **Condensible species heat capacities:**
  
  * For each condensible species, look up vapor and condensed heat capacities using temperature-dependent functions:
    
    * ``cp_v[i]``: Vapor-phase heat capacity (e.g., ``H2OHeat()``, ``NH3Heat()``)
    * ``cp_c[i]``: Condensed-phase heat capacity (e.g., ``H2O_liquid_heat_capacity()``, ``NH3_liquid_heat_capacity()``)
  
  * Get cloud retention factor: ``alpha[i] = get_global_alpha_value(lay, i)`` (layer-dependent, accounts for rainout/sedimentation)

**Step 4: Calculate latent heat and beta parameter**

* **For each condensible species:**
  
  * Calculate latent heat: ``latent[i] = ms_latent(species_id, tl[lay])`` (temperature-dependent)
  
  * **Beta parameter calculation** (critical for lapse rate):
    
    * Compute partial pressure: ``partial_pressure = Xv[i] * pl[lay]``
    * Get critical temperature ``T_crit`` for species (e.g., H2O: 647.1 K, NH3: 405.5 K)
    * **If ``partial_pressure ≥ psat[i]`` AND ``tl[lay] < T_crit``:**
      
      * Condensation is occurring: ``beta[i] = latent[i] / (R_GAS * tl[lay])``
      * This accounts for latent heat release in lapse rate
    
    * **Otherwise** (undersaturated or above critical temperature):
      
      * No condensation: ``beta[i] = 0.0``
      * Species treated as dry (no latent heat effect)

**Step 5: Calculate adiabatic lapse rate** (Graham et al. 2021, Equation 1)

* **Lapse rate numerator:** ``lapse_num = Xd + Σ Xv[i]`` (dry gas + all vapor phases)

* **Lapse rate denominator:** 
  
  * Calculate ``sum_beta_xv = Σ(beta[i] * Xv[i])`` (latent heat contribution)
  * Calculate ``big_sum_denom_num_left_term = cp_d * Xd`` (dry gas heat capacity contribution)
  * Calculate ``big_sum_denom_num_right_term = Σ(Xv[i] * (cp_v[i] - R_GAS*beta[i] + R_GAS*beta[i]²) + alpha[i] * Xc[i] * cp_c[i])`` (vapor and condensed heat capacity contributions)
  * Lapse rate denominator: ``lapse_denom = Xd * (big_sum_denom_num_left_term + big_sum_denom_num_right_term) / (R_GAS * (Xd + sum_beta_xv)) + sum_beta_xv``

* **Final lapse rate:** ``lapse[lay] = lapse_num / lapse_denom``

**Step 6: Calculate heat capacity**

* **Heat capacity numerator:** ``cp_num = cp_d * Xd + Σ(Xv[i] * cp_v[i] + alpha[i] * Xc[i] * cp_c[i])``
* **Heat capacity denominator:** ``cp_denom = Xd + Σ Xv[i]`` (only dry gas and vapor, not condensed)
* **Return heat capacity:** ``*cp = cp_num / cp_denom``

**Step 7: Update global arrays** (only if not frozen)

* **If ``clouds_frozen == 0``:**
  
  * Update vapor abundances: ``xx[lay][species] = Xv[i] * MM[lay]``
  * Update cloud abundances: ``clouds[lay][species] = Xc[i] * MM[lay]``
  
* **Calculate saturation ratios:** ``saturation_ratios[i] = (Xv[i] * pl[lay]) / psat[i]`` for diagnostic output

**Cloud Freezing:**

* When ``FREEZE_CLOUD`` is enabled and freeze condition is met, the function reads from frozen arrays instead of calculating condensation
* This preserves cloud state across iterations while still allowing correct lapse rate calculations
* Global arrays are not modified when frozen, ensuring consistency

calculate_cloud_properties()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~847: Calculate cloud particle sizes and settling velocities using Hu et al. (2019) microphysics**

**Inputs:** Gravitational acceleration, temperature, pressure, mean molecular mass, species ID, eddy diffusion coefficient, layer index

**Outputs:** Particle radii ``r0``, ``r1``, ``r2``, volume, settling velocity, scale height, mass, number density

**Step-by-Step Algorithm:**

**Step 1: Get species-specific properties**

* Look up material properties: particle density ``rho`` (kg/m³), accommodation coefficient ``acc``, molecular mass ``molecular_mass_condensible`` (AMU)
* Properties depend on species and temperature (e.g., H2O: liquid 1000 kg/m³, ice 917 kg/m³)

**Step 2: Calculate atmospheric properties**

* **Atmospheric scale height:** ``H = k_B * T / (mean_molecular_mass * AMU * g)`` (m)
* **Dimensionless fall parameter:** ``u = Kzz / H`` (ratio of diffusion to settling)
* **Atmospheric viscosity:** Calculate using Sutherland's formula with composition-dependent constants (H2, Air, CO2, N2, etc.)
* **Mean free path:** ``λ = 2μ / (P * √(8*mean_molecular_mass/(πRT)))`` (m)
* **Excess number density:** ``deltan = DELTA_P / (k_B * T)`` (molecules/m³) - supersaturation driving condensation

**Step 3: Iterative solution for equilibrium particle size**

* **Initialize:** Set Cunningham slip correction ``Cc0 = 1.0``, ventilation factor ``fa = 1.0``, distribution width ``sig = 2.0``
* **Iterate until convergence** (up to 1000 iterations):
  
  * **Condensation term:** ``cc = -48^(1/3) * π^(2/3) * D * molecular_mass * fa * deltan / rho * exp(-ln²(σ))`` (growth rate)
  * **Settling term:** ``aa = rho * g / (μ * 162^(1/3) * π^(2/3) * H) * Cc * exp(-ln²(σ))`` (settling velocity coefficient)
  * **Diffusion term:** ``bb = -u / H`` (turbulent mixing opposes settling)
  * **Solve quadratic:** ``V = [(-bb + √(bb² - 4*aa*cc)) / (2*aa)]^(3/2)`` for equilibrium volume
  * **Handle updraft-dominated case:** If no real solution, use asymptotic solution ``V = 39.9 * [μ*u*exp(ln²σ)/(ρ*g*Cc)]^(3/2)``
  * **Update slip corrections:**
    
    * Calculate Knudsen number: ``Kn = λ / d`` where ``d = (6V/π)^(1/3) * exp(-ln²(σ))``
    * Update Cunningham correction: ``Cc1 = 1 + Kn*(1.257 + 0.4*exp(-1.1/Kn))``
    * Update ventilation factor: ``fa1 = (1 + Kn) / (1 + 2*Kn*(1+Kn)/acc)``
  
  * **Check convergence:** If ``|Cc1 - Cc0| + |fa1 - fa| < 0.001``, exit loop

**Step 4: Calculate particle size distribution moments**

* **Mode radius (r0):** ``r0 = (3V/(4π))^(1/3) * exp(-1.5*ln²(σ)) * 1e6`` (μm) - smallest particles, nucleation radius
* **Surface-area-weighted radius (r1):** ``r1 = (3V/(4π))^(1/3) * exp(-ln²(σ)) * 1e6`` (μm) - used for cloud optics
* **Volume-weighted radius (r2):** ``r2 = (3V/(4π))^(1/3) * exp(-0.5*ln²(σ)) * 1e6`` (μm) - largest particles, used for sedimentation

**Step 5: Calculate settling velocity**

* **Hu+2019 formula:** ``v_fall = (ρ*g*Cc) / ((162π²)^(1/3) * μ) * V^(2/3) * exp(-ln²(σ))`` (m/s)
* **Apply correction:** ``v_d = max(v_fall - u, 0)`` (accounts for turbulent mixing)
* **Alternative Stokes law:** ``v_settle = 2*r²*ρ*g*Cc/(9*μ)`` (verification, gives same result)

**Step 6: Calculate particle number density**

* **Mass per particle:** ``mass_per_particle = V * rho`` (kg)
* **Molecules per particle:** ``molecules_per_particle = mass_per_particle / (molecular_mass * AMU)``
* **Particle number density:** ``n_density = (clouds[layer][species] / molecules_per_particle) * 1e6`` (particles/m³)

**Step 7: Return outputs**

* Return particle radii ``r0``, ``r1``, ``r2`` (μm), volume ``VP`` (cm³), settling velocity (m/s), scale height (m), mass (kg), number density (particles/m³)

cloud_redistribution_none() and exponential_cloud()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~1951: Redistribute cloud particles vertically**

**cloud_redistribution_none()** ``INCLUDE_CLOUD_PHYSICS = 1``:

* **No vertical redistribution** - clouds remain exactly where they condense
* For each layer and condensible species:
  
  * Call ``calculate_cloud_properties()`` to compute particle sizes and settling velocities
  * Store particle radii in ``particle_r0[][]``, ``particle_r1[][]``, ``particle_r2[][]`` arrays
  * Store particle number density in ``particle_number_density[][]`` array
  * Cloud abundances remain unchanged from condensation calculation

**exponential_cloud()** ``INCLUDE_CLOUD_PHYSICS = 2``:

* **Hybrid A&M (2001) + Hu+2019 transport physics:**
  
  * **Step 1:** Store original condensed distribution from ``condensation_and_lapse_rate()``
  * **Step 2:** For each condensible species independently:
    
    * Find cloud bottom layer (highest pressure with significant condensation)
    * Calculate particle properties using ``calculate_cloud_properties()`` at each layer
    * Apply A&M transport physics: redistribute condensate with exponential decay based on settling velocity and eddy diffusion
    * Calculate cloud scale height: ``H_cloud = Kzz / v_settle``
    * Redistribute: ``n_cloud(z) = n_cloud(z_bottom) * exp(-(z - z_bottom) / H_cloud)``
  
  * **Step 3:** Conserve total mass - return excess condensate to vapor phase in upper layers
  * **Step 4:** Update ``clouds[][]`` and ``xx[][]`` arrays with new distribution

ms_conv_check() and ms_temp_adj()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~588: Identify and adjust convectively unstable layers**

**ms_conv_check()**:

* **For each layer:**
  
  * Calculate actual temperature gradient: ``dT/dP`` from current temperature profile
  * Compare to adiabatic lapse rate ``lapse[j]`` calculated by ``condensation_and_lapse_rate()``
  * **If ``dT/dP > lapse[j]``:** Mark layer as convectively unstable (``isconv[j] = 1``)
  * **Otherwise:** Mark as radiative (``isconv[j] = 0``)

* **Group adjacent unstable layers:** Identify contiguous convective regions for efficient adjustment

**ms_temp_adj()**:

* **For each convective region:**
  
  * **Dry adiabat** (cloud-free regions): ``T_new = T_ref * (P_new / P_ref)^(R / (cp * μ_mean))``
  * **Moist adiabat** (cloudy regions): Adjust temperature to follow adiabatic profile accounting for latent heat release
  
  * The adiabatic profile is calculated from ``lapse[j]`` which already accounts for:
    
    * Dry gas heat capacity
    * Vapor-phase heat capacity
    * Condensed-phase heat capacity (with retention factor alpha)
    * Latent heat release (beta parameter)

* **Iterate:** Repeat convective adjustment until no new unstable layers are found (``deltaconv == 0``)

.. _cloud_optics:

Cloud Optics ``cloud_optics.c``
--------------------------------

The ``cloud_optics.c`` module handles reading Mie scattering lookup tables and computing cloud optical properties for radiative transfer.

read_cloud_optical_tables_mie()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~472: Load Mie scattering lookup tables**

* Determine format: EPACRIS format ``USE_EPACRIS_FORMAT = 1`` or LX-Mie format ``USE_EPACRIS_FORMAT = 0``
* **EPACRIS format:** Read three files per species: ``Albedo.dat``, ``Cross.dat``, ``Geo.dat``
* **LX-Mie format:** Scan directory for ``r*.dat`` files, extract particle radius from filename, read optical properties
* Store tables in ``cloud_mie_optics[]`` array for each cloud species in ``CLOUD_SPECIES_LIST``

calculate_cloud_opacity_arrays()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~788: Compute cloud optical properties**

* For each cloud species and layer:
  
  * If cloud density very small, set opacity to zero and continue
  * Get particle radius ``particle_r2`` and particle number density

* For each wavelength:
  
  * Interpolate in particle size dimension using log-space interpolation
  * Interpolate in wavelength dimension linearly
  * Handle out-of-range values by clamping to table limits

* Compute opacity: ``c[j][i] = n * σ_ext`` where ``n`` is particle number density and ``σ_ext`` is extinction cross-section
* Store in species-specific arrays: ``cH2O[][]``, ``aH2O[][]``, ``gH2O[][]``, etc.

**Units:**

* Extinction coefficient ``c``: cm⁻¹
* Albedo ``a``: dimensionless (0 = pure absorption, 1 = pure scattering)
* Asymmetry parameter ``g``: dimensionless (-1 = backward, 0 = isotropic, 1 = forward scattering)

Global Variables Reference
---------------------------

This section provides a quick reference for important global variables declared in ``epacris_main.c`` and used throughout the codebase.

Atmospheric Structure Arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``Tdoub[]``, ``Pdoub[]``, ``MMdoub[]``, ``zdoub[]``: Double-resolution grid arrays for non-isothermal layers.
* ``zl[]``, ``pl[]``, ``tl[]``: Altitude (km), pressure (Pa), and temperature (K) at layer centers.
* ``MM[]``: Number density at layer centers (molecules/cm³).
* ``meanmolecular[]``: Mean molecular mass for each layer (atomic mass units).
* ``GA``: Gravitational acceleration (m/s²).

Radiative Transfer Arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``wavelength[]``: Wavelength grid (nm).
* ``solar[]``: Stellar flux (W/m²/nm).
* ``crossr[]``: Rayleigh scattering cross-sections (cm²).
* ``opacXXX[][]``: Gas opacity arrays for each molecular species [layer][wavelength] (cm²/molecule).
* ``cH2O[][]``, ``aH2O[][]``, ``gH2O[][]``: Cloud optical properties for H₂O [layer][wavelength] (extinction in cm⁻¹, albedo dimensionless, asymmetry dimensionless).
* ``cNH3[][]``, ``aNH3[][]``, ``gNH3[][]``: Cloud optical properties for NH₃.

Collision-Induced Absorption Arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``H2H2CIA[][]``, ``H2HeCIA[][]``, ``N2N2CIA[][]``, etc.: CIA opacities [layer][wavelength] (cm²/molecule).

Chemistry Arrays
~~~~~~~~~~~~~~~~

* ``xx[][]``: Number density for each species at each layer [layer][species] (molecules/cm³).
* ``clouds[][]``: Cloud particle abundances [layer][species] (molecules/cm³).
* ``ReactionR[][]``, ``ReactionM[][]``, ``ReactionP[][]``, ``ReactionT[][]``: Chemical reaction stoichiometry arrays.

Cloud Physics Arrays
~~~~~~~~~~~~~~~~~~~~

* ``particle_r0[][]``: Mode radius (nucleation radius) [layer][condensible] (μm).
* ``particle_r1[][]``: Surface-area-weighted radius [layer][condensible] (μm) - used for cloud optics.
* ``particle_r2[][]``: Volume-weighted radius [layer][condensible] (μm).
* ``particle_number_density[][]``: Particle number density [layer][condensible] (particles/m³).
* ``fall_velocity_ms[][]``: Gravitational settling velocity [layer][condensible] (m/s).
* ``cloud_retention[][]``: Cloud retention factor [layer][condensible] (dimensionless).
