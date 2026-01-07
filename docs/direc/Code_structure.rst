Code Structure
==============

This section describes the structure of EPACRIS, the main program flow, and the organization of key routines. The code is written in C and follows a single-file compilation approach where the main file ``epacris_main.c`` includes all auxiliary C files. The program iteratively solves for the atmospheric temperature-pressure profile and chemical composition by coupling radiative transfer, convective adjustment, and chemistry modules until convergence is achieved.

Program Architecture
--------------------

EPACRIS consists of three main computational levels:

1. **Main driver** ``epacris_main.c``: Initializes atmospheric grids, loads opacity and chemistry data, manages climate-chemistry coupling iterations, and coordinates the overall simulation flow.

2. **Climate solver** ``climate.c``: Computes radiative-convective equilibrium by iteratively solving radiative transfer and applying convective adjustments. This module calls cloud physics and condensation routines.

3. **Supporting modules**: Provide specialized functionality including cloud physics ``conv_cond_funcs.c``, cloud optics ``cloud_optics.c``, opacity reading ``readcross.c``, ``readcia.c``, chemistry ``chemequil.c``, and utilities ``Interpolation.c``, ``Convert.c``.

Main Program Flow (epacris_main.c)
-----------------------------------

The following sections describe each major step in the execution sequence.

Initialization and Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~100: Parameter setup and information display**

* Initiliaze global variables and import C files
* Print run configuration information including planet properties (mass, radius, orbital distance), solver settings, and cloud physics mode.
* Initialize loop counters and local variables.
* Calculate global values like planet surface gravity ``GA`` and convert irradiation angle to radians.

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
  
  * Define pressure boundaries (``PTOP``, ``PBOTTOM``, ``PMIDDLE``, etc.).
  * Calculate equilibrium temperature at top of atmosphere from stellar irradiation.
  * Use parametric temperature profile with specified values at key pressure levels.
  * Initializes an isothermal profile with automatic equilibrium temperature if ``TTOP = 0``
  * Calculate atmospheric scale height and altitude grid from hydrostatic equilibrium.

* **TPMODE = 1**: Import T-P profile from file specified by ``TPLIST``
  
  * Read altitude (km), pressure (log₁₀(Pa)), and temperature (K) data from file.
  * Skip header lines (starting with ``#``).
  * Convert log₁₀(P) to natural log(P) for internal use.
  * Interpolate temperature onto model pressure grid.
  * Calculate layer-center values (``tl[]``, ``pl[]``, ``zl[]``).

**Altitude Grid Calculation:**

For each layer from bottom to top:

* Calculate scale height
* Set altitude grid: ``z[j] = z[j-1] - H * ln(P[j]/P[j-1])``.
* Layer centers at midpoints: ``zl[j] = (z[j] + z[j-1]) / 2``.

**Double Grid Construction:**

EPACRIS uses a double-resolution grid (``Tdoub``, ``Pdoub``, ``MMdoub``, ``zdoub``) for accurate radiative transfer in non-isothermal layers:

* Even indices (``2*j``): Layer boundaries.
* Odd indices (``2*j-1``): Layer centers.

This allows the radiative transfer solver to account for temperature variations within layers.

Chemistry Setup
~~~~~~~~~~~~~~~

**Line ~400: Load chemical species and reactions**

**Species List Import:**

* Read species file specified by ``SPECIES_LIST`` (e.g., ``Library/SpeciesList/species_HNCSO.dat``).
* Parse species properties: name, type (X/F/C/A), standard number, molecular mass, initial mixing ratio, boundary conditions.
* Classify species:
  
  * **X (Variable)**: Species solved with full continuity equation (``numx`` species).
  * **F (Fast)**: Species in photochemical equilibrium (``numf`` species).
  * **C (Constant)**: Fixed mixing ratio species (``numc`` species).
  * **A (Aerosol)**: Aerosol/haze species treated as variable (subset of X).

**Reaction List Import (Not used in this version):**

* Read reaction file specified by ``REACTION_LIST``.
* Parse reaction types:
  
  * **R**: Bimolecular reactions (``numr`` reactions).
  * **M**: Termolecular reactions (``numm`` reactions).
  * **P**: Photolysis reactions (``nump`` reactions).
  * **T**: Thermal dissociation reactions (``numt`` reactions).

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

* Allocate 2D arrays (``dmatrix``) for each molecular opacity: ``opacH2O[1..zbin][0..NLAMBDA-1]``.
* Allocate cloud optical property arrays: ``cH2O``, ``aH2O``, ``gH2O`` (and NH₃).
* Each array requires ~15 MB for typical grid sizes.

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
* Read tables for H₂-H₂, H₂-He, N₂-N₂, N₂-H₂, CO₂-CO₂ (species currently hardcoded).
* Interpolate onto model grid and store in ``XXXCIA[][]`` arrays.

**Cloud Optical Property Loading:**

* Call ``read_cloud_optical_tables_mie()`` from ``cloud_optics.c`` (see below section on cloud_optics.c).
* Load Mie scattering lookup tables for condensible species (H₂O, NH₃, etc.) specified in ``CLOUD_SPECIES_LIST``.
* Tables contain extinction cross-section, single-scattering albedo, and asymmetry parameter as function of particle size and wavelength.
* Support two formats: EPACRIS format (``USE_EPACRIS_FORMAT=1``) or LX-Mie format (``USE_EPACRIS_FORMAT=0``).

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
* ``GreyTemp()`` (``GreyTemp.c``): Simple grey atmosphere radiative equilibrium (not recommended for detailed studies).
* ``ms_Climate()`` (``climate.c``): Full non-grey radiative-convective solver with clouds (see :ref:`climate_solver` section).
* Returns converged temperature profile in ``tempeq[]``.
* Copy to ``Tnew[]`` for use in subsequent iterations.

Climate-Chemistry Coupling Loop (NMAX iterations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~1050: Iterative convergence loop**

For ``i = 0`` to ``NMAX-1``:

a. **Check Convergence:**
   
   * Calculate mean absolute temperature change: ``TVARTOTAL``.
   * If ``TVARTOTAL < TVARTOTAL_TOL`` (typically 1 K), convergence achieved → exit loop.

b. **Reset Arrays:**
   
   * Zero out ``clouds[][]`` array to prevent accumulation between iterations.
   * Reset cloud optical properties to transparent state (``c = 0``, ``a = 0``, ``g = 0``).

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

* Return 0 (successful completion).

.. _climate_solver:

Climate Solver (climate.c)
--------------------------

The ``ms_Climate()`` function is the core of EPACRIS and performs radiative-convective equilibrium calculations. This section describes its internal structure and operation.

Function Overview
~~~~~~~~~~~~~~~~~

**Inputs:**

* ``tempeq[]``: Initial temperature profile (overwritten with solution).
* ``P[]``: Pressure grid (unchanged).
* ``T[]``: Reference temperature profile (used for initialization).
* ``Tint``: Internal heat flux temperature (K).
* ``outnewtemp[]``, ``outrtdiag[]``, ``outrcdiag[]``, ``outcondiag[]``: Output file paths.
* ``nmax_iteration``: Current NMAX iteration number (for diagnostics).

**Outputs:**

* Updated temperature profile in ``tempeq[]``.
* Diagnostic files with flux profiles, heating rates, cloud properties.
* ``clouds[][]`` array populated with condensed species abundances.

Radiative-Convective Iteration Loop
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Lines ~182-850: Main RC iteration loop**

The function performs ``NMAX_RC`` radiative-convective iterations. Each iteration consists of:

1. **Radiative Transfer Phase** (lines ~211-422)
2. **Convective Adjustment Phase** (lines ~427-700)
3. **Convergence Check** (lines ~810-850)

**Radiative Transfer Phase:**

For each RC iteration ``i``:

* **Step 1: Initialize RT loop** (lines ~214-247)
  
  * Set RT step limit: ``RTsteplimit = NMAX_RT`` for first iteration, ``NRT_RC`` for subsequent iterations.
  * Reset radiative flux array ``Rflux[]``.
  * Initialize convergence status structure.

* **Step 2: RT iteration loop** (lines ~249-422)
  
  * While not converged and ``RTstepcount < RTsteplimit``:
    
    * Increment ``RTstepcount``.
    * **Reinterpolate opacities** (every 10 steps): Call ``reinterpolate_all_opacities()`` and ``reinterpolate_all_cia_opacities()`` to update gas opacities for current T/P profile.
    * **Detect condensibles** (if ``CONDENSATION_MODE = 1`` or ``2``): Call ``detect_condensibles_atmosphere()`` to identify species that can condense.
    * **Call radiative transfer solver**: ``ms_RadTrans()`` computes upward and downward fluxes at all wavelengths and layers.
    * **Update temperature**: ``tempbnew[] = tempb[] + dt[] * heating_rate`` where ``dt[]`` is time step and heating rate is computed from flux divergence.
    * **Check convergence**: Compute ``Rfluxmax`` (maximum residual flux), ``dRfluxmax`` (maximum flux gradient), and ``radiationO`` (net outgoing flux). Check if all convergence criteria are met.

* **Step 3: RT convergence** (lines ~400-422)
  
  * Print convergence diagnostics.
  * Proceed to convective adjustment even if RT converged (convection may still be needed).

**Convective Adjustment Phase:**

* **Step 1: Condensation and lapse rate calculation** (lines ~440-453)
  
  * For each layer ``j``:
    
    * Call ``condensation_and_lapse_rate()`` to:
      
      * Calculate saturation vapor pressures for all condensible species.
      * Determine condensed mass (``clouds[j][species]``) and vapor abundance (``xx[j][species]``).
      * Compute adiabatic lapse rate using Graham et al. (2021) formulation.
      * Calculate heat capacity ``cp[j]`` accounting for latent heat release.

* **Step 2: Cloud freezing** (lines ~455-466)
  
  * If ``FREEZE_CLOUD`` is enabled and freeze condition is met:
    
    * Call ``freeze_cloud_state()`` to store current cloud and gas abundances for condensible species.
    * Subsequent iterations will use frozen values instead of recalculating.

* **Step 3: Cloud physics** (lines ~498-520)
  
  * If ``INCLUDE_CLOUD_PHYSICS > 0``:
    
    * Call ``cloud_redistribution_none()`` or ``exponential_cloud()`` to compute particle sizes (``particle_r0``, ``particle_r1``, ``particle_r2``) and settling velocities.
    * Call ``calculate_cloud_opacity_arrays()`` to interpolate cloud optical properties from Mie tables.

* **Step 4: Convective instability check** (lines ~522-600)
  
  * Call ``ms_conv_check()`` to identify convectively unstable layers where ``dT/dP > (dT/dP)_adiabatic``.
  * Mark unstable layers in ``isconv[]`` array.
  * Group adjacent unstable layers into convective regions.

* **Step 5: Convective adjustment** (lines ~600-700)
  
  * Call ``ms_temp_adj()`` to adjust temperature in unstable regions to follow adiabatic profile:
    
    * Dry adiabat in cloud-free regions: ``dT/dP = g / (cp * ρ)``.
    * Moist adiabat in cloudy regions: accounts for latent heat release (see :ref:`conv_cond` section).

* **Step 6: Iterate until stable** (lines ~432-433)
  
  * Repeat convective adjustment until no new unstable layers are found (``deltaconv == 0``).

**Convergence Check:**

* **Step 1: Temperature change** (lines ~810-850)
  
  * Compute maximum temperature change between RC iterations.
  * If change is below tolerance, convergence achieved.

* **Step 2: Final output** (lines ~850-881)
  
  * Write diagnostic files.
  * Copy final temperature to ``tempeq[]``.

Radiative Transfer Solver (ms_RadTrans)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``ms_RadTrans()`` function (in ``ms_radtrans_test.c``) solves the wavelength-dependent radiative transfer equation:

**Algorithm:**

* Uses either Toon et al. (1989) delta-2-stream method (``TWO_STR_SOLVER = 0``) or Heng et al. (2018) non-isothermal 2-stream method (``TWO_STR_SOLVER = 1``).
* For each wavelength:
  
  * Compute optical depth: ``τ = (wa + ws) / (MM * μ_mean * g) * ΔP`` where ``wa`` is absorption coefficient, ``ws`` is scattering coefficient, ``MM`` is number density, ``μ_mean`` is mean molecular mass, ``g`` is gravity, and ``ΔP`` is pressure difference.
  * Solve two-stream equations for upward and downward fluxes.
  * Account for gas absorption (molecular opacities), Rayleigh scattering, collision-induced absorption, and cloud absorption/scattering.

* Compute net flux: ``F_net = F_up - F_down``.
* Calculate heating rate: ``dT/dt = -1/(ρ * cp) * dF_net/dz``.

.. _conv_cond:

Convection and Condensation (conv_cond_funcs.c)
-------------------------------------------------

The ``conv_cond_funcs.c`` module implements condensation thermodynamics, cloud microphysics, and adiabatic lapse rate calculations. This section describes its key functions.

condensation_and_lapse_rate()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose:** Calculate condensation equilibrium, cloud abundances, and adiabatic lapse rate for a single atmospheric layer.

**Inputs:**

* ``lay``: Layer index.
* ``lapse[]``: Output array for adiabatic lapse rate (dT/dP).
* ``xxHe``: Helium number density.
* ``cp``: Output pointer for heat capacity.
* ``saturation_ratios[]``: Output array for saturation ratios (P_vapor / P_sat).

**Algorithm:**

1. **Saturation vapor pressure calculation:**
   
   * For each condensible species ``i``, compute saturation vapor pressure ``P_sat[i]`` using species-specific functions (``ms_psat_h2o()``, ``ms_psat_nh3()``, etc.).
   * These functions use Antoine equations or polynomial fits.

2. **Condensation equilibrium:**
   
   * Compute partial pressure: ``P_vapor[i] = (xx[lay][i] / MM[lay]) * pl[lay]``.
   * If ``P_vapor[i] > P_sat[i]``, condense excess: ``clouds[lay][i] = (P_vapor[i] - P_sat[i]) / P_sat[i] * xx[lay][i]``.
   * Update vapor abundance: ``xx[lay][i] = xx[lay][i] - clouds[lay][i]``.
   * Compute mole fractions: ``Xv[i] = xx[lay][i] / MM[lay]`` (vapor), ``Xc[i] = clouds[lay][i] / MM[lay]`` (condensed).

3. **Latent heat calculation:**
   
   * For each condensible species, compute latent heat: ``L[i] = ms_latent(species_id, T)``.
   * Compute beta parameter: ``β[i] = L[i] / (R * T)`` where ``R`` is gas constant.

4. **Adiabatic lapse rate (Graham et al. 2021):**
   
   * Compute dry gas mole fraction: ``Xd = 1 - Σ Xv[i] - Σ Xc[i]``.
   * Compute heat capacities: ``cp_v[i]`` (vapor), ``cp_c[i]`` (condensed), ``cp_d`` (dry gas).
   * Compute lapse rate numerator: ``num = g * (Xd * cp_d + Σ Xv[i] * cp_v[i] + Σ Xc[i] * cp_c[i])``.
   * Compute lapse rate denominator: ``denom = Xd * cp_d + Σ Xv[i] * cp_v[i] + Σ Xc[i] * cp_c[i] + Σ β[i] * Xv[i] * cp_v[i]``.
   * Lapse rate: ``dT/dP = num / denom``.

5. **Heat capacity:**
   
   * ``cp = denom / (Xd + Σ Xv[i] + Σ Xc[i])``.

**Cloud Freezing:**

* If ``clouds_frozen == 1`` (clouds are frozen):
  
  * Read ``Xc[i]`` and ``Xv[i]`` from frozen arrays instead of calculating.
  * Ensure ``Xv[i] = Xv_sat[i]`` if ``Xc[i] > 0`` for correct latent heat calculation.
  * Do not update global ``xx[][]`` or ``clouds[][]`` arrays.

* If ``clouds_frozen == 0``:
  
  * Perform normal condensation calculation.
  * Update global arrays.

calculate_cloud_properties()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose:** Calculate cloud particle sizes and settling velocities using microphysics from Hu et al. (2019).

**Inputs:**

* ``g``: Gravitational acceleration.
* ``T``, ``P``: Temperature and pressure.
* ``mean_molecular_mass``: Mean molecular mass of atmosphere.
* ``condensible_species_id``: Species ID (e.g., 7 for H₂O).
* ``Kzz``: Eddy diffusion coefficient.
* ``layer``: Layer index.

**Outputs:**

* ``r0``: Mode radius (nucleation/monomer radius) in μm.
* ``r1``: Surface-area-weighted radius in μm (used for cloud optics).
* ``r2``: Volume-weighted radius in μm.
* ``VP``: Particle volume in cm³.
* ``effective_settling_velocity``: Gravitational settling velocity in m/s.
* ``scale_height``: Cloud scale height.
* ``mass_per_particle``: Particle mass in kg.
* ``n_density``: Particle number density in particles/m³.

**Algorithm:**

1. **Get particle properties:**
   
   * Look up material density, accommodation coefficient, and molecular mass for species.

2. **Particle growth equation:**
   
   * Solve balance between growth (condensation) and evaporation: ``dr/dt = α * v_th * (P_vapor - P_sat) / (4 * ρ_particle)`` where ``α`` is accommodation coefficient and ``v_th`` is thermal velocity.
   * Compute equilibrium radius ``r0`` from growth/sedimentation balance.

3. **Size distribution moments:**
   
   * Assuming log-normal distribution, compute ``r1`` and ``r2`` from ``r0`` and distribution width.

4. **Settling velocity:**
   
   * Compute Stokes drag: ``v_settle = (2/9) * g * ρ_particle * r² / (η * f_Knudsen)`` where ``η`` is dynamic viscosity and ``f_Knudsen`` is Knudsen correction factor.

5. **Cloud retention:**
   
   * Compute retention factor from balance of sedimentation and mixing: ``H_cloud / H_gas = Kzz / (v_settle * H_gas)``.

cloud_redistribution_none() and exponential_cloud()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose:** Redistribute cloud particles vertically based on transport physics.

**cloud_redistribution_none()** (``INCLUDE_CLOUD_PHYSICS = 1``):

* No vertical redistribution - clouds remain where they condense.
* Simply computes particle sizes using ``calculate_cloud_properties()``.

**exponential_cloud()** (``INCLUDE_CLOUD_PHYSICS = 2``):

* Redistributes clouds with exponential decay: ``n_cloud(z) = n_cloud(z_condense) * exp(-(z - z_condense) / H_cloud)``.
* Conserves total condensed mass.
* Returns excess condensate to vapor phase in upper layers.

ms_conv_check() and ms_temp_adj()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose:** Identify convectively unstable layers and adjust temperature to adiabatic profile.

**ms_conv_check()**:

* For each layer, compare ``dT/dP`` to adiabatic lapse rate ``lapse[j]``.
* If ``dT/dP > lapse[j]``, mark layer as convective (``isconv[j] = 1``).
* Group adjacent convective layers into regions.

**ms_temp_adj()**:

* For each convective region, adjust temperature to follow adiabatic profile:
  
  * ``T_new = T_ref * (P_new / P_ref)^(R / (cp * μ_mean))`` for dry adiabat.
  * For moist adiabat, account for latent heat release (computed in ``condensation_and_lapse_rate()``).

Cloud Optics (cloud_optics.c)
------------------------------

The ``cloud_optics.c`` module handles reading Mie scattering lookup tables and computing cloud optical properties for radiative transfer.

read_cloud_optical_tables_mie()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose:** Load Mie scattering lookup tables for all cloud species specified in ``CLOUD_SPECIES_LIST``.

**Algorithm:**

1. **Determine format:**
   
   * If ``USE_EPACRIS_FORMAT = 1``: Read EPACRIS format files (``Albedo.dat``, ``Cross.dat``, ``Geo.dat``).
   * If ``USE_EPACRIS_FORMAT = 0``: Read LX-Mie format files (``r0.010000.dat``, ``r0.012589.dat``, etc.).

2. **EPACRIS format:**
   
   * For each species, read three files:
     
     * ``Albedo.dat``: Single-scattering albedo [particle_size][wavelength].
     * ``Cross.dat``: Extinction cross-section (cm²) [particle_size][wavelength].
     * ``Geo.dat``: Asymmetry parameter [particle_size][wavelength].
   
   * Particle sizes and wavelengths are specified in file headers.

3. **LX-Mie format:**
   
   * Scan directory for files matching pattern ``r*.dat``.
   * Extract particle radius from filename.
   * Read each file: wavelength (μm), size_param, extinction (cm²), scattering (cm²), absorption (cm²), albedo, asymmetry_g.
   * Sort files by particle radius.

4. **Store in global arrays:**
   
   * Store tables in ``cloud_mie_optics[]`` array (one entry per cloud species).
   * Store species IDs in ``cloud_species_ids[]`` array.

calculate_cloud_opacity_arrays()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose:** Compute cloud optical properties (extinction coefficient, albedo, asymmetry parameter) for all layers and wavelengths using Mie tables.

**Algorithm:**

1. **Loop over cloud species:**
   
   * For each loaded cloud species (e.g., H₂O, NH₃).

2. **Loop over layers:**
   
   * For each layer ``j``:
     
     * If ``clouds[j][species_id] < 1e-12``: Set opacity to zero, albedo to 1.0, asymmetry to 0.0, continue.
     * Get particle radius: ``r = particle_r2[j][species_idx]`` (volume-weighted radius).
     * Get particle number density: ``n = particle_number_density[j][species_idx]`` (particles/m³).

3. **Loop over wavelengths:**
   
   * For each wavelength ``i`` in RT grid:
     
     * **Interpolate in particle size dimension:**
       
       * Find bounding particle sizes in Mie table: ``r_low <= r <= r_high``.
       * Interpolate extinction cross-section, albedo, and asymmetry parameter using log-space interpolation.
       * Handle out-of-range values by clamping to table limits.
     
     * **Interpolate in wavelength dimension:**
       
       * Convert RT wavelength (nm) to Mie wavelength (μm).
       * Find bounding wavelengths in Mie table.
       * Interpolate optical properties linearly in wavelength.
       * Handle out-of-range values by using nearest table value.

4. **Compute opacity:**
   
   * Extinction coefficient: ``c[j][i] = n * σ_ext`` where ``σ_ext`` is extinction cross-section (cm²) and ``n`` is converted to particles/cm³.
   * Albedo: ``a[j][i] = albedo_interpolated``.
   * Asymmetry: ``g[j][i] = asymmetry_interpolated``.

5. **Write to output arrays:**
   
   * Store in species-specific arrays: ``cH2O[][]``, ``aH2O[][]``, ``gH2O[][]``, etc.

**Units:**

* Extinction coefficient ``c``: cm⁻¹ (opacity per unit path length).
* Albedo ``a``: dimensionless (0 = pure absorption, 1 = pure scattering).
* Asymmetry parameter ``g``: dimensionless (-1 = backward scattering, 0 = isotropic, 1 = forward scattering).

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
