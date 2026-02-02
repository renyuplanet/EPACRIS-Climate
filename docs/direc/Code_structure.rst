Code Structure
==============

This section describes the structure of EPACRIS, the main program flow, and the organization of key routines. The code is written in C and follows a single-file compilation approach where the main file ``epacris_main.c`` includes all auxiliary C files. The program iteratively solves for the atmospheric temperature-pressure profile and chemical composition by coupling radiative transfer, convective adjustment, and chemistry modules until convergence is achieved.

For setup instructions, see the :doc:`Getting_started` section. For configuration parameter reference, see the :doc:`Configuration_file` section.

.. _program_architecture:

Program Architecture
--------------------

EPACRIS consists of three main computational levels:

1. **Main driver** ``epacris_main.c``: Initializes atmospheric grids, loads opacity and chemistry data, manages climate-chemistry coupling iterations, and coordinates the overall simulation flow.

2. **Climate solver** ``climate.c``: Computes radiative-convective equilibrium by iteratively solving radiative transfer and applying convective adjustments. This module calls cloud physics and condensation routines. See the :ref:`climate_solver` section for details.

3. **Supporting modules**: Provide specialized functionality including cloud physics ``conv_cond_funcs.c`` (see :ref:`conv_cond`), cloud optics ``cloud_optics.c`` (see :ref:`cloud_optics`), opacity reading ``readcross.c``, ``readcia.c``, chemistry ``chemequil.c``, and utilities ``Interpolation.c``, ``Convert.c``.

.. _main_program_flow:

Main Program Flow ``epacris_main.c``
--------------------------------------

The following sections describe each major step in the execution sequence. For configuration details, see the :doc:`Configuration_file` section. For setup instructions, see the :doc:`Getting_started` section.

.. _initialization_config:

Initialization and Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~0: Parameter setup and information display**

* Initialize global variables and import headers/C files
* Print run configuration information including:
  
  * Planet properties: mass, radius, orbital distance
  * Solver settings: radiative-convective solver type, iteration limits
  * Cloud physics mode: condensation mode, cloud physics inclusion
  * Chemistry settings: species list, reaction list
  
* Initialize loop counters and local variables
* Calculate global values:
  
  * Planet surface gravity: ``GA = GRAVITY * MASS_PLANET / RADIUS_PLANET²`` (m/s²)
  * Convert irradiation angle: ``THETAREF = THETAANGLE * π/180`` (radians)

.. _wavelength_grid:

Wavelength Grid Setup
~~~~~~~~~~~~~~~~~~~~~~

**Line ~200: Construct wavelength grid**

* Generate log-spaced wavelength array used for all radiative transfer calculations from ``LAMBDALOW`` to ``LAMBDAHIGH`` with ``NLAMBDA`` points.
* Wavelength stored in nanometers for consistency with opacity databases.

.. _rayleigh_scattering:

Rayleigh Scattering
~~~~~~~~~~~~~~~~~~~

**Line ~205: Calculate Rayleigh cross-sections**

* For each wavelength in the grid:
  
  * Compute refractive index for atmospheric gas (H2, N2, CO2, etc.) using functions from ``RefIdx.c``
  * Calculate wavelength-dependent Rayleigh scattering cross-section using refractive index
  * Store in ``crossr[]`` array [wavelength] (cm²) for use in radiative transfer calculations

.. _stellar_spectrum:

Stellar Spectrum
~~~~~~~~~~~~~~~~

**Line ~220: Load and process stellar spectrum**

* Read stellar spectrum file specified by ``STAR_SPEC`` config parameter
* Parse wavelength and flux data from file
* Interpolate stellar flux onto model wavelength grid using linear interpolation
* Scale flux from 1 AU to planet's orbital distance: ``F_planet = F_1AU / ORBIT²``
* Apply ``FaintSun`` factor to simulate planetary albedo or stellar variability: ``F_final = F_planet * FaintSun``
* Extrapolate to longer wavelengths using Rayleigh-Jeans tail if needed: ``F(λ) ∝ λ⁻⁴``
* Store final stellar flux in ``solar[]`` array [wavelength] (W/m²/nm)

.. _tp_profile_init:

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

For each layer from bottom to top (starting from surface at ``z[0] = 0``):

* Calculate scale height: ``H = k_B * T / (mean_molecular_mass * AMU * g)`` (km)
* Integrate altitude from hydrostatic equilibrium: ``z[j] = z[j-1] - H * ln(P[j]/P[j-1])``
* Calculate layer center altitudes: ``zl[j] = z[j-1] - H * ln(pl[j]/P[j-1])`` where ``pl[j]`` is mid-layer pressure

**Double Grid Construction:**

EPACRIS uses a double-resolution grid ``Tdoub``, ``Pdoub``, ``MMdoub``, ``zdoub`` for accurate radiative transfer in non-isothermal layers:

* Even indices ``2*j``: Layer boundaries.
* Odd indices ``2*j-1``: Layer centers.

This allows the radiative transfer solver to account for temperature variations within layers.

.. _chemistry_setup:

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

.. _initial_composition:

Initial Atmospheric Composition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~790: Initialize atmospheric composition**

Three initialization modes based on ``IMODE`` config setting:

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

* Sum total number density: ``n_total = Σ xx[j][i]`` for all species
* Sum total mass: ``m_total = Σ(xx[j][i] * m_i)`` where ``m_i`` is molecular mass
* Calculate Helium abundance to fill remainder: ``heliumnumber[j] = MM[j] - n_total`` (if negative, set to 0)
* Compute mean molecular mass: ``meanmolecular[j] = m_total / n_total`` (atomic mass units)
* Calculate ``MMZ[j] = meanmolecular[j] * AMU`` (kg/mol) for use in ideal gas law

.. _opacity_loading:

Opacity Loading
~~~~~~~~~~~~~~~

**Line ~990: Initialize opacity arrays and load data**

**Memory Allocation:**

* Allocate 2D arrays for each molecular opacity: ``opacXXX[zbin+1][NLAMBDA]`` where XXX is species name
* Allocate cloud optical property arrays for each defined cloud species in ``CLOUD_SPECIES_LIST``:
  
  * ``cH2O[][]``, ``aH2O[][]``, ``gH2O[][]``: H₂O cloud properties [layer][wavelength]
  * ``cNH3[][]``, ``aNH3[][]``, ``gNH3[][]``: NH₃ cloud properties [layer][wavelength]
  * Additional species arrays allocated as needed

**Gas Opacity File Reading:**

* Call ``read_all_opacities()`` from ``readcross.c`` to load gas opacities
* For each species in ``OPACITY_SPECIES_LIST``:
  
  * Build file path: ``[CROSSHEADING]/opac[species].dat`` (e.g., "../Opacity/MayTest/opacH2O.dat")
  * Read 2D opacity table: wavelength × temperature-pressure grid points
  * Cache opacity data in memory for efficient reinterpolation during RT iterations
  * Interpolate onto model's wavelength, temperature, and pressure grid
  * Store in corresponding ``opacXXX[][]`` array [layer][wavelength] (cm²/molecule)

**CIA Opacity Loading:**

* Call ``readcia()`` from ``readcia.c`` to load collision-induced absorption
* Read tables for H₂-H₂, H₂-He, H₂-H, N₂-H₂, N₂-N₂, CO₂-CO₂ (species currently hardcoded)
* Interpolate onto model grid and store in ``XXXCIA[][]`` arrays [layer][wavelength] (cm²/molecule)
* Calculate mean CIA values: ``MeanXXXCIA[]`` and solar-weighted means ``SMeanXXXCIA[]`` [layer]

**Cloud Optical Property Loading:**

* Call ``read_cloud_optical_tables_mie()`` from ``cloud_optics.c`` (see :ref:`cloud_optics` section for details)
* Load Mie scattering lookup tables for condensible species (H₂O, NH₃, etc.) specified in ``CLOUD_SPECIES_LIST``
* Tables contain extinction cross-section, single-scattering albedo, and asymmetry parameter as function of particle size and wavelength
* Support two formats: EPACRIS format (``USE_EPACRIS_FORMAT=1``) or LX-Mie format (``USE_EPACRIS_FORMAT=0``)
* Store tables in ``cloud_mie_optics[]`` structure array for later interpolation

.. _output_preparation:

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

.. _initial_climate:

Initial Climate Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Lines ~1040: First radiative-convective solve**

* Call ``ms_Climate()`` or ``GreyTemp()`` depending on ``RadConv_Solver`` flag (see :doc:`Configuration_file` for ``RadConv_Solver`` parameter).
* ``GreyTemp()`` ``GreyTemp.c``: Simple grey atmosphere radiative equilibrium (not recommended for detailed studies).
* ``ms_Climate()`` ``climate.c``: Full non-grey radiative-convective solver with clouds (see :ref:`climate_solver` section for details).
* Returns converged temperature profile in ``tempeq[]``.
* Copy to ``Tnew[]`` for use in subsequent iterations.

.. _climate_chemistry_loop:

Climate-Chemistry Coupling Loop (NMAX iterations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~1050: Iterative convergence loop**

After the initial radiative-convective convergence is achieved, perform ``NMAX`` climate-chemistry coupling iterations. For each iteration ``i = 0`` to ``NMAX-1``:

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

.. _finalization:

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

.. _rt_convergence_check:

Radiative Transfer Convergence Check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~30: ``check_rt_convergence()`` function**

The ``check_rt_convergence()`` function evaluates convergence criteria for the radiative transfer iterations within the climate solver. It checks three conditions:

* **Flux convergence**: Maximum radiative flux residual (``Rfluxmax``) must be below either an absolute tolerance (``tol_rc``) or a relative tolerance scaled by the internal heat flux (``tol_rc_r * σTint⁴``).

* **Gradient convergence**: Maximum radiative flux gradient (``dRfluxmax``) must be below the gradient tolerance (``Tol_RC_gradient``).

* **Net flux convergence**: Net outgoing flux at the top of atmosphere (``radiationO``) must satisfy the same flux tolerance criteria.

The function returns an ``RTConvergenceStatus`` structure indicating which criteria are satisfied. Overall convergence is achieved when flux convergence is satisfied along with either gradient or net flux convergence.

**Line ~50: ``ms_Climate()`` main function**

The ``ms_Climate()`` function is the core of EPACRIS, performing radiative-convective equilibrium calculations. This section describes its internal structure and operation.

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

* Define local variables and arrays for radiative-convective iterations
* Copy temperature profile to working array: ``tempb[] = T[]``
* Recalculate Helium abundance for each layer: ``heliumnumber[j] = MM[j] - Σ xx[j][i]`` for all species
* Initialize condensibles mode based on ``CONDENSATION_MODE`` setting:
  
  * ``CONDENSATION_MODE = 0``: Use predefined ``CONDENSIBLES[]`` list
  * ``CONDENSATION_MODE = 1``: Automatic detection based on saturation
  * ``CONDENSATION_MODE = 2``: Hybrid (manual list + automatic validation)
  
* Call ``initialize_condensibles_mode()`` to set up condensation detection
* Detect condensible species if ``CONDENSATION_MODE = 1`` or ``2``: call ``detect_condensibles_atmosphere()``
* Allocate memory for saturation ratios array: ``saturation_ratios[zbin+1][MAX_CONDENSIBLES]``
* Initialize convective layer arrays: ``isconv[] = 0`` (all radiative initially)
* Initialize convergence variables: ``RTstepcount = 0``, convergence status structure
* Set up RT stop file path: ``[OUT_DIR]/rt.stop`` (allows early termination)
* Create live plotting directory if ``LIVE_PLOTTING`` enabled: ``[OUT_DIR]/live_plot/``
* Calculate altitude grid ``znew[]`` from hydrostatic equilibrium using current temperature profile:
  
  * ``znew[0] = 0.0`` (surface)
  * For each layer: ``znew[j] = znew[j-1] - scaleheight * ln(P[j]/P[j-1])``
  * Scale height: ``scaleheight = k_B * T / (meanmolecular * AMU * GA) / 1000`` (km)

Radiative-Convective Iteration Loop
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Lines ~180: Main Radiative-Convective iteration loop ``ms_Climate``**

The function performs ``NMAX_RC`` radiative-convective iterations. Each iteration consists of:

1. **Radiative Transfer Phase**
2. **Convective Adjustment Phase** 
3. **Convergence Check**

**Radiative Transfer Phase:**

For each RC iteration ``i``:

* **Set RT step limit:**
  
  * First iteration (``i == 1``): ``RTsteplimit = NMAX_RT`` (typically 800 steps)
  * Subsequent iterations (``i > 1``): ``RTsteplimit = NRT_RC`` (This value will significantly impact the results if the radiative solver is not converging!)
  
* Reset radiative flux array: ``Rflux[] = 0.0``
* Initialize convergence status structure: ``RTConvergenceStatus status = {0}``

**Line ~245: RT iteration loop**

* While not converged and ``RTstepcount < RTsteplimit`` and RT stop file doesn't exist:
  
  * **Opacity reinterpolation** (every 10 steps):
    
    * Call ``reinterpolate_all_opacities()`` to update gas opacities for current T/P profile
    * Call ``reinterpolate_all_cia_opacities()`` to update CIA opacities
    * Ensures opacities match current atmospheric state
  
  * **Dynamic condensation detection** (if ``CONDENSATION_MODE = 1`` or ``2``):
    
    * Call ``detect_condensibles_atmosphere()`` if ``CONDENSATION_TIMING = 1`` or every ``NRT_RC`` steps
    * Updates ``NCONDENSIBLES`` and ``CONDENSIBLES[]`` array if new condensible species detected
  
  * **Radiative transfer calculation:**
    
    * Call ``ms_RadTrans()`` to compute upward and downward fluxes at all wavelengths and layers
    * Updates temperature profile ``tempbnew[]`` from flux divergence
    * Returns diagnostic fluxes: ``radiationI0`` (TOA incoming), ``radiationI1`` (BOA incoming), ``radiationO`` (TOA outgoing)
  
  * **Convergence check:**
    
    * Calculate maximum radiative flux residual: ``Rfluxmax = max(|Rflux[j]|)`` for radiative layers
    * Calculate maximum flux gradient: ``dRfluxmax = max(|Rflux[j-1] - Rflux[j]|)``
    * Call ``check_rt_convergence()`` to evaluate convergence criteria
    * Exit loop if converged or step limit reached

* After radiative transfer phase completes (converged or step limit reached), proceed to convective adjustment phase.

**Convective Adjustment Phase:**

**Line ~420: Convection loop**

* **Condensation and lapse rate calculation:**
  
  * For each layer, call ``condensation_and_lapse_rate()`` from ``conv_cond_funcs.c`` (see :ref:`condensation_lapse_rate` section):
    
    * Calculate saturation vapor pressures for all condensible species using temperature-dependent functions
    * Determine condensed mass and vapor abundance based on saturation equilibrium
    * Apply cold trapping if ``ENABLE_COLD_TRAP = 1`` to limit vapor above condensation regions
    * Compute adiabatic lapse rate using Graham et al. (2021) formulation accounting for:
      
      * Dry gas heat capacity
      * Vapor-phase heat capacity
      * Condensed-phase heat capacity (with retention factor alpha)
      * Latent heat release (beta parameter)
    
    * Calculate heat capacity accounting for latent heat release
    * Update global arrays ``xx[][]`` and ``clouds[][]`` (unless clouds are frozen)
    * Store saturation ratios in ``saturation_ratios[][]`` array for diagnostics

**Line ~450: Cloud freezing**

* If ``FREEZE_CLOUD`` enabled and freeze condition met, ``freeze_cloud_state()`` stores:
  
  * Cloud abundances (``clouds[j][i]``) for all species
  * Gas phase abundances (``xx[j][i]``) for all species
  
* In subsequent iterations, ``restore_frozen_clouds()`` restores frozen values only for condensible species, preventing their cloud and gas abundances from evolving further

**Line ~463: Recalculate helium after condensation changes**

* After condensation calculation updates ``xx[][]`` array, recalculate Helium abundance:
  
  * ``heliumnumber[j] = MM[j] - Σ xx[j][i]`` for all species
  * Ensures Helium abundance accounts for material moved from gas phase to clouds

**Line ~490: Cloud physics calculation**

* If ``INCLUDE_CLOUD_PHYSICS > 0``:
  
  * Call ``store_cloud_properties()`` to compute particle sizes and number density for cloud optics
    
    * Calculates particle radii (r0, r1, r2), volume, mass, and number density
    * Stores ``particle_r2[][]`` and ``particle_number_density[][]`` in global arrays
    * Note: ``exponential_cloud()`` exists but is incomplete and not currently used
  
  * Call ``calculate_cloud_opacity_arrays()`` to interpolate cloud optical properties from Mie tables (see :ref:`cloud_optics` section for details)
    
    * Interpolates Mie table data to current particle sizes and RT wavelength grid
    * Computes extinction opacity, albedo, and asymmetry parameter for each layer and wavelength
    * Stores results in ``cH2O[][]``, ``aH2O[][]``, ``gH2O[][]``, etc. arrays

**Line ~515: Convective instability check and adjustment**

* **Convective instability detection:**
  
  * Call ``ms_conv_check()`` to identify convectively unstable layers (see :ref:`convective_adjustment` section for details)
    
    * For each layer, compare actual temperature gradient ``dT/dP`` to adiabatic lapse rate ``lapse[j]``
    * If ``dT/dP > lapse[j]``: Mark layer as convectively unstable (``isconv[j] = 1``)
    * Group adjacent unstable layers into contiguous convective regions
  
  * Count convective layers: ``ncl`` (convective), ``nrl`` (radiative)

* **Temperature adjustment:**
  
  * Call ``ms_temp_adj()`` to adjust temperature in unstable regions to follow adiabatic profile (see :ref:`convective_adjustment` section)
    
    * **Dry adiabat** (cloud-free regions): ``T_new = T_ref * (P_new / P_ref)^(R/(cp*μ_mean))``
    * **Moist adiabat** (cloudy regions): Adjust temperature using ``lapse[j]`` which accounts for:
      
      * Dry gas heat capacity
      * Vapor-phase heat capacity  
      * Condensed-phase heat capacity (with retention factor alpha)
      * Latent heat release (beta parameter)
    
    * Adjust temperature to maintain adiabatic profile within each convective region
    * Calculate potential temperature for each convective region to ensure consistency

* **Iteration:**
  
  * Repeat convective adjustment loop until no new unstable layers are found (``deltaconv == 0``)
  * This ensures all unstable regions are properly adjusted

**Line ~620: Rainout code (not fully implemented or tested)**

* Some rainout code exists in ``simulate_rainout()`` function which can selectively remove material from the atmosphere
* The alpha parameter from adiabatic calculations can readjust itself self-consistently
* **Status:** This code is not finished or tested and may not work correctly

**Line ~715: Calculation of atmospheric properties for diagnostics and plotting**

* Calculate properties needed for live plotting (if ``LIVE_PLOTTING`` enabled):
  
  * Update saturation ratios for all condensible species
  * Calculate potential temperature for convective regions
  * Prepare diagnostic strings with current RT step information

**Line ~790: Temperature smoothing**

* Apply optional temperature smoothing to reduce numerical oscillations
* Update double-resolution grid ``Tdoub[]`` with smoothed temperatures

**Line ~810: Convergence check**

* Check if radiative-convective iteration has converged:
  
  * Compare current temperature profile to previous iteration
  * Check if convective boundaries have stabilized
  * Determine if additional RC iterations are needed

**Line ~825: End of Radiative-Convective iteration loop and final output**

* After ``NMAX_RC`` iterations or convergence:
  
  * Write final temperature profile to output file
  * Write radiative transfer diagnostics
  * Write convective adjustment diagnostics
  * Write condensation diagnostics (if enabled)


.. _radiative_transfer_solver:

Radiative Transfer Solver ``ms_RadTrans``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~30: ``ms_RadTrans()`` in ``ms_radtrans_test.c``**

Called from ``climate.c`` during radiative transfer iterations (see :ref:`climate_solver` section) to compute radiative fluxes and heating rates.

**Function Overview:**

* Solves two-stream radiative transfer equations for upward and downward fluxes
* Accounts for multiple opacity sources: gas absorption, Rayleigh scattering, CIA, and cloud scattering/absorption
* Computes heating rates from flux divergence to update temperature profile

**Algorithm:**

For each wavelength in the grid:

1. **Compute optical properties:**
   
   * Gas absorption: Sum over all species ``κ_gas = Σ(xx[layer][species] * opacXXX[layer][wavelength])``
   * Rayleigh scattering: ``σ_rayleigh = crossr[wavelength] * MM[layer]``
   * CIA: Sum collision-induced absorption contributions ``κ_CIA = Σ(xx[layer][species1] * xx[layer][species2] * XXXCIA[layer][wavelength])``
   * Cloud scattering/absorption: ``κ_cloud = cH2O[layer][wavelength] + cNH3[layer][wavelength] + ...``
   * Total extinction: ``κ_total = κ_gas + κ_CIA + κ_cloud + σ_rayleigh``

2. **Solve two-stream equations:**
   
   * **Method selection:**
     
     * ``TWO_STR_SOLVER = 0``: Toon et al. (1989) delta-2-stream method (isothermal layers)
     * ``TWO_STR_SOLVER = 1``: Heng et al. (2018) non-isothermal 2-stream method (accounts for temperature gradients)
   
   * Compute optical depth: ``τ = κ_total * Δz`` where ``Δz`` is layer thickness
   * Solve two-stream equations for upward ``F↑`` and downward ``F↓`` fluxes
   * Apply boundary conditions:
     
     * Top of atmosphere: Incoming stellar flux ``F↓(TOA) = solar[wavelength] * cos(THETAREF)``
     * Bottom: Surface reflection/emission based on ``PSURFAB`` and ``PSURFEM``

3. **Compute heating rates:**
   
   * Calculate net flux: ``F_net = F↑ - F↓``
   * Compute flux divergence: ``dF/dz = (F_net[j] - F_net[j-1]) / Δz``
   * Calculate heating rate: ``dT/dt = -(1/ρcp) * dF/dz`` where ``ρ`` is density and ``cp`` is heat capacity
   * Update temperature: ``T_new = T_old + dT/dt * dt`` where ``dt`` is time step

**Outputs:**

* Updated temperature profile ``tempbnew[]``
* Radiative flux array ``Rflux[]`` for convergence checking
* Diagnostic outputs: ``radiationI0`` (TOA incoming), ``radiationI1`` (BOA incoming), ``radiationO`` (TOA outgoing)

.. _conv_cond:

Convection and Condensation ``conv_cond_funcs.c``
---------------------------------------------------

The ``conv_cond_funcs.c`` module implements condensation thermodynamics, cloud microphysics, and adiabatic lapse rate calculations.

Function Overview
~~~~~~~~~~~~~~~~~~

**Main Functions:**

* ``condensation_and_lapse_rate()``: Calculates condensation equilibrium and adiabatic lapse rate for each layer
* ``calculate_cloud_properties()``: Computes cloud particle sizes and settling velocities using Hu et al. (2019) microphysics
* ``store_cloud_properties()``: Stores particle properties for cloud optics calculations
* ``ms_conv_check()``: Identifies convectively unstable layers
* ``ms_temp_adj()``: Adjusts temperature in unstable regions to follow adiabatic profile

**Supporting Functions:**

* ``get_particle_properties()``: Retrieves material properties (density, accommodation coefficient) for cloud species
* ``ms_psat_XXX()``: Species-specific saturation vapor pressure functions
* ``ms_latent()``: Latent heat of vaporization/sublimation functions
* Cloud freezing functions: ``freeze_cloud_state()``, ``restore_frozen_clouds()``
* Condensation detection functions: ``detect_condensibles_atmosphere()``, ``check_species_condensible()``

.. _condensation_lapse_rate:

condensation_and_lapse_rate()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~138: Main function to calculate condensation equilibrium and adiabatic lapse rate**

**Inputs:**

* ``lay``: Layer index
* ``lapse[]``: Output array for adiabatic lapse rate (dT/dP)
* ``xxHe``: Helium number density
* ``cp``: Output pointer for heat capacity
* ``saturation_ratios[]``: Output array for saturation ratios (P_vapor / P_sat)

**Algorithm:**

**Step 1: Initialize variables and check cloud freezing state**

* Initialize arrays for mole fractions:
  
  * ``Xv[]``: Condensible vapor mole fractions
  * ``Xc[]``: Condensed vapor mole fractions
  * Heat capacities: ``cp_v[]`` (vapor), ``cp_c[]`` (condensed)
  * Latent heat: ``latent[]``
  * Beta parameters: ``beta[]`` (latent heat contribution to lapse rate)
  
* Set ``beta[] = 0.0`` for all species initially to prevent false latent heat effects
* Initialize dry gas mole fraction: ``Xd = 1.0`` (will be reduced by condensible fractions)
* Check if clouds are frozen: if ``clouds_frozen == 1``, proceed to Step 2a; otherwise proceed to Step 2b

**Step 2a: Frozen cloud mode** if ``clouds_frozen == 1``

* Read frozen cloud and gas abundances directly from ``frozen_clouds[][]`` and ``frozen_xx[][]`` arrays
  
  * Frozen state set by ``freeze_cloud_state()`` function (see :ref:`climate_solver` section)
  
* Convert to mole fractions: ``Xc[i] = frozen_clouds[lay][species] / MM[lay]``, ``Xv[i] = frozen_xx[lay][species] / MM[lay]``
* Calculate dry gas mole fraction: ``Xd -= Xv[i] + Xc[i]``
* Skip condensation calculation and proceed to Step 3
  
  * Cloud abundances remain fixed, but lapse rate still calculated correctly

**Step 2b: Normal condensation mode** if ``clouds_frozen == 0``

* **For each condensible species ``i``:**
  
  * **Saturation vapor pressure:** Compute ``psat[i]`` using species-specific functions:
    
    * ``ms_psat_h2o()``, ``ms_psat_nh3()``, ``ms_psat_co()``, ``ms_psat_ch4()``, etc.
    * Based on layer temperature ``tl[lay]``
    * Temperature-dependent (e.g., Antoine equation or empirical fits)
  
  * **Current state:** Read existing cloud and vapor abundances:
    
    * ``Xc[i] = clouds[lay][species] / MM[lay]`` (condensed mole fraction)
    * ``Xv[i] = xx[lay][species] / MM[lay]`` (vapor mole fraction)
  
  * **Total condensible:** Calculate ``Xtotal = Xv[i] + Xc[i]`` (total condensible material available)
  
  * **Cold trapping** (if ``ENABLE_COLD_TRAP = 1``, see :doc:`Configuration_file`):
    
    * Check all layers below current layer for condensation of this species
    * If condensation found below, find minimum gas-phase abundance in condensing layers
    * Limit ``Xtotal`` to this minimum value to simulate efficient removal by settling/rainout
    * Apply "ghost cold trap fix": if current layer has no condensation but VMR is lower than layer below, set ``Xtotal = gas_below`` to fix artificial reductions
    * Store original ``Xtotal`` using ``store_original_xtotal()`` for tracking
    
    * Cold trapping prevents unrealistic vapor buildup above condensation regions
  
  * **Equilibrium partitioning:**
    
    * Calculate saturation mole fraction: ``Xv_sat = psat[i] / pl[lay]``
    * **If ``Xtotal > Xv_sat`` (supersaturated):**
      
      * Set vapor to saturation: ``Xv[i] = Xv_sat``
      * Condense excess: ``Xc[i] = Xtotal - Xv_sat``
      * Condensation is occurring (needed for beta calculation)
    
    * **If ``Xtotal ≤ Xv_sat`` (undersaturated):**
      
      * All material in vapor phase: ``Xv[i] = Xtotal``, ``Xc[i] = 0.0``
      * No condensation (beta will be zero)
  
  * **Dry gas fraction:** Update ``Xd -= Xv[i] + Xc[i]`` (remove condensible from dry gas)

**Step 3: Calculate heat capacities**

* **Dry species heat capacity:**
  
  * Initialize ``cpxx_dry = 0.0`` and ``MM_dry = 0.0``
  * Add Helium contribution: ``cpxx_dry += xxHe * HeHeat(tl[lay])``, ``MM_dry += xxHe``
  * For each non-condensible species with heat capacity functions (CO, CH₄, N₂, etc. if not in condensibles list), add their contributions:
    
    * ``cpxx_dry += xx[lay][species] * SpeciesHeat(tl[lay])``
    * ``MM_dry += xx[lay][species]``
  
  * Calculate average dry heat capacity: ``cp_d = cpxx_dry / MM_dry`` (J/(kg·K))

* **Condensible species heat capacities:**
  
  * For each condensible species, look up vapor and condensed heat capacities using temperature-dependent functions:
    
    * ``cp_v[i]``: Vapor-phase heat capacity (e.g., ``H2OHeat()``, ``NH3Heat()``) [J/(kg·K)]
    * ``cp_c[i]``: Condensed-phase heat capacity (e.g., ``H2O_liquid_heat_capacity()``, ``NH3_liquid_heat_capacity()``) [J/(kg·K)]
  
  * Read cloud retention factor for each species: ``alpha[i] = get_global_alpha_value(lay, i)``
    
    * Layer-dependent factor (0-1) accounting for rainout/sedimentation
    * ``alpha = 1.0``: All condensed material contributes to heat capacity
    * ``alpha < 1.0``: Reduced contribution due to material removal
    * Controlled by ``ALPHA_RAINOUT`` config parameter (see :doc:`Configuration_file`)

**Step 4: Calculate latent heat and beta parameter**

* **For each condensible species:**
  
  * Calculate latent heat: ``latent[i] = ms_latent(species_id, tl[lay])`` (temperature-dependent, J/kg)
  
  * **Beta parameter calculation** (critical for lapse rate):
    
    * Compute partial pressure of condensible species: ``partial_pressure = Xv[i] * pl[lay]`` (Pa)
    * Get critical temperature ``T_crit`` for species (e.g., H₂O: 647.1 K, NH₃: 405.5 K)
    * **If ``partial_pressure ≥ psat[i]`` AND ``tl[lay] < T_crit``:**
      
      * Condensation is occurring (supersaturated and below critical point)
      * Set beta parameter: ``beta[i] = latent[i] / (R_GAS * tl[lay])`` (dimensionless)
      * This accounts for latent heat release in lapse rate calculation
    
    * **Otherwise** (undersaturated or above critical temperature):
      
      * No condensation: ``beta[i] = 0.0``
      * Species treated as dry (no latent heat effect on lapse rate)

**Step 5: Calculate adiabatic lapse rate** (Following Graham et al. 2021, Equation 1)

* **Lapse rate numerator:** ``lapse_num = Xd + Σ Xv[i]`` (dry gas + all vapor phases, dimensionless)

* **Lapse rate denominator:** 
  
  * Calculate latent heat contribution: ``sum_beta_xv = Σ(beta[i] * Xv[i])``
  * Calculate dry gas heat capacity contribution: ``big_sum_denom_num_left_term = cp_d * Xd``
  * Calculate vapor and condensed heat capacity contributions:
    
    * ``big_sum_denom_num_right_term = Σ(Xv[i] * (cp_v[i] - R_GAS*beta[i] + R_GAS*beta[i]²) + alpha[i] * Xc[i] * cp_c[i])``
    * Accounts for vapor-phase heat capacity with latent heat corrections
    * Accounts for condensed-phase heat capacity weighted by retention factor
  
  * Lapse rate denominator: ``lapse_denom = Xd * (big_sum_denom_num_left_term + big_sum_denom_num_right_term) / (R_GAS * (Xd + sum_beta_xv)) + sum_beta_xv``

* **Final lapse rate:** ``lapse[lay] = lapse_num / lapse_denom`` (dimensionless, dT/dP)
  
  * Represents temperature change per pressure change in adiabatic process
  * Larger values indicate steeper temperature gradients

**Step 6: Calculate heat capacity**

* **Heat capacity numerator:** ``cp_num = cp_d * Xd + Σ(Xv[i] * cp_v[i] + alpha[i] * Xc[i] * cp_c[i])``
  
  * Weighted sum of dry gas, vapor-phase, and condensed-phase heat capacities
  * Condensed-phase contribution reduced by retention factor ``alpha[i]``

* **Heat capacity denominator:** ``cp_denom = Xd + Σ Xv[i]`` (only dry gas and vapor, not condensed)
  
  * Normalizes by total gas-phase mole fraction

* **Return heat capacity:** ``*cp = cp_num / cp_denom`` (J/(kg·K))
  
  * Effective heat capacity accounting for all phases and retention factors

**Step 7: Update global abundance arrays** (only if clouds are not frozen via ``FREEZE_CLOUD``)

* **If ``clouds_frozen == 0``:**
  
  * Update vapor abundances: ``xx[lay][species] = Xv[i] * MM[lay]`` (molecules/cm³)
  * Update cloud abundances: ``clouds[lay][species] = Xc[i] * MM[lay]`` (molecules/cm³)
  
* **Calculate saturation ratios:** ``saturation_ratios[i] = (Xv[i] * pl[lay]) / psat[i]``
  
  * Ratio of actual vapor pressure to saturation vapor pressure
  * ``saturation_ratio = 1.0``: Saturated
  * ``saturation_ratio > 1.0``: Supersaturated (condensation occurring)
  * ``saturation_ratio < 1.0``: Undersaturated (all vapor)
  * Used for diagnostic output and condensation detection

**Cloud Freezing:**

* When ``FREEZE_CLOUD`` is enabled and freeze condition is met, the function reads from frozen arrays instead of calculating condensation
* This preserves cloud state across iterations while still allowing correct lapse rate calculations
* Global arrays are not modified when frozen, ensuring consistency

**Line ~590: ``simulate_rainout()`` function**

* **Purpose:** Simulates removal of condensible material from the atmosphere through rainout/sedimentation
* **Status:** Not fully tested - may or may not work correctly
* **Functionality:** 
  
  * Calculates mass loss ratio for each layer based on settling velocities
  * Removes condensed material from upper layers
  * Adjusts alpha parameter (cloud retention factor) self-consistently


calculate_cloud_properties()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_particle_properties()
~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~775: Get material properties for cloud species**

* **Purpose:** Retrieves species-specific material properties needed for cloud microphysics calculations
* **Inputs:** Species ID, temperature
* **Outputs:** Particle density ``rho`` (kg/m³), accommodation coefficient ``acc`` (dimensionless), molecular mass (AMU)
* **Properties:** Temperature-dependent (e.g., H₂O liquid density ~1000 kg/m³, ice density ~917 kg/m³)

.. _calculate_cloud_properties:

calculate_cloud_properties()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~845: Calculate cloud particle sizes and settling velocities using Hu et al. (2019) microphysics**

**Inputs:**

* ``T``: Temperature (K)
* ``P``: Pressure (Pa)
* ``mean_molecular_mass``: Mean molecular mass of atmosphere (AMU)
* ``condensible_species_id``: Species ID (e.g., 7 for H₂O, 9 for NH₃)
* ``Kzz``: Eddy diffusion coefficient (m²/s)
* ``layer``: Layer index

**Outputs:**

* ``r0``, ``r1``, ``r2``: Particle radii (μm) - mode, surface-area-weighted, volume-weighted
* ``VP``: Particle volume (cm³)
* ``effective_settling_velocity``: Gravitational settling velocity (m/s)
* ``scale_height``: Atmospheric scale height (m)
* ``mass_per_particle``: Particle mass (kg)
* ``n_density``: Particle number density (particles/m³)

**Algorithm:**

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
* **Surface-area-weighted radius (r1):** ``r1 = (3V/(4π))^(1/3) * exp(-ln²(σ)) * 1e6`` (μm) - intermediate size
* **Volume-weighted radius (r2):** ``r2 = (3V/(4π))^(1/3) * exp(-0.5*ln²(σ)) * 1e6`` (μm) - largest particles, used for cloud optics and opacity interpolation (see ``cloud_optics.c`` line ~850)

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

store_cloud_properties()
~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~1103: Calculate and store cloud particle properties for optics**

Called from ``climate.c`` (see :ref:`climate_solver` section) after condensation calculation (see :ref:`condensation_lapse_rate` section) to prepare particle properties for cloud optics.

**store_cloud_properties()** ``INCLUDE_CLOUD_PHYSICS = 1``:

* **Purpose:** Calculate and store cloud particle physical properties (sizes, number density) needed for cloud optics calculations. This function does NOT redistribute clouds - it only stores properties based on the current cloud distribution from condensation.

* For each layer and condensible species:
  
  * Call ``calculate_cloud_properties()`` to compute particle sizes (r0, r1, r2), volume, mass, and number density
  * Store particle radii in ``particle_r2[][]`` array (used by cloud optics)
  * Store particle number density in ``particle_number_density[][]`` array (used by cloud optics)
  * Cloud abundances (``clouds[][]``) remain unchanged from condensation calculation
  * Particle properties are stored in global arrays for use by ``calculate_cloud_opacity_arrays()``

**exponential_cloud()** ``INCLUDE_CLOUD_PHYSICS = 2``:

* **Status:** The exponential cloud function exists in the codebase but is currently incomplete and not fully integrated. It is located at the end of ``conv_cond_funcs.c`` and needs to be properly linked to the cloud physics and optics calculation pipeline. It should produce an exponential cloud but the optical properties will not be correct.

* **Intended functionality:** When fully implemented, this function will apply Ackerman & Marley (2001) transport physics to redistribute cloud particles vertically based on settling velocities and eddy diffusion, accounting for particle growth and sedimentation.
  
* **Note:** Worth developing is custom cloud shapes are needed.

.. _convective_adjustment:

ms_conv_check() and ms_temp_adj()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~1170: Identify and adjust convectively unstable layers**

Called from ``climate.c`` (see :ref:`climate_solver` section) during convective adjustment phase to ensure atmospheric stability.

**ms_conv_check()**:

* **For each layer:**
  
  * Calculate actual temperature gradient: ``dT/dP`` from current temperature profile
  * Compare to adiabatic lapse rate ``lapse[j]`` calculated by ``condensation_and_lapse_rate()`` (see :ref:`condensation_lapse_rate` section)
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

The ``cloud_optics.c`` module handles reading Mie scattering lookup tables and computing cloud optical properties for radiative transfer. It supports two Mie table formats: EPACRIS format (separate Albedo.dat, Cross.dat, Geo.dat files) and LX-Mie format (radius-keyed files).

Function Overview
~~~~~~~~~~~~~~~~~~

**Main Functions:**

* ``read_cloud_optical_tables_mie()``: Entry point for loading Mie tables (called from ``epacris_main.c``)
* ``calculate_cloud_opacity_arrays()``: Interpolates optical properties to RT wavelength grid and computes layer opacities

**Supporting Functions:**

* ``read_cloud_optical_tables_epacris()``: Reads EPACRIS format tables
* ``read_epacris_species_tables()``: Reads tables for a single species in EPACRIS format
* ``read_epacris_mie_file()``: Reads individual EPACRIS files (Albedo.dat, Cross.dat, Geo.dat)
* ``read_mie_file()``: Reads LX-Mie format files
* ``get_cloud_output_arrays()``: Maps species ID to output opacity arrays
* Interpolation functions: ``interp1_log()``, ``interp1_loglog()``, ``interp1_wavelength()``

read_cloud_optical_tables_mie()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~472: Load Mie scattering lookup tables**

Called once during initialization from ``epacris_main.c``. Determines which format to use based on ``USE_EPACRIS_FORMAT`` config setting:

* If ``USE_EPACRIS_FORMAT = 1``: Calls ``read_cloud_optical_tables_epacris()``
* If ``USE_EPACRIS_FORMAT = 0``: Reads LX-Mie format files

**LX-Mie Format (``USE_EPACRIS_FORMAT = 0``):**

* Reads cloud species list from ``CLOUD_SPECIES_LIST`` config macro
* For each species:
  
  * Builds directory path: ``CLOUD_MIE_DIRECTORY_LXMIE/<species_name>/``
  * Scans directory for radius-keyed files (``r0.010000.dat``, ``r0.012589.dat``, etc.)
  * Extracts particle radius from filename using ``extract_radius_from_filename()``
  * Sorts files by radius using ``compare_particle_files()``
  * Reads each file using ``read_mie_file()`` which parses:
    
    * Wavelength (μm), size parameter, extinction cross-section (cm²), scattering (cm²), absorption (cm²), albedo, asymmetry parameter g
  
  * Stores data in ``CloudOpticalTableMie`` structure arrays

**EPACRIS Format (``USE_EPACRIS_FORMAT = 1``):**

* Calls ``read_cloud_optical_tables_epacris()`` which:
  
  * Reads cloud species list from ``CLOUD_SPECIES_LIST`` config macro
  * For each species, calls ``read_epacris_species_tables()``:
    
    * Builds directory path: ``CLOUD_MIE_DIRECTORY_EPACRIS/<species_name>/``
    * Reads three files: ``Albedo.dat``, ``Cross.dat``, ``Geo.dat``
    * Each file contains particle radii in first column, optical properties for 1387 wavelengths in subsequent columns
    * Uses ``read_epacris_mie_file()`` to parse each file
    * Generates wavelength grid using ``generate_epacris_wavelength_grid()``
    * Verifies all three files have matching particle size counts
    * Stores data in ``CloudOpticalTableMie`` structure

**Output:**

* Populates global arrays ``cloud_mie_optics[]`` and ``cloud_species_ids[]``
* Sets ``n_cloud_species_loaded`` counter

calculate_cloud_opacity_arrays()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Line ~788: Compute cloud optical properties for radiative transfer**

Called from ``climate.c`` (see :ref:`climate_solver` section) after cloud particle properties are calculated via ``store_cloud_properties()`` (see :ref:`store_cloud_properties` section). Interpolates Mie table data to compute layer-by-layer cloud opacities.

**Inputs (global arrays):**

* ``clouds[][]``: Cloud number density [molecules/cm³]
* ``particle_r2[][]``: Volume-weighted particle radius [μm] 
* ``particle_number_density[][]``: Particle number density [particles/m³]
* ``wavelength[]``: RT wavelength grid [nm]

**Outputs (global arrays):**

* ``cH2O[][]``, ``cNH3[][]``: Extinction cross-section [cm⁻¹] per layer and wavelength
* ``aH2O[][]``, ``aNH3[][]``: Single-scattering albedo [dimensionless]
* ``gH2O[][]``, ``gNH3[][]``: Asymmetry parameter [dimensionless]

**Algorithm:**

For each loaded cloud species and each atmospheric layer:

1. **Check for clouds:** If ``clouds[layer][species] < 1e-12``, set opacity to zero and continue

2. **Get particle properties:**
   
   * Particle radius: ``particle_r2[layer][species_idx]`` [μm]
   * Particle number density: ``particle_number_density[layer][species_idx]`` [particles/m³] → convert to [particles/cm³]

3. **Interpolate in particle size dimension:**
   
   * For each Mie table wavelength:
     
     * Extract extinction cross-section, albedo, asymmetry across all particle sizes
     * Interpolate to current particle radius:
       
       * Extinction: log-log interpolation (``interp1_loglog()``)
       * Albedo: linear in log-radius space (``interp1_log()``)
       * Asymmetry: linear in log-radius space (``interp1_log()``)
     
     * Store interpolated values for this wavelength

4. **Interpolate in wavelength dimension:**
   
   * For each RT wavelength:
     
     * Interpolate extinction, albedo, asymmetry from Mie wavelengths to RT wavelength using ``interp1_wavelength()``
     * Compute extinction opacity: ``c_out[layer][wavelength] = particle_number_density_cm3 * extinction_cross`` [cm⁻¹]
     * Store albedo and asymmetry: ``a_out[layer][wavelength]``, ``g_out[layer][wavelength]``

**Interpolation Methods:**

* ``interp1_loglog()``: Log-log interpolation for extinction cross-section
* ``interp1_log()``: Linear interpolation in log-radius space for albedo and asymmetry
* ``interp1_wavelength()``: Linear interpolation in wavelength space
* Out-of-range values are clamped to table boundaries automatically

.. _global_variables:

Global Variables Reference
---------------------------

This section provides a quick reference for important global variables declared in ``epacris_main.c`` and used throughout the codebase.

Atmospheric Structure Arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``Tdoub[]``, ``Pdoub[]``, ``MMdoub[]``, ``zdoub[]``, ``TAUdoub[]``: Double-resolution grid arrays for non-isothermal layers (temperature, pressure, number density, altitude, optical depth).
* ``zl[]``, ``pl[]``, ``tl[]``: Altitude (km), pressure (Pa), and temperature (K) at layer centers.
* ``MM[]``: Number density at layer centers (molecules/cm³).
* ``MMZ[]``: Mean molecular mass at layer centers [kg/mol].
* ``meanmolecular[]``: Mean molecular mass for each layer (atomic mass units).
* ``GA``: Gravitational acceleration (m/s²).
* ``scaleheight``: Atmospheric scale height (km) - calculated per layer during hydrostatic equilibrium.
* ``Tnew[]``, ``Pnew[]``: Updated temperature and pressure profiles during climate-chemistry iterations.
* ``mkv[]``: Temperature array used in chemistry calculations.

Radiative Transfer Arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``wavelength[]``: Wavelength grid (nm).
* ``solar[]``: Stellar flux (W/m²/nm).
* ``crossr[]``: Rayleigh scattering cross-sections (cm²).
* ``crossa[3][]``, ``sinab[3][]``, ``asym[3][]``: Aerosol cross-sections, single-scattering albedo, and asymmetry parameter [wavelength] (for up to 3 aerosol types).
* ``opacXXX[][]``: Gas opacity arrays for each molecular species [layer][wavelength] (cm²/molecule). Species include: CO₂, O₂, SO₂, H₂O, OH, H₂CO, H₂O₂, HO₂, H₂S, CO, O₃, CH₄, NH₃, C₂H₂, C₂H₄, C₂H₆, HCN, CH₂O₂, HNO₃, N₂O, N₂, NO, NO₂, OCS, HF, HCl, HBr, HI, ClO, HClO, HBrO, PH₃, CH₃Cl, CH₃Br, DMS, CS₂.
* ``cH2O[][]``, ``aH2O[][]``, ``gH2O[][]``: Cloud optical properties for H₂O [layer][wavelength] (extinction in cm⁻¹, albedo dimensionless, asymmetry dimensionless).
* ``cNH3[][]``, ``aNH3[][]``, ``gNH3[][]``: Cloud optical properties for NH₃.
* ``THETAREF``: Slant path angle for radiative transfer [radians] (converted from ``THETAANGLE`` in degrees).
* ``new_ttop``: Calculated equilibrium temperature at top of atmosphere (K).
* ``CROSSHEADING_STR[]``: String used to build opacity file paths (e.g., "../Opacity/MayTest/opacH2O.dat").
* ``species[]``: Array of opacity species names from ``OPACITY_SPECIES_LIST`` config macro.
* ``NUM_SPECIES``: Number of opacity species loaded.

Collision-Induced Absorption Arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``H2H2CIA[][]``, ``H2HeCIA[][]``, ``H2HCIA[][]``, ``N2H2CIA[][]``, ``N2N2CIA[][]``, ``CO2CO2CIA[][]``: CIA opacities [layer][wavelength] (cm²/molecule).
* ``MeanH2H2CIA[]``, ``MeanH2HeCIA[]``, ``MeanH2HCIA[]``, ``MeanN2H2CIA[]``, ``MeanN2N2CIA[]``, ``MeanCO2CO2CIA[]``: Mean CIA values [layer] (cm²/molecule).
* ``SMeanH2H2CIA[]``, ``SMeanH2HeCIA[]``, ``SMeanH2HCIA[]``, ``SMeanN2H2CIA[]``, ``SMeanN2N2CIA[]``, ``SMeanCO2CO2CIA[]``: Solar-weighted mean CIA values [layer] (cm²/molecule).

Chemistry Arrays
~~~~~~~~~~~~~~~~

* ``xx[][]``: Number density for each species at each layer [layer][species] (molecules/cm³).
* ``clouds[][]``: Cloud particle abundances [layer][species] (molecules/cm³).
* ``ReactionR[][]``, ``ReactionM[][]``, ``ReactionP[][]``, ``ReactionT[][]``: Chemical reaction stoichiometry arrays for kinetic reactions, modified reactions, photolysis reactions, and thermal reactions.
* ``numr``, ``numm``, ``numt``, ``nump``: Counters for number of kinetic, modified, thermal, and photolysis reactions.
* ``numx``, ``numc``, ``numf``, ``numa``: Counters for number of species, condensibles, frozen species, and aerosol species.
* ``waternum``, ``waterx``: Water-related reaction counters.
* ``MeanXXX[]``: Mean mixing ratios [layer] for each species (dimensionless).
* ``SMeanXXX[]``: Solar-weighted mean mixing ratios [layer] for each species (dimensionless).

Cloud Physics Arrays
~~~~~~~~~~~~~~~~~~~~

* ``particle_r0[][]``: Mode radius (nucleation radius) [layer][condensible] (μm).
* ``particle_r1[][]``: Surface-area-weighted radius [layer][condensible] (μm) - stored but not currently used.
* ``particle_r2[][]``: Volume-weighted radius [layer][condensible] (μm) - used for cloud optics interpolation.
* ``particle_VP[][]``: Particle volume [layer][condensible] (cm³) - stored but not currently used.
* ``particle_mass[][]``: Particle mass [layer][condensible] (kg) - stored but not currently used.
* ``particle_number_density[][]``: Particle number density [layer][condensible] (particles/m³) - used for cloud optics.
* ``fall_velocity_ms[][]``: Gravitational settling velocity [layer][condensible] (m/s).
* ``cloud_retention[][]``: Cloud retention factor [layer][condensible] (dimensionless).

Cloud Optics Arrays
~~~~~~~~~~~~~~~~~~~

* ``cloud_mie_optics[]``: Mie scattering lookup tables for each cloud species (structure array).
* ``cloud_species_ids[]``: Array of cloud species IDs corresponding to loaded Mie tables.
* ``n_cloud_species_loaded``: Number of cloud species with loaded Mie tables.

Dynamic Condensibles Management
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``NCONDENSIBLES``: Number of condensible species detected (updated dynamically).
* ``CONDENSIBLES[]``: Array of condensible species IDs [condensible_index].

Radiative Transfer Control Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``RTstepcount``: Radiative transfer step counter (incremented during RT iterations).
* ``rt_drfluxmax_init``: Initial radiative flux maximum for scaling convergence checks.
