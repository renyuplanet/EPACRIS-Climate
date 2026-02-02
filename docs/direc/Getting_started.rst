Getting started
===============

This guide provides step-by-step instructions for setting up and running EPACRIS.

Prerequisites
-------------

Before running EPACRIS, ensure you have downloaded the opacities and Mie tables using this `link <https://zenodo.org/records/17888800>`_.

Make sure you have:

* A C compiler (e.g., ``gcc``)
* Opacity and Mie tables directory setup from the above download link
* (If using live plotting) Python 3 with the following libraries: ``matplotlib``, ``numpy``, and ``pandas``

Setting up a simple case
------------------------

The parameters here are for a high albedo K2-18 b case using a H-N-C-O-S network and increased 10x solar metallicity. 


All configuration is done through the ``config.h`` file in the root directory. The config setup is copied to the output location for each run. The following steps explain key parameters for a simple setup. For a complete parameter reference, see the :doc:`Configuration_file` section. For details on how these parameters are used in the code, see the :doc:`Code_structure` section.

**Step 1: Set output directory and run name**

At the top of ``config.h``, set the run name and output directory:

.. code-block:: c

   #define IN_FILE_NAME "K2-18b"           // Run identifier
   #define OUT_DIR "Results/K2-18b/"      // Output directory is created automatically

**Step 2: Choose opacity species**

Set the opacity species list you want to include:

.. code-block:: c

   #define OPACITY_SPECIES_LIST \
       "H2O", "NH3", "CO2", "CO", "CH4", \
       "H2S", "N2", "OH", "C2H6", "CH2O2", \
       "HNO3", "N2O", "SO2", "NO2", "NO", \
       "O2", "O3", "OCS", "HCN", "HO2", \
       "C2H2", "C2H4", "H2CO", "H2O2"

Currently, CIA opacities are hardcoded in the readcia.c file and can be added or removed there.

Set the directory where opacity files are located, these can be downloaded `here <https://zenodo.org/records/17888800>`_.

.. code-block:: c

   #define CROSSHEADING "../Opacity/main_opacities/"

If adding your own line-by-line opacities, ensure opacity files exist for each species in the format ``opac<SPECIES>.dat`` (e.g., ``opacH2O.dat``, ``opacNH3.dat``). The format of opacities is explained in the :doc:`Code_structure` section (see :ref:`Code_structure:opacity_loading`).

**Step 3: Configure cloud physics**

If you want to include cloud microphysics and Mie scattering, enable the ``INCLUDE_CLOUD_PHYSICS`` flag and set the ``USE_EPACRIS_FORMAT`` to ``1`` for the default EPACRIS Mie table format or ``0`` for the LX-Mie format. ``INCLUDE_CLOUD_PHYSICS 1`` is recommended for a self-consistent model. Included molecule IDs in ``CLOUD_SPECIES_LIST`` will tell the code which condensed species to process for cloud microphysics, while ``KZZ`` will impact cloud particle sizes. Note, chemical species in EPACRIS are referenced with integer IDs, which are determined in the species file such as ``Library/SpeciesList/species_HNCSO.dat`` in the ``Standard Number`` column. For details on cloud physics implementation, see the :doc:`Code_structure` section (specifically :ref:`Code_structure:store_cloud_properties` and :ref:`Code_structure:cloud_optics`).

.. code-block:: c

   #define INCLUDE_CLOUD_PHYSICS 1  // 0 = no clouds, 1 = clouds without redistribution, 2 = clouds with redistribution
   #define USE_EPACRIS_FORMAT 1     //  Mie tables format: 0 = LX-Mie (HELIOS) format, 1 = EPACRIS format
   #define CLOUD_MIE_DIRECTORY_EPACRIS "../Opacity/Clouds/EPACRIS_MIE"
   #define CLOUD_SPECIES_LIST 7, 9  // Species IDs: 7=H2O, 9=NH3, 20=CO, 21=CH4, 52=CO2, etc.
   #define KZZ 1.0E+8               // cm²/s (eddy diffusion coefficient); affects cloud particle sizes


Mie folder structure for the EPACRIS format should follow ``EPACRIS_MIE/H2O/``, ``EPACRIS_MIE/NH3/``, etc. with the folders containing the files ``Albedo.dat``, ``Cross.dat``, ``Geo.dat``.

For LX-Mie format, the folder structure should follow ``LXMieOuput/H2O/``, ``LXMieOuput/NH3/`` and should contain output files ``r0.010000.dat``, ``r0.012589.dat``, ``r0.015849.dat``, etc., as described in the  `HELIOS documentation <https://heliosexo.readthedocs.io/en/latest/sections/tutorial.html#including-clouds>`_

**Step 4: Set planet and stellar properties**

Set the usual planet parameters (planet mass, radius, semi-major axis).

.. code-block:: c

   #define MASS_PLANET 5.15401e+25     // kg (planet mass)
   #define RADIUS_PLANET 1.66470e+07   // m (planet radius)
   #define ORBIT 0.1120                // AU (semi-major axis)

Assign the location of the stellar spectrum file and adjust stellar properties. If required, adjust the ``FaintSun`` parameter to reduce the incoming stellar flux, mimicking albedo:

.. code-block:: c

   #define STAR_SPEC "Library/Star/gj876.txt"  // Path to stellar spectrum file
   #define STAR_RADIUS 0.44                    // Solar radius
   #define STAR_TEMP 3457                      // K (stellar effective temperature)
   #define FaintSun 0.3                        // (1 - albedo) factor

The stellar spectrum file should contain two columns: wavelength (nm) and flux (W/m²/nm at 1 AU).

**Step 5: Set the chemistry mode and relevant paths**

Set the chemistry mode to ``IMODE 0`` for chemical equilibrium. Otherwise, a different value can be chosen, such as ``IMODE 1``, to import a predetermined molecular composition from the species list file (see the :doc:`Configuration_file` section for a full explanation of this parameter). The relevant paths for the elemental budget and molecular species list can also be set here. For details on how chemistry is initialized, see the :doc:`Code_structure` section (specifically :ref:`Code_structure:initial_composition`).

.. code-block:: c

   #define IMODE 0
   #define ELE_ABUN "Library/elemental_abundance_files/new_x10Solar.dat" // Elemental budget file
   #define SPECIES_LIST "Library/SpeciesList/species_HNCSO.dat" // Molecular species list

**Step 6: Set advection and internal heat flux**

Don't forget to adjust heat redistribution and slant path angle of incoming flux, as well as the internal heat flux.

.. code-block:: c

   #define FADV 0.25          // Advection factor (0.25 = uniformly distributed, 0.6667 = no advection)
   #define THETAANGLE 0.      // Slant path angle in degrees (60 = global average, 30 = hemispheric)
   #define TINTSET 50.0       // Internal heat flux temperature (K)



Running the model
-----------------

1. Ensure all required data files are in place:

   * Opacity files in ``Opacity/main_opacities/``
   * Mie tables in ``Opacity/Clouds/EPACRIS_MIE/`` or ``Opacity/Clouds/LXMieOuput/``


2. Edit ``config.h`` with your planet parameters (see above).

3. Compile the code:
   
   .. code-block:: bash
   
      gcc epacris_main.c -lm -o epacris

4. Run the executable from the EPACRIS root directory:
   
   .. code-block:: bash
   
      ./epacris

The model will:

* Print initialization information to the terminal
* Create the output directory if it doesn't exist
* Copy ``config.h`` to the output directory as ``config_<IN_FILE_NAME>.txt``
* Write output files to the output directory (see :ref:`Code_structure:output_preparation` for file descriptions)
* Run the climate solver (see :ref:`Code_structure:climate_solver` for details)
* (If enabled) Generate live debugging plots in ``<OUT_DIR>/live_plot/``


Reading output files
--------------------

The model generates several output files in the ``OUT_DIR`` directory (see :ref:`Code_structure:output_preparation` for details on output file setup):

Main output files
~~~~~~~~~~~~~~~~~

* **NewTemperature.dat**: Final temperature-pressure profile

  * Columns: Layer number, Altitude (km), Pressure (log10(Pa)) Temperature (K)
  
* **ConcentrationSTD_T.dat**: Final atmospheric composition

  * Columns: Altitude grids (km), Temperature (K), Pressure (log10(Pa)), Species mixing ratios (one column per species ID)
  
* **Diagnostic_RT.dat**: Radiative transfer diagnostics

  * Contains energy balance information
  
* **Diagnostic_RC.dat**: Radiative-convective diagnostics

  * Contains information about convective adjustment and cloudy layers
  
* **Diagnostic_condens.dat**: Condensation diagnostics (if clouds enabled) [testing]

  * Contains cloud information such as particle sizes and number densities

Configuration file
~~~~~~~~~~~~~~~~~~

* **config_<IN_FILE_NAME>.txt**: Copy of ``config.h`` used for the run

  * Useful for reproducing results and tracking parameter settings

Live plots (if enabled)
~~~~~~~~~~~~~~~~~~~~~~~

If ``LIVE_PLOTTING = 1``, the model generates temperature-pressure profile evolution plots with other diagnostic information in ``<OUT_DIR>/live_plot/``. To adjust live plot formatting, edit the ``Tools/live_plot.py`` script.




For more detailed information about all the configuration parameters, see the :doc:`Configuration_file` section.
For a general description of the code structure and functionality, see the :doc:`Code_structure` section.
