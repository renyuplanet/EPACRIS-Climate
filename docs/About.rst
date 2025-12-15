.. EPACRIS-cloud documentation master file, created by
   sphinx-quickstart on Mon Dec  8 10:31:47 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

EPACRIS-Climate Cloud v0.1 Documentation
========================================


EPACRIS (ExoPlanet Atmospheric Chemistry & Radiative Interaction Simulator) is a one-dimensional atmospheric structure model that solves for the temperature-pressure profile and chemical composition of planetary atmospheres in radiative-convective equilibrium. The model iteratively couples radiative transfer, convective adjustment, equilibrium chemistry, and cloud microphysics to determine the steady-state atmospheric structure and composition, based on fundamental principles and user-specified initial and boundary conditions.

Current version overview
------------------------

EPACRIS-Climate Cloud v0.1 includes newly added multi-species condensation and cloud microphysics, which is incorporated in the radiatiative-transfer routines. The code accounts for wavelength-dependant opacities of gases (line-by-line), including collision-induced, as well as absorption and scattering of condensed cloud particles. Current version supports arbitrary equilibrium chemistry compositions and can simulate both dilute and non-dilute atmospheres, including cases where condensates significantly affect the atmospheric structure and energy balance.

The cloud physics module implements multi-species condensation and cloud microphysics. The moist adiabatic lapse rate follows the formulation of Graham et al. (2021). The model tracks multiple condensible species simultaneously, computing their saturation vapor pressures and partitioning between vapor and condensed phases at each atmospheric layer. Cloud particle sizes are calculated using a log-normal size distribution, with particle properties determined from the balance between condensation growth, gravitational settling, and turbulent mixing. Cloud optical properties are computed via Mie scattering theory using pre-computed lookup tables, ensuring consistency between the particle size distribution used for optical property calculations and the physical properties computed from the cloud microphysics.

Authors
-------

* `Renyu Hu <https://renyuplanet.github.io/>`_ (Jet Propulsion Laboratory, California Institute of Technology)
* Markus Scheucher (Jet Propulsion Laboratory, California Institute of Technology)
* Mantas Zilinskas (Jet Propulsion Laboratory, California Institute of Technology)

Quick user guide
----------------

* To compile: ``gcc epacris_main.c -lm``
* To run: ``./a.out``
* Config file is ``config.h`` in the root folder.

Acknowledgement
---------------

The research was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration (80NM0018D0004).

License
-------

Copyright © 2025, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.

This software may be subject to U.S. export control laws. By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. User has the responsibility to obtain export licenses, or other export authority as may be required before exporting such information to foreign countries or providing access to foreign persons.

Licensed under the Apache License, Version 2.0 (the "Licence");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
`http://www.apache.org/licenses/LICENSE-2.0 <http://www.apache.org/licenses/LICENSE-2.0>`_

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Reporting Issues
----------------

For any issues and bugs please send an e-mail at `renyu.hu@jpl.nasa.gov <mailto:renyu.hu@jpl.nasa.gov>`_, or submit an issue through the Github system.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   direc/Getting_started
   direc/Configuration_file
   direc/Code_structure
   direc/License