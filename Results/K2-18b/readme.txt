Six example cases represent varying parameters relating to clouds/condensation switched on/off:

Example 1: All condensation is disabled. Clouds cannot form and the TP profile is not affected.
#define INCLUDE_CLOUD_PHYSICS 0
#define ENABLE_COLD_TRAP  0
#define FREEZE_CLOUD 0
#define CONDENSATION_MODE 0
#define NCONDENSIBLES_MANUAL 0  //how many potentially condensing species for manual mode
#define CONDENSIBLES_MANUAL (int[]){} //H2O=7; NH3=9; CO=20; CH4=21; CO2=52; H2=53; O2=54; N2=55
#define NRT_RC      50

Example 2: Condensation is switched on, H2O condenses, but cold traps are disabled, resulting in a block of condensed atmosphere.
#define INCLUDE_CLOUD_PHYSICS 0
#define ENABLE_COLD_TRAP  0
#define FREEZE_CLOUD 0
#define CONDENSATION_MODE 1
#define NRT_RC      50

Example 3: Condensation is on and cold traps are enabled, resulting in a cloud profile. Cloud physics/opacities are turned off here.
#define INCLUDE_CLOUD_PHYSICS 0
#define ENABLE_COLD_TRAP  1
#define FREEZE_CLOUD 0
#define CONDENSATION_MODE 1
#define NRT_RC      50

Example 4: This is the standard mode where clouds form and cloud scattering and opacities now affect the TP profile. Though cloud shape evolves drastically through iterations, which reduces its effect on the TP profile.
#define INCLUDE_CLOUD_PHYSICS 1
#define ENABLE_COLD_TRAP  1
#define FREEZE_CLOUD 0
#define CONDENSATION_MODE 1
#define NRT_RC      50

Example 5: Clouds form and cloud opacities are included but the cloud shape is frozen to its initial state. This results in a more severe effect on the TP profile
#define INCLUDE_CLOUD_PHYSICS 1
#define ENABLE_COLD_TRAP  1
#define FREEZE_CLOUD 1
#define CONDENSATION_MODE 1
#define NRT_RC      50

Example 6: Same as Example 5 but the radiative transfer iterations is set to 200 between each convective adjustment, resulting in an even more severe effect from clouds.
#define INCLUDE_CLOUD_PHYSICS 1
#define ENABLE_COLD_TRAP  1
#define FREEZE_CLOUD 0
#define CONDENSATION_MODE 1
#define NRT_RC      200