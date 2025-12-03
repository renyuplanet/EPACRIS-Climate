/**
 * climate.h
 * 
 * Header file for radiative-convective climate solver
 * Contains function declarations, type definitions, and macros
 */

#ifndef CLIMATE_H
#define CLIMATE_H

#include <stdbool.h>
#include "config.h"

// Forward declarations for types defined in config.h
// (zbin, etc. are defined there)

/**
 * Convergence status structure for radiative transfer
 * Tracks individual convergence criteria and overall status
 */
typedef struct {
    bool flux_converged;        // Rfluxmax converged
    bool gradient_converged;    // dRfluxmax converged
    bool net_flux_converged;    // radiationO converged
    bool overall_converged;     // Overall RT convergence achieved
} RTConvergenceStatus;

/**
 * Helper macro: Print diagnostics every step for jacobian solver (TIME_STEPPING==0),
 * otherwise every PRINT_ITER steps
 */
#define SHOULD_PRINT_DIAG(step_count) \
    ((TIME_STEPPING == 0) ? 1 : ((step_count) % PRINT_ITER == 0))

/**
 * Check radiative transfer convergence criteria
 * 
 * @param Rfluxmax Maximum radiative flux residual
 * @param dRfluxmax Maximum radiative flux gradient
 * @param Tint Interior temperature [K]
 * @param tol_rc Absolute tolerance [W/m²]
 * @param tol_rc_r Relative tolerance factor (multiplied by internal heat flux)
 * @param radiationO Net outgoing flux at TOA [W/m²]
 * @return RTConvergenceStatus structure with convergence flags
 */
RTConvergenceStatus check_rt_convergence(double Rfluxmax, double dRfluxmax, 
                                        double Tint, double tol_rc, double tol_rc_r, 
                                        double radiationO);

/**
 * Main radiative-convective climate solver
 * 
 * Solves for atmospheric temperature profile using iterative radiative transfer
 * and convective adjustment until convergence is achieved.
 * 
 * @param tempeq Output: Equilibrium temperature profile [K]
 * @param P Pressure profile [Pa]
 * @param T Input: Initial temperature profile [K]
 * @param Tint Interior temperature [K]
 * @param outnewtemp Output filename for temperature profile
 * @param outrtdiag Output filename for RT diagnostics
 * @param outrcdiag Output filename for RC diagnostics
 * @param outcondiag Output filename for condensation diagnostics
 * @param nmax_iteration Current NMAX iteration number (for diagnostics)
 */
void ms_Climate(double tempeq[], double P[], double T[], double Tint, 
                char outnewtemp[], char outrtdiag[], char outrcdiag[], 
                char outcondiag[], int nmax_iteration);

#endif // CLIMATE_H

