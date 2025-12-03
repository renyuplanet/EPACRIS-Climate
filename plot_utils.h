#ifndef PLOT_UTILS_H
#define PLOT_UTILS_H

#include <stdio.h>
#include "constant.h"
#include "global_temp.h"

// Function to handle live plotting with diagnostics
void write_live_plot_data(int total_step_count, char *live_plot_dir,
                         double *tempb, double *P, char *diagnostic_header,
                         double Tint, double tol_rc, double tol_rc_r,
                         double *cp_array, double *lapse_array, int *isconv_array, double **saturation_ratio_array,
                         double *pot_temp_array, int nmax_iteration, double (*particle_number_density)[MAX_CONDENSIBLES]);

#endif // PLOT_UTILS_H 