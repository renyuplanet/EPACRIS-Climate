#include "plot_utils.h"
#include <math.h>
#include <stdlib.h>

void write_live_plot_data(int total_step_count, char *live_plot_dir,
                         double *tempb, double *P, char *diagnostic_header,
                         double Tint, double tol_rc, double tol_rc_r,
                         double *cp_array, double *lapse_array, int *isconv_array, double **saturation_ratio_array,
                         double *pot_temp_array, int nmax_iteration, double (*particle_number_density)[MAX_CONDENSIBLES]) {
    
    FILE *temp_data = fopen("Tools/temp_data.txt", "w");
    
    // Write header with all species numbers, cloud counterparts, heat capacity, lapse rate, potential temperature, saturation ratios, and particle sizes
    fprintf(temp_data, "# Temperature Pressure");
    for (int ispec=1; ispec<=NSP; ispec++) {
        fprintf(temp_data, " %d", ispec);
    }
    // Add cloud species headers with 'c' prefix
    for (int ispec=1; ispec<=NSP; ispec++) {
        fprintf(temp_data, " c%d", ispec);
    }
    // Add heat capacity header
    fprintf(temp_data, " cp");
    // Add lapse rate header
    fprintf(temp_data, " lapse");
    // Add convective layer flag
    fprintf(temp_data, " isconv");
    // Add potential temperature header
    fprintf(temp_data, " pot_temp");
    // Add saturation ratio headers for each condensible species
    for (int i=0; i<NCONDENSIBLES; i++) {
        fprintf(temp_data, " sat_ratio_%d", CONDENSIBLES[i]);
    }
    // Add particle number density headers for each condensible species
    for (int i=0; i<NCONDENSIBLES; i++) {
        fprintf(temp_data, " particle_number_density_%d", CONDENSIBLES[i]);
    }
    fprintf(temp_data, "\n");
    
    // Add diagnostic header
    fprintf(temp_data, "# %s\n", diagnostic_header);
    fprintf(temp_data, "# CONVERGENCE_TOLERANCES: Tol_RC=%.2e Tol_RC_R=%.2e BBflux=%.3e\n", 
            tol_rc, tol_rc_r, SIGMA * pow(Tint, 4.0));

    for (int j=zbin; j>=1; j--) {  // Skip j=0 (surface boundary) since cp[0] doesn't exist
        fprintf(temp_data, "%f %e", tempb[j], P[j]/1e5);
        
        // First write gas phase abundances
        for (int ispec=1; ispec<=NSP; ispec++) {
            fprintf(temp_data, " %e", xx[j][ispec]/MM[j]);  // Convert to mixing ratio
        }
        
        // Then write cloud abundances
        for (int ispec=1; ispec<=NSP; ispec++) {
            // Most species will have zero clouds
            double cloud_vmr = 0.0;
            
            // Check if this is a condensible species
            for (int k=0; k<NCONDENSIBLES; k++) {
                if (ispec == CONDENSIBLES[k]) {
                    cloud_vmr = clouds[j][ispec]/MM[j];
                    break;
                }
            }
            fprintf(temp_data, " %e", cloud_vmr);
        }
        
        // Add heat capacity data (cp[j] exists for all atmospheric layers j>=1)
        fprintf(temp_data, " %e", cp_array[j]);
        
        // Add lapse rate data (lapse[j] exists for all atmospheric layers j>=1)
        fprintf(temp_data, " %e", lapse_array[j]);
        
        // Add convective layer flag
        fprintf(temp_data, " %d", isconv_array[j]);
        
        // Add potential temperature data (pot_temp[j] exists for all atmospheric layers j>=1)
        fprintf(temp_data, " %e", pot_temp_array[j]);
        
        // Add saturation ratio data for each condensible species
        for (int i=0; i<NCONDENSIBLES; i++) {
            fprintf(temp_data, " %e", saturation_ratio_array[j][i]);
        }
        
        // Add particle number density data for each condensible species (particles/m³)
        for (int i=0; i<NCONDENSIBLES; i++) {
            fprintf(temp_data, " %e", particle_number_density[j][i]);
        }
        
        fprintf(temp_data, "\n");
    }
    fclose(temp_data);
    
    // Execute Python script with step count, NMAX iteration, and species list arguments
    char python_command[2048];
    sprintf(python_command, "python3 Tools/live_plot.py %d %s %d %s", total_step_count, live_plot_dir, nmax_iteration, SPECIES_LIST);
    system(python_command);
} 