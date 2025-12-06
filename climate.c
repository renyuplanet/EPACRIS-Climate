// Standard library headers
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>

// Project headers
#include "config.h"
#include "constant.h"
#include "global_temp.h"  // Include global variable declarations
#include "cloud_optics.h"  // Cloud optical property functions
#include "ms_functions.h"  // Contains filleq, fillmi, etc.
#include "climate.h"  // This file's header (contains RTConvergenceStatus, SHOULD_PRINT_DIAG, function declarations)
#include "conv_cond_funcs.h"
#include "readcross.h"
#include "readcia.h"

// Import C files
#include "ludcmp.c"
#include "lubksb.c"
#include "BTridiagonal.c"
#include "ms_radtrans_test.c"
#include "conv_cond_funcs.c"
#include "plot_utils.c"

// ============================================================================
// Function Implementations
// ============================================================================

RTConvergenceStatus check_rt_convergence(double Rfluxmax, double dRfluxmax, 
                                       double Tint, double tol_rc, double tol_rc_r, double radiationO) {
    RTConvergenceStatus status = {0};
    
    // Calculate reference values
    // Internal flux check makes little sense, but keeping it
    double bb_flux = SIGMA * pow(Tint, 4.0);
    
    // Check individual criteria
    status.flux_converged = (Rfluxmax < tol_rc) || (Rfluxmax < tol_rc_r * bb_flux);
    status.gradient_converged = (dRfluxmax < Tol_RC_gradient);
    status.net_flux_converged = (fabs(radiationO) < tol_rc) || (fabs(radiationO) < tol_rc_r * bb_flux);

    // Return converged if Rflux in layers is converged with either gradient_converged or net_flux_converged
    status.overall_converged = (status.flux_converged && status.gradient_converged) || (status.flux_converged && status.net_flux_converged);
    return status;
}

void ms_Climate(double tempeq[], double P[], double T[], double Tint, char outnewtemp[],char outrtdiag[],char outrcdiag[],char outcondiag[], int nmax_iteration)
{
    // Initialize variables
    int i,j,k,radcount,jabove,iconden;
    double tempb[zbin+1], tempbnew[zbin+1], tempc[zbin+1]; //temperature profile arrays
    double znew[zbin+1]; // For printing altitude to files
    double lapse[zbin+1], cp[zbin+1];
    double scaleheight; // For calculating altitude

    // Copy temperature to a new array
    for (j=0; j<=zbin; j++) {
		tempb[j] = T[j]; 
		//tempb[j] = Tdoub[2*j]; //bottom-up! // This line originalyl commented out
	}

    /* Helium calculation 
    The left over mass in the atmosphere is assumed to be helium.
    If the value is negative, it is set to zero.*/ 
    double heliumnumber[zbin+1];
    for (j=1; j<=zbin; j++) {
        heliumnumber[j] = MM[j];
        for (i=1; i<=NSP; i++) {
            heliumnumber[j] -= xx[j][i];
        }
        if (heliumnumber[j]<0.0) {
            heliumnumber[j] = 0.0;
        }
    }
    
    // Condensible detection 
    // Initialize condensibles mode based on CONDENSATION_MODE setting
    initialize_condensibles_mode();
    
    // If using automatic or hybrid mode, detect condensibles across atmosphere
    if (CONDENSATION_MODE == 1 || CONDENSATION_MODE == 2) {
        // Following is run once to detect before loop if TIMING = 0
        detect_condensibles_atmosphere();
    }
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //double Rflux[zbin+1];
    double *Rflux;
    
    // Allocate 2D array for saturation ratios [layer][condensible_species]
    // Use MAX_CONDENSIBLES to handle dynamic changes in NCONDENSIBLES
    double **saturation_ratios;
    saturation_ratios = (double **)malloc((zbin+1) * sizeof(double *));
    for (int j=0; j<=zbin; j++) {
        saturation_ratios[j] = (double *)malloc(MAX_CONDENSIBLES * sizeof(double));
        // Initialize all to zero (including unused slots)
        for (int i=0; i<MAX_CONDENSIBLES; i++) {
            saturation_ratios[j][i] = 0.0;
        }
    }
    
    // Use global particle_number_density array for plotting (already populated by cloud_redistribution_none)
    // No need to allocate - use global array directly
    double Rfluxmax,dRfluxmax;
    int isconv[zbin+1], ncl=0, nrl=zbin; //convective layers setup
    int isconv_old[zbin+1], prevconv[zbin+1];
    int isconv_sum, deltaconv;
    double t_old[zbin+1],tl_old[zbin+1]; //for convective adj test
    int isequil[zbin+1], sumisequil;
    double dt[zbin+2]; //time step if TIME_STEPPING
    char RTstopfile[1024]; 
    
    // Variables for realistic rainout pressure adjustment
    double P_original[zbin+1];  // Store original pressure profile
    double cumulative_mass_loss[zbin+1];  // Track mass loss per layer
    int rainout_events_completed = 0;  // Track rainout events for this NMAX iteration
    for (j=0; j<=zbin; j++) {
        P_original[j] = P[j];
        cumulative_mass_loss[j] = 1.0;  // Initialize to no mass loss
    }
    
    // Radiation flux variables for diagnostics
    double radiationI0, radiationI1, radiationO; // TOA incoming, BOA incoming, TOA outgoing
    
    int ncreg[zbin+1] = {0}; //convective regions
    double pot_temp[zbin+1] = {0.0}; //potential temperature per convective region
    int nclouds[zbin+1] = {0}; //number of condensed species per layer
    strcpy(RTstopfile,OUT_DIR);
    strcat(RTstopfile,"rt.stop");
    //printf("%s %s\n", "RTstopfile = ", RTstopfile);
//atexit(pexit);exit(0); //ms debugging mode

    for (j=1;j<=zbin+1;j++) dt[j] = 1.0e+0; //for TIME_STEPPING: init in some random unit to ...
    for (j=1; j<=zbin; j++) isconv[j] = 0;
//    for (j=1; j<=zbin; j++) ms_adiabat(j,lapse,heliumnumber[j],&cp[j]); //we only need initial cp[] here if dt depends on it

    Rflux = dvector(0,nrl);

//--------------------------------------------------------------------- 
// *** Initialize znew with initial hydrostatic equilibrium ***
//--------------------------------------------------------------------- 
    znew[0] = 0.0;
    for (j=1; j<=zbin; j++) {
        tempc[j] = (tempb[j] + tempb[j-1]) / 2.0;
        //printf("DEBUG j=%d: meanmolecular[%d]=%.6e, tempc[%d]=%.6e\n", j, j, meanmolecular[j], j, tempc[j]);
        //printf("DEBUG GA=%.6e, AMU=%.6e, KBOLTZMANN=%.6e\n", GA, AMU, KBOLTZMANN);
        scaleheight = KBOLTZMANN * tempc[j] / meanmolecular[j] / AMU / GA / 1000;
        //printf("DEBUG scaleheight = %.6e\n", scaleheight);
        znew[j] = znew[j-1] - scaleheight * log(P[j] / P[j-1]);
        //printf("DEBUG znew[%d] = %.6e\n", j, znew[j]);;
    }

//--------------------------------------------------------------------- 
// *** Live plot debugg ***
//--------------------------------------------------------------------- 
    int total_step_count = 0; // For live plot, to keep persistence
#if LIVE_PLOTTING
    // Create live_plot directory path in the output folder
    char live_plot_dir[1024];
    sprintf(live_plot_dir, "%slive_plot/", OUT_DIR);
    
    // Create the live_plot directory if it doesn't exist
    char mkdir_cmd[1100];
    #ifdef _WIN32
        sprintf(mkdir_cmd, "mkdir \"%s\" 2>nul", live_plot_dir);
    #else
        sprintf(mkdir_cmd, "mkdir -p \"%s\" 2>/dev/null", live_plot_dir);
    #endif
    system(mkdir_cmd);
#else
    // Dummy variable to avoid compilation errors when LIVE_PLOTTING is disabled
    char live_plot_dir[1] = {0};
#endif

//*************************************************************
// main iterative Rad-Conv loop:
//*************************************************************
    // For number of radiative-convective iterations (defined in config file)
    for (i=1; i<=NMAX_RC; i++) {
        
        /* store last RC boudary */
        for (j=0; j<=zbin; j++) t_old[j] = tempb[j];
        for (j=1; j<=zbin; j++) {
            //tl_old[j] = tl[j]; //ms23: seems weird to overwrite here
            isconv_old[j] = isconv[j];
            prevconv[j] = isconv[j]; //for convective adjustment iteration
            isconv[j] = 0; //Reset to all radiative for RT (RT expects this)
            // Not resetting isconv to 0 causes the code to crash
        }

        // Resets number of convective layers
        ncl=0;
        // Resets the number of radiative layers (all assumed radiative)
        nrl=zbin;

        // Print total radiative and convective and isconv[j] layers
        // printf("Total radiative layers: %d\n", nrl);
        // printf("Total convective layers: %d\n", ncl);
        // printf("isconv[j] layers: ");
        // for (j=1; j<=zbin; j++) {
        //     printf("%d ", isconv[j]);
        // }
        // printf("\n");




        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //Radiative Transfer iteration
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        int pevery=20, pcount=0;//Added for printing every few steps
        double Tsaved[zbin+1];
        int saveevery=20; //save Told every n steps for dt progression
        int RTsteplimit;

        for (j=0;j<=nrl;j++) Rflux[j] = 0.0; 
        
        RTstepcount=0; //initialized in main
        if (i==1) RTsteplimit = NMAX_RT;
        if (i>1) RTsteplimit = NRT_RC;

        for (j=0;j<=zbin;j++) Tsaved[j] = tempb[j]; //store T for dt progression
        for (j=0;j<=zbin;j++) isequil[j] = 0; //store T for dt progression

        // Write initial profile for live plot before radiative loop
        total_step_count++; // Increment for initial profile plot
        char initial_diag[512];
        sprintf(initial_diag, "RADIATIVE_DIAGNOSTICS: RTstep=%d Rfluxmax=%.3e dRfluxmax=%.3e equilib_layers=%d ncl=%d nrl=%d radiationI0=%.3e radiationI1=%.3e radiationO=%.3e", 
                0, 0.0, 0.0, 0, ncl, nrl, 0.0, 0.0, 0.0);
                
        // Initialize saturation ratios to zero for initial profile
        for (int j=1; j<=zbin; j++) {
            tl[j] = 0.5 * (tempb[j]+tempb[j-1]);
            for (int k=0; k<NCONDENSIBLES; k++) {
                saturation_ratios[j][k] = 0.0;
            }
        }

#if LIVE_PLOTTING
        write_live_plot_data(total_step_count, live_plot_dir, tempb, P, initial_diag, Tint, Tol_RC, Tol_RC_R, cp, lapse, isconv, saturation_ratios, pot_temp, nmax_iteration, particle_number_density);
#endif

        // Main radiative transfer loop using clean convergence checking
        RTConvergenceStatus status = {0};  // Initialize status
        
        while ((!status.overall_converged || RTstepcount == 0) && RTstepcount < RTsteplimit && access(RTstopfile, F_OK) != 0) {
        //while ( (((Rfluxmax >= Tol_RC_R * SIGMA*pow(Tint,4.0))  || (Rfluxmax >= Tol_RC)) || RTstepcount==0 ) && rtstepcount<1000 ) ;}
            RTstepcount +=1;
            pcount +=1;
            sumisequil=0;

            // Reinterpolation of opacities here, before RadTrans is called
            if (RTstepcount % 10 == 0) {
                reinterpolate_all_cia_opacities();
                reinterpolate_all_opacities();
            }
            
            // DYNAMIC CONDENSATION DETECTION - refresh every NRT_RC steps during radiative transfer
            if (CONDENSATION_MODE == 1 || CONDENSATION_MODE == 2) {
                if (CONDENSATION_TIMING == 1 || (CONDENSATION_TIMING == 2 && SHOULD_PRINT_DIAG(total_step_count))) {
                    int old_ncondensibles = NCONDENSIBLES;
                    //printf("\n=== RT STEP %d: DYNAMIC CONDENSATION DETECTION ===\n", RTstepcount);
                    detect_condensibles_atmosphere();
                    if (NCONDENSIBLES != old_ncondensibles) {
                        printf("CONDENSIBLES COUNT CHANGED: %d → %d during RT step %d\n", 
                               old_ncondensibles, NCONDENSIBLES, RTstepcount);
                    }
                }
            }

            ms_RadTrans(Rflux, tempbnew, P, ncl, isconv, lapse, tempb, Tint, cp, dt, isequil, &radiationI0, &radiationI1, &radiationO);
            total_step_count++; //for live plot, keeps persistence
            
            
        // Live plotting during radiative transfer
#if LIVE_PLOTTING
        if ((i == 1 && SHOULD_PRINT_DIAG(RTstepcount)) || (i > 1 && SHOULD_PRINT_DIAG(total_step_count))) {
            char rt_diag[512];
            sprintf(rt_diag, "RADIATIVE_DIAGNOSTICS: RTstep=%d Rfluxmax=%.3e dRfluxmax=%.3e equilib_layers=%d ncl=%d nrl=%d radiationI0=%.3e radiationI1=%.3e radiationO=%.3e", 
                    RTstepcount, Rfluxmax, dRfluxmax, sumisequil, ncl, nrl, radiationI0, radiationI1, radiationO);
            // Saturation ratios will be updated during next convective adjustment
            write_live_plot_data(total_step_count, live_plot_dir, tempbnew, P, rt_diag, Tint, Tol_RC, Tol_RC_R, cp, lapse, isconv, saturation_ratios, pot_temp, nmax_iteration, particle_number_density);
        }
#endif
            
            // Print status
            //printf("%s %d\n", "Radiative loop RTstepCount = ", RTstepcount);
            //printf("%s\t%s\t%s\t%s\t%s\t%s\n","==== i","P[:]","T_old[:]","T_new[:]","Layer_equil","dt[:] ====");
            for (j=zbin; j>=0; j--) {
                //printf("%d\t%e\t%f\t%f\t%d\t%e\n",j,P[j],tempb[j],tempbnew[j],isequil[j],dt[zbin+1-j]);
                sumisequil += isequil[j];
            }
            
            if (RTstepcount%saveevery==1) {//store T every n steps for dt progression
                for (j=0;j<=zbin;j++) {
                    //if(isequil[j]==0){//only advance dt for non-equilibrium layers
                        if (fabs(tempbnew[j] - Tsaved[j]) < 0.5*saveevery*(fabs(tempbnew[j]-tempb[j]))) {
                            dt[zbin+1-j] /=1.5; //as in HELIOS
                        } else {
                            dt[zbin+1-j] *= 1.1; //as in HELIOS
                        }
                    //}
                    Tsaved[j] = tempb[j];
                }
            }


            // Updating temperature profile for boundary and mid-point layers
            for (j=0; j<=zbin; j++) tempb[j]=tempbnew[j];
            for (j=1; j<=zbin; j++) tl[j] = 0.5* (tempb[j]+tempb[j-1]);

            // Calculating max radiative flux and its gradient for convergence check
            Rfluxmax=0.0;
            dRfluxmax=0.0;
            radcount = 0;
            for (k=1; k<=zbin; k++) {
                if (isconv[k] == 0) { //only continue if layer is radiative
                    radcount += 1;
                    Rfluxmax = fmax(Rfluxmax,fabs(Rflux[radcount]));
                    dRfluxmax = fmax(dRfluxmax,fabs((Rflux[radcount-1] - Rflux[radcount])));
                }
            }

            // if (RTstepcount % PRINT_ITER == 0) {
            //     printf("%s\t%s\t%s\t%s\n","Conv_loop","Rad_loop","Residual flux max","dFnet max");
            //     printf("%d\t\t%d\t\t%f\t%f\n", i, RTstepcount, Rfluxmax, dRfluxmax);
            //     printf("%s\n",filleq);
            // }

            // Convergence check!!! Convergence is reached if:
            // 1. The max radiative flux is less than the tolerance times the blackbody flux at the surface temperature
            // 2. The max radiative flux is less than the tolerance
            // 3. The max gradient of the radiative flux is less than 20% of the tolerance
            // 4. The number of radiative steps is greater than the limit
            
            // Check convergence conditions
            status = check_rt_convergence(Rfluxmax, dRfluxmax, Tint, Tol_RC, Tol_RC_R, radiationO);
            
            // DIAGNOSTIC: Print convergence info every step for jacobian solver, otherwise every PRINT_ITER steps
            if ((i == 1 && SHOULD_PRINT_DIAG(RTstepcount)) || (i > 1 && SHOULD_PRINT_DIAG(total_step_count))) {
                double bb_flux = SIGMA * pow(Tint, 4.0);
                printf("Diagnostic from the climate.c file:\n");
                printf("RT CONVERGENCE DIAGNOSTIC - Step %d:\n", RTstepcount);
                printf("  Rfluxmax = %.4e W/m² (target: < %.2e or < %.2e (Tol_RC_R * bb_flux))\n", Rfluxmax, Tol_RC, Tol_RC_R * bb_flux);
                printf("  dRfluxmax = %.4e (target: < %.2e (Tol_RC_gradient))\n", dRfluxmax, Tol_RC_gradient);
                printf("  radiationO = %.4e (target: < %.2e or < %.2e (Tol_RC_R * bb_flux))\n", radiationO, Tol_RC, Tol_RC_R * bb_flux);
                printf("  Converged: flux=%s, gradient=%s, net_flux=%s, OVERALL=%s\n", 
                       status.flux_converged ? "YES" : "NO", 
                       status.gradient_converged ? "YES" : "NO",
                       status.net_flux_converged ? "YES" : "NO",
                       status.overall_converged ? "CONVERGED" : "ITERATING");
                printf("\n");
            }
            
            // Exit if converged or step limit reached
            if ( status.overall_converged || RTstepcount >= RTsteplimit ) {
                /* Generate final plot */
#if LIVE_PLOTTING
                char final_rt_diag[512];
                sprintf(final_rt_diag, "RADIATIVE_DIAGNOSTICS: RTstep=%d Rfluxmax=%.3e dRfluxmax=%.3e equilib_layers=%d ncl=%d nrl=%d radiationI0=%.3e radiationI1=%.3e radiationO=%.3e", 
                        RTstepcount, Rfluxmax, dRfluxmax, sumisequil, ncl, nrl, radiationI0, radiationI1, radiationO);
                // Saturation ratios will be updated during next convective adjustment
                write_live_plot_data(total_step_count, live_plot_dir, tempbnew, P, final_rt_diag, Tint, Tol_RC, Tol_RC_R, cp, lapse, isconv, saturation_ratios, pot_temp, nmax_iteration, particle_number_density);
#endif
                break;
            }

            // Print-out new TP profile every step for jacobian solver, otherwise every pevery steps
            if ((TIME_STEPPING == 0) || (pcount == pevery))
            {
                /* Print-out new TP profile*/
                //ms23: also print diagnostics
                FILE *fp,*frt;
                fp=fopen(outnewtemp,"w");
                frt=fopen(outrtdiag,"w");
                fprintf(frt,"%s\t%8s\t%8s\t%10s\t%10s\t%8s\t%s\n","Layer","Height","log(P)","Temp_bound.","Rflux","isequil","dt");
                for (j=0; j<=zbin; j++) 
                {
                    fprintf(fp, "%f\t%f\t%f\n", znew[j], log10(P[j]), tempb[j]); //temperature file
                    fprintf(frt, "%d\t%8f\t%8f\t%10f\t%10f\t%d\t%e\n",j, znew[j], log10(P[j]), tempb[j], Rflux[j], isequil[j], dt[j]); //RadTrans diagnostics
                }
                fclose(fp);
                fclose(frt);
                //pcount = 0; nneded for convection later
            }
        }



        //============================================================
        //ms22: END of radiative transfer iteration
        //============================================================

        // This is if you dont want to do convective adjustment
        if(access(RTstopfile, F_OK )==0) printf("\n\n%s\n\n","===== RTstopfile found. RT loop aborted! =====");

        // DIAGNOSTIC: Check if radiative transfer converged
        // Note: We still need to check convection even if RT converged
        if (status.overall_converged) {
            double bb_flux = SIGMA * pow(Tint, 4.0);
            printf("\n=== RADIATIVE TRANSFER CONVERGED ===\n");
            printf("RT achieved convergence after %d steps in RC iteration %d\n", RTstepcount, i);
            printf("Final Rfluxmax = %.3e W/m² (target: < %.2e or < %.2e (Tol_RC_R * bb_flux))\n", Rfluxmax, Tol_RC, Tol_RC_R * bb_flux);
            printf("Final dRfluxmax = %.3e W/m² (target: < %.2e (Tol_RC_gradient))\n", dRfluxmax, Tol_RC_gradient);
            printf("Final radiationO = %.3e W/m² (target: < %.2e or < %.2e (Tol_RC_R * bb_flux))\n\n", radiationO, Tol_RC, Tol_RC_R * bb_flux);
            printf("=== PROCEEDING TO CONVECTIVE ADJUSTMENT ===\n");
            // Do NOT break here - we need to check convection and clouds
        }

        // If RT didn't converge, print diagnostic
        if (!status.overall_converged) {
            double bb_flux = SIGMA * pow(Tint, 4.0);
            printf("\n=== RADIATIVE TRANSFER DID NOT CONVERGE ===\n");
            printf("RT stopped after %d steps (limit) in RC iteration %d\n", RTstepcount, i);
            printf("Final Rfluxmax = %.3e W/m² (target: < %.2e or < %.2e)\n", Rfluxmax, Tol_RC, Tol_RC_R * bb_flux);
            printf("Final dRfluxmax = %.3e W/m² (target: < %.2e)\n", dRfluxmax, Tol_RC_gradient);
            printf("Final radiationO = %.3e W/m² (target: < %.2e or < %.2e)\n", radiationO, Tol_RC, Tol_RC_R * bb_flux);
            printf("Proceeding to convective adjustment...\n");
        }




        //============================================================
        // CONVECTION and CONDENSATION BEGINS HERE
        //============================================================

        deltaconv = zbin; // Start with all layers convective
        while (deltaconv > 0) // Continue until no convective layers are found
        {
            // Reset potential temperature and convective region for each layer
            for (j=1; j<=zbin; j++) {
                pot_temp[j] = 0.0; // Reset potential temperature
                ncreg[j] = 0; // Reset convective region
            } 

            // Calculate condensation equilibrium, cold traps, lapse rate
            // heat capacity lapse rate and heat capacity for each layer
            for (j=1; j<=zbin; j++) 
            {
                // 1. Calculate condensation equilibrium and lapse rate
                condensation_and_lapse_rate(j,lapse,heliumnumber[j],&cp[j],saturation_ratios[j]); //calculate condensation equilibrium, lapse rate and heat capacity
                
                // Record number of cloud layers
                nclouds[j] = 0;
                for (k=1;k<=NSP;k++)
                {
                    if(clouds[j][k] > 0.0) nclouds[j] += 1;
                }
            }
            
            // Recalculate heliumnumber after condensation (xx changes when gas condenses)
            // This ensures heliumnumber accounts for material moved from gas (xx) to clouds
            // Use 'k' instead of 'i' to avoid overwriting the RC loop variable
            for (j=1; j<=zbin; j++) {
                heliumnumber[j] = MM[j];
                for (k=1; k<=NSP; k++) {
                    heliumnumber[j] -= xx[j][k];
                }
                if (heliumnumber[j]<0.0) {
                    heliumnumber[j] = 0.0;
                }
            }

            // printf("CLOUD ABUNDANCES:\n");
            // //print cloud abundances
            // for (j=1; j<=zbin; j++) {
            //     for (k=1; k<=NSP; k++) {
            //         if (clouds[j][k] > 0.0) {
            //             printf("cloud[%d][%d] = %.6e\n", j, k, clouds[j][k]);
            //         }
            //     }
            //     printf("\n");
            // }
            

            // Calculate cloud properties
            if (INCLUDE_CLOUD_PHYSICS == 0) {
                // No cloud physics calculation
            } else if (INCLUDE_CLOUD_PHYSICS == 1) {
                cloud_redistribution_none(GA, P); // Calculate physics without redistribution
                // particle_number_density is populated globally, no copy needed
            } else if (INCLUDE_CLOUD_PHYSICS == 2) {
                // exponential_cloud populates global particle_number_density array
                exponential_cloud(GA, P, NULL); // Parameter unused, uses global array
            }
            
            // Calculate cloud optical properties from particle sizes and densities
            // This interpolates Mie scattering data to RT wavelength grid
            // Uses global arrays: clouds, particle_r2, particle_r0, particle_VP, particle_mass, wavelength
            // Writes to global arrays: cH2O, aH2O, gH2O, cNH3, aNH3, gNH3
            // Only calculate if cloud physics is enabled (clouds array will be populated)
            if (INCLUDE_CLOUD_PHYSICS > 0) {
                calculate_cloud_opacity_arrays();
            }

            ncl = 0;
            nrl = 0;

            /* determine convection, and record convective layers */
            ms_conv_check(tempb, P, lapse, isconv, &ncl, &nrl); //ms22: check convective layers

            // // Check if single layers are causing issues with spiking
            // //          printf("%s %d %s %d\n", "ncl", ncl, "nrl", nrl);
            //  // remove single radiative layer between convective layers
            // for (j=2; j<zbin; j++) { //rh->ms 2021
            //         if (isconv[j]==0 && (isconv[j+1]==1 && isconv[j-1]==1)) {
            //             isconv[j]=1;
            //             ncl = ncl+1;
            //             nrl = nrl-1;
            //         }
            //     }

        
            // remove single convective layer between radiative layers
    /*          for (j=2; j<zbin; j++) { //rh->ms 2021
                    if (isconv[j]==1 && (isconv[j+1]==0 && isconv[j-1]==0)) {
                        isconv[j]=0;
                        ncl = ncl-1;
                        nrl = nrl+1;
                    }
                }
    */        
    

            if(ncl>=0) ms_temp_adj(tempb, P, lapse, isconv, cp, ncreg, pot_temp); //ms22: assign convective regimes, adjust temperatures

            
            // UNUSED CODE START
            // // Track cloud mass changes for diagnostic purposes - we're at the end of a convective adjustment step
            // static double prev_total_cloud = -1.0; // Track between iterations, -1 indicates first run
            // // static double prev_total_gas = -1.0;   // Track between iterations
            // // double curr_total_cloud = 0.0;
            // // double curr_total_gas = 0.0;

            // // Calculate total cloud and gas mass
            // // for (j=0; j<=zbin; j++) {
            // //     for (k=0; k<NCONDENSIBLES; k++) {
            // //         int cspec = CONDENSIBLES[k];
            // //         curr_total_cloud += clouds[j][cspec];
            // //         curr_total_gas += xx[j][cspec];
            // //     }
            // // }
            
            // // Report cloud mass changes if this isn't the first run
            // // if (prev_total_cloud >= 0.0) {
            // //     double cloud_change = curr_total_cloud - prev_total_cloud;
            // //     double gas_change = curr_total_gas - prev_total_gas;
                
            // //     if (fabs(cloud_change) > 1.0e-10 || fabs(gas_change) > 1.0e-10) {
            // //         printf("\nCONVECTIVE STEP %d CLOUD DIAGNOSTICS:\n", total_step_count);
            // //         printf("  Total cloud mass:  %.6e → %.6e (change: %.6e)\n", 
            // //                prev_total_cloud, curr_total_cloud, cloud_change);
            // //         printf("  Total gas mass:    %.6e → %.6e (change: %.6e)\n", 
            // //                prev_total_gas, curr_total_gas, gas_change);
            // //         printf("  Total mass change: %.6e (should be near zero)\n", 
            // //                cloud_change + gas_change);
            // //         printf("\n");
            // //     }
            // // }
            
            // // Store current values for next iteration
            // // prev_total_cloud = curr_total_cloud;
            // // prev_total_gas = curr_total_gas;
            //UNUSED CODE END




            // Add live plotting during convection - plot at every step since there are few convective steps
#if LIVE_PLOTTING
            total_step_count++;
            char conv_diag[512];
            sprintf(conv_diag, "RADIATIVE_DIAGNOSTICS: RTstep=%d Rfluxmax=%.3e dRfluxmax=%.3e equilib_layers=%d ncl=%d nrl=%d radiationI0=%.3e radiationI1=%.3e radiationO=%.3e", 
                    RTstepcount, Rfluxmax, dRfluxmax, sumisequil, ncl, nrl, radiationI0, radiationI1, radiationO);
           
            // Saturation ratios are updated by ms_adiabat calls above
            write_live_plot_data(total_step_count, live_plot_dir, tempb, P, conv_diag, Tint, Tol_RC, Tol_RC_R, cp, lapse, isconv, saturation_ratios, pot_temp, nmax_iteration, particle_number_density);
#else
            total_step_count++; // Still increment for consistency
#endif
            //printf("%s %d\n", "Convection step count = ", total_step_count);


            
            //ms22: check if convective regions were fully addressed:
            deltaconv = 0;
            for (j=1;j<=zbin;j++)
            {
                if (abs(prevconv[j]-isconv[j])>0) deltaconv += 1;
                prevconv[j] = isconv[j];
            }
            for (j=1; j<=zbin; j++) tl[j] = 0.5* (tempb[j]+tempb[j-1]);
            //printf("%s %d\n","deltaconv=",deltaconv);

            
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        } // END of convective adjustment iteration
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        // APPLY RAINOUT BASED ON EVENT MODE
        // NOT TESTED!!!!
        // Not sure if this is necessary anymore, but can be modiified to have some kind of realistic removal of material
        int apply_rainout = 0;  // Flag to determine if rainout should be applied
        
        if (RAINOUT_MODE == 0) {
            // no rainout
            apply_rainout = 0;
        } else if (RAINOUT_MODE == 1) {
            // Single event mode: rainout only once at specified iteration
            apply_rainout = (i == RAINOUT_TRIGGER_ITERATION && rainout_events_completed == 0);
        } else if (RAINOUT_MODE == 2) {
            // Multiple events mode: rainout at specified iterations up to max events
            apply_rainout = (i == RAINOUT_TRIGGER_ITERATION && rainout_events_completed < MAX_RAINOUT_EVENTS);
            // Could extend this to have multiple trigger iterations
        }
        
        if (apply_rainout) {
            rainout_events_completed++;

            printf("\n");
            printf("---------------------------------------------------------\n");
            printf("*** CONVECTIVE STEP %d: Rainout EVENT %d ***\n", total_step_count, rainout_events_completed);
            printf("---------------------------------------------------------\n");
            printf("\n");

            for (j=1; j<=zbin; j++) {
                double layer_mass_loss = 1.0;  // Track mass loss for this layer
                simulate_rainout(j, &layer_mass_loss); // Apply rainout physics
                
                // Store mass loss for pressure adjustment
                cumulative_mass_loss[j] *= layer_mass_loss;
            }        

        } else {
            if (RAINOUT_MODE > 0) {
                printf("RC iteration %d: No rainout (events completed: %d/%d)\n", 
                       i, rainout_events_completed, MAX_RAINOUT_EVENTS);
            }
        }

        //--------------------------------------------------------
        // Some code I was trying to use to do pressure adjustment for realistic rainout approach
        // Not used now, but possibly something similar can be implemented in the future
        //--------------------------------------------------------
        // /* PRESSURE ADJUSTMENT for realistic rainout approach */
        // if (PRESSURE_CONSERVATION == 1) {
        //     // Calculate total mass loss across all layers during this convective step
        //     double total_mass_loss = 0.0;
        //     int layers_with_loss = 0;
        //     double max_mass_loss = 0.0;
            
        //     for (j=1; j<=zbin; j++) {
        //         double layer_loss = 1.0 - cumulative_mass_loss[j];
        //         if (layer_loss > 1.0e-6) {
        //             total_mass_loss += layer_loss;
        //             layers_with_loss++;
        //             max_mass_loss = fmax(max_mass_loss, layer_loss);
        //         }
        //     }
            
        //     // Apply pressure adjustments while maintaining hydrostatic equilibrium
        //     if (layers_with_loss > 0 && total_mass_loss > 1.0e-6) {
        //         printf("PRESSURE ADJUSTMENT: %d layers with mass loss, total loss = %.6f, max loss = %.6f\n", 
        //                layers_with_loss, total_mass_loss, max_mass_loss);
                
        //         // Adjust pressure profile from top down to maintain hydrostatic equilibrium
        //         for (j=zbin; j>=1; j--) {
        //             // Apply cumulative mass loss to pressure
        //             P[j] *= cumulative_mass_loss[j];
                    
        //             // Recalculate hydrostatic equilibrium
        //             if (j > 1) {
        //                 // Update pressure at lower boundary based on updated mass and temperature
        //                 tempc[j] = (tempb[j]+tempb[j-1])/2.0;
        //                 scaleheight = KBOLTZMANN * tempc[j] / meanmolecular[j] / AMU / GA / 1000;
                        
        //                 // Maintain hydrostatic consistency
        //                 double dz = znew[j] - znew[j-1];
        //                 P[j-1] = P[j] * exp(dz / scaleheight);
        //             }
        //         }
                
        //         // Update boundary conditions
        //         pl[0] = P[0];  // Surface pressure
        //         for (j=1; j<=zbin; j++) {
        //             pl[j] = 0.5 * (P[j] + P[j-1]);  // Mid-layer pressures
        //         }
                
        //         printf("PRESSURE UPDATE: Surface pressure: %.3e Pa (change: %.3f%%)\n", 
        //                P[0], 100.0*(P[0]/P_original[0] - 1.0));
        //     }
        // }        
        //--------------------------------------------------------
        //--------------------------------------------------------
        
        /* hydrostatic equilibrium */
        //printf("=== Computing hydrostatic equilibrium ===\n");
        znew[0]=0.0;
        for (j=1; j<=zbin; j++) {
            tempc[j] = (tempb[j]+tempb[j-1])/2.0;
            scaleheight = KBOLTZMANN * tempc[j] / meanmolecular[j] / AMU / GA / 1000 ;
            znew[j] = znew[j-1] - scaleheight*log(P[j]/P[j-1]);

            //ms23: double grid
            Tdoub[2*j] = tempb[j];
            Tdoub[2*j-1] = tempc[j];
        }
        Tdoub[0] = tempb[0];

        // printf("=== Convective layer diagnostics ===\n");
        // printf("%s %d %s %d\n", "Number of convective layers", ncl, "Number of radiative layers", nrl);
        /* compare old and new RC boundary */
        isconv_sum=0;
        for (j=1; j<=zbin; j++) {
            isconv_sum += abs(isconv_old[j]-isconv[j]);
        }
        // printf("%s %d\n", "Change in convective layer is", isconv_sum);
        
        // DIAGNOSTIC: Show convection status
        printf("\n--- CONVECTIVE ADJUSTMENT SUMMARY ---\n");
        printf("RC iteration %d: Found %d convective layers, %d radiative layers\n", i, ncl, nrl);
        printf("Convective boundary changes: %d layers\n", isconv_sum);
        
        // Show current condensibles state
        if (CONDENSATION_MODE == 1 || CONDENSATION_MODE == 2) {
            printf("Current condensibles count: %d species [", NCONDENSIBLES);
            for (int k=0; k<NCONDENSIBLES; k++) {
                printf("%d (%s)", CONDENSIBLES[k], get_species_name(CONDENSIBLES[k]));
                if (k < NCONDENSIBLES-1) printf(", ");
            }
            printf("]\n");
        }
        
        // Check if we're stuck in non-converging radiative state
        if (ncl == 0 && !status.overall_converged) {
            printf("WARNING: No convection needed but radiative transfer not converged!\n");
            printf("This suggests NRT_RC=%d may be too small for radiative convergence.\n", NRT_RC);
            printf("Consider increasing NRT_RC or relaxing radiative tolerances.\n");
        }
        printf("----------------------------------------------------------------------\n");

        if ((TIME_STEPPING == 0) || (pcount == PRINT_ITER))
        {
        //ms23: print convection diagnostics
            FILE *frc,*fcon;
            frc=fopen(outrcdiag,"w");
            fprintf(frc,"%s\t%10s\t%8s\t%8s\t%5s\t%s\t%s\t%s\t%10s\t%s\n","Layer","Height","log(P)","cp","lapse","nClouds","isconv","Conv_Region","Pot_Temp","Temp_center");
            fcon=fopen(outcondiag,"w");
            fprintf(fcon,"%s\t%8s\t%s","Layer","log(P)","nClouds");
            // Print headers for current condensibles (NCONDENSIBLES can change during iterations)
            for (k=0;k<NCONDENSIBLES;k++) {
                fprintf(fcon,"\t%s[%d]\t%s[%d]","Molec",CONDENSIBLES[k],"Cloud",CONDENSIBLES[k]);
            }
            fprintf(fcon,"\n");

            for (j=1; j<=zbin; j++) 
            {
                fprintf(frc,"%d\t%10f\t%8f\t%8f\t%8f\t%d\t%d\t%d\t%10e\t%f\n",j,znew[j],log10(P[j]),cp[j],lapse[j],nclouds[j],isconv[j],ncreg[j],pot_temp[j],tempc[j]); //Rad-Conv diagnostics
                fprintf(fcon,"%d\t%f\t%d",j,log10(P[j]),nclouds[j]); //Condensibles diagnostics
                // Print data for current condensibles (NCONDENSIBLES reflects current state)
                for (k=0;k<NCONDENSIBLES;k++) {
                    fprintf(fcon,"\t%.4e\t%.4e",xx[j][CONDENSIBLES[k]]/MM[j],clouds[j][CONDENSIBLES[k]]/MM[j]);
                }
                fprintf(fcon,"\n");
            }
            fclose(frc);
            fclose(fcon);
            pcount = 0;
        }

        // Smooth temperature profile
        // printf("=== Smoothing temperature profile ===\n");
        for (j=1; j<zbin; j++) tempeq[j] = tempb[j-1]/3. + tempb[j]/3. + tempb[j+1]/3.;
        //ms23: different smoothing for TOA and BOA
        tempeq[0] = tempeq[1] + 1.0/3.0 * (tempeq[1] - tempeq[4]); //trend of smoothed T 1 through T 4
        tempeq[zbin] = tempeq[zbin-1] + 1.0/3.0 * (tempeq[zbin-1] - tempeq[zbin-4]); //trend of smoothed T N-1 through T N-4
        //tempeq[0] = tempb[0]/2. + tempb[1]/2.; //surf
        //tempeq[zbin] = tempb[zbin-1]/2. + tempb[zbin]/2.; //TOA
        
        for (j=1; j<=zbin; j++) tl[j] = 0.5* (tempb[j]+tempb[j-1]);

        // LIVE PLOT AFTER TEMPERATURE SMOOTHING (END OF RC ITERATION)
#if LIVE_PLOTTING
        char final_rc_diag[512];
        sprintf(final_rc_diag, "RADIATIVE_DIAGNOSTICS: RTstep=%d Rfluxmax=%.3e dRfluxmax=%.3e equilib_layers=%d ncl=%d nrl=%d radiationI0=%.3e radiationI1=%.3e radiationO=%.3e RC_ITERATION_END: RC=%d isconv_sum=%d", 
                RTstepcount, Rfluxmax, dRfluxmax, sumisequil, ncl, nrl, radiationI0, radiationI1, radiationO, i, isconv_sum);
        write_live_plot_data(total_step_count, live_plot_dir, tempeq, P, final_rc_diag, Tint, Tol_RC, Tol_RC_R, cp, lapse, isconv, saturation_ratios, pot_temp, nmax_iteration, particle_number_density);
#endif

        // DETERMINE IF CONVERGED
        if (isconv_sum == 0 && i!=1 && sumisequil>=zbin) {
            printf("%s\n",fillmi); 
            printf("Climate converged\n");
            printf("%s\n",fillmi); 
            

            break;
        }
    
        if (ncl == 0) {
            printf("%s\n",fillmi); 
            printf("No convective layers found. Rad-Conv-loop done, but not all layers meet selected radiative equilibrium requirement\n");
            printf("%s\n",fillmi); 
            

            break;
        }
//*************************************************************
    }// END main iterative Rad-Conv loop:
//*************************************************************

    /* Print-out new TP profile */
    //ms23: also diagnostic files
    FILE *fp,*frt,*frc,*fcon;
    fp=fopen(outnewtemp,"w");
    // temperature file header
    fprintf(fp, "#Height(km)\tPressure(log10(Pascal))\tTemperature(K)\n");


    frt=fopen(outrtdiag,"w");
    fprintf(frt,"%s\t%8s\t%8s\t%10s\t%10s\t%8s\t%s\n","Layer","Height","log(P)","Temp_bound.","Rflux","isequil","dt");
    frc=fopen(outrcdiag,"w");
    fprintf(frc,"%s\t%10s\t%8s\t%8s\t%8s\t%s\t%s\t%s\t%10s\t%s\n","Layer","Height","log(P)","cp","lapse","nClouds","isconv","Conv_Region","Pot_Temp","Temp_center");
    fcon=fopen(outcondiag,"w");
    fprintf(fcon,"%s\t%8s\t%s","Layer","log(P)","nClouds");
    // Print headers for final condensibles state (NCONDENSIBLES reflects final converged state)
    for (k=0;k<NCONDENSIBLES;k++) fprintf(fcon,"\t%s %3d\t%s %4d","Molec.",CONDENSIBLES[k],"Cloud",CONDENSIBLES[k]); //Condensibles diagnostics
    fprintf(fcon,"\n");

    for (j=0; j<=zbin; j++) 
    {
        fprintf(fp, "%f\t%f\t%f\n", znew[j], log10(P[j]), tempeq[j]); // temperature file output
        fprintf(frt, "%d\t%8f\t%8f\t%10f\t%10f\t%d\t%e\n",j, znew[j], log10(P[j]), tempb[j], Rflux[j], isequil[j], dt[j]); //RadTrans diagnostics
        if (j>0)
        {
            fprintf(frc,"%d\t%10f\t%8f\t%8f\t%5f\t%d\t%d\t%d\t%10e\t%f\n",j,znew[j],log10(P[j]),cp[j],lapse[j],nclouds[j],isconv[j],ncreg[j],pot_temp[j],tempc[j]); //Rad-Conv diagnostics
            fprintf(fcon,"%d\t%f\t%d",j,log10(P[j]),nclouds[j]); //Condensibles diagnostics
            // Print data for final condensibles state (NCONDENSIBLES reflects final converged state)
            for (k=0;k<NCONDENSIBLES;k++) fprintf(fcon,"\t%.4e\t%.4e",xx[j][CONDENSIBLES[k]]/MM[j],clouds[j][CONDENSIBLES[k]]/MM[j]); //Condensibles diagnostics
            fprintf(fcon,"\n");
        }
    }
    fclose(fp);
    fclose(frt);
    fclose(frc);
    fclose(fcon);
    
    // Free allocated memory
    for (int j=0; j<=zbin; j++) {
        free(saturation_ratios[j]);
    }
    free(saturation_ratios);
    // particle_number_density is a global array, no need to free
    
    free_dvector(Rflux,0,nrl);
}
//atexit(pexit);exit(0); //ms debugging mode
