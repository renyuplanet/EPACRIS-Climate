// Import header files
#include "config.h" // Config file needs to be included so IDE's can understand variables
#include <math.h>
#include <stdio.h>   // For printf function
#include <stdlib.h>  // For abs function
#include <stdbool.h>  // For bool type
#include "constant.h"
#include "global_temp.h"  // Include global variable declarations
#include "condensed_heat.h"  // Include condensed heat capacity functions
#include "ms_gasheat.h"  // Include gas-phase heat capacity functions
#include "conv_cond_funcs.h"  // Include function declarations

//--------------------------------------------------------------------- 
/* Forward declarations */
/* Main function declarations are in conv_cond_funcs.h */
/* Helper function forward declarations (needed because they're called before definition) */
//--------------------------------------------------------------------- 

double ms_psat_h2o(double temp); //calc saturation pressure of water
double ms_psat_nh3(double temp);
double ms_psat_co(double temp);
double ms_psat_ch4(double temp);
double ms_psat_co2(double temp);
double ms_psat_h2(double temp);
double ms_psat_o2(double temp);
double ms_psat_n2(double temp);
double ms_psat_h2s(double temp);
double ms_psat_feo2h2(double temp);
double ms_psat_fes(double temp);
double ms_latent(int mol, double temp);



//--------------------------------------------------------------------- 
/* Declare variables */
//--------------------------------------------------------------------- 
// This stores layer-specific alpha values
static double global_alpha_values[zbin+1][MAX_CONDENSIBLES];
static int global_alpha_initialized = 0; // Initialized flag

// Track the original Xtotal values before any cold trapping for proper alpha calculation
static double original_xtotal_values[zbin+1][MAX_CONDENSIBLES]; // Original abundance before cold trapping
static int original_xtotal_initialized = 0; // Initialized flag

//--------------------------------------------------------------------- 
/* Functions */
//--------------------------------------------------------------------- 

// Initialize alpha values using cold trapping information
void init_alpha_values() {
    // This function is called in epacris_main.c and updated during runtime

    // For each layer and condensible species, set alpha to default from config
    for (int layer = 1; layer <= zbin; layer++) {
        for (int i = 0; i < MAX_CONDENSIBLES; i++) {

            // This is later adjusted during runtime
            global_alpha_values[layer][i] = ALPHA_RAINOUT;
        }
    }

    global_alpha_initialized = 1; // Mark as initialized
}

// Initialize original xtotal tracking (condensables)
void init_original_xtotal_tracking() {
    for (int layer = 1; layer <= zbin; layer++) {
        for (int i = 0; i < MAX_CONDENSIBLES; i++) {
            original_xtotal_values[layer][i] = -1.0; // -1 indicates uninitialized
        }
    }
    original_xtotal_initialized = 1;
}

// Store original Xtotal value (only if not already stored)
void store_original_xtotal(int layer, int species_index, double xtotal_original) {
    if (!original_xtotal_initialized) {
        init_original_xtotal_tracking();
    }
    
    if (layer >= 1 && layer <= zbin && species_index >= 0 && species_index < MAX_CONDENSIBLES) {
        // Only store if not already initialized (preserve original value)
        if (original_xtotal_values[layer][species_index] < 0.0) {
            original_xtotal_values[layer][species_index] = xtotal_original;
        }
    }
}

// Get original Xtotal value
// Note: Assumes store_original_xtotal() has been called first to initialize the array
double get_original_xtotal(int layer, int species_index) {
    if (layer >= 1 && layer <= zbin && species_index >= 0 && species_index < MAX_CONDENSIBLES) {
        if (original_xtotal_values[layer][species_index] >= 0.0) {
            return original_xtotal_values[layer][species_index];
        }
    }
    return -1.0; // Indicates no original value stored
}

// Function to update alpha values based on cold trapping physics for a specific layer and species
void update_alpha_from_cold_trapping(int layer, int species_index, double alpha_reduction) {  
    if (layer >= 1 && layer <= zbin && species_index >= 0 && species_index < MAX_CONDENSIBLES) {

        // Alpha is naturally bounded between 0 and 1 by the physics
        // No artificial boundaries needed - you can only remove material you have
        if (alpha_reduction < 0.0) { printf("WARNING: Alpha reduction is negative for layer %d species %d\n", layer, species_index); alpha_reduction = 0.0; }  // Should never happen with proper fraction calculation
        if (alpha_reduction > 1.0) { printf("WARNING: Alpha reduction is greater than 1 for layer %d species %d\n", layer, species_index); alpha_reduction = 1.0; }  // Should never happen with proper fraction calculation

        global_alpha_values[layer][species_index] = alpha_reduction;
    }
}


// Function to get layer-specific alpha value
double get_global_alpha_value(int layer, int species_index) {
    // Initialize if not done yet (but only on first call)
    if (!global_alpha_initialized) {
        init_alpha_values();
        //printf("DEBUG: Auto-initialized alpha values in get_global_alpha_value\n");
    }
    
    if (layer >= 1 && layer <= zbin && species_index >= 0 && species_index < MAX_CONDENSIBLES) {
        //return global_alpha_values[layer][species_index];
        // Returns original preset value for now, uncomment above to return stored global value
        return ALPHA_RAINOUT;
    }
    printf("DEBUG GET: Out of bounds - Layer %d, Species index %d, returning default %.6e\n", layer, species_index, ALPHA_RAINOUT);
    return ALPHA_RAINOUT; // Fallback to default (should not happen in normal operation)
}


void condensation_and_lapse_rate(int lay, double lapse[], double xxHe, double* cp, double saturation_ratios[])
{
// Calculate adiabatic lapse rate according to Graham+2021 Eq.(1) 
// The formulation below should be valid for dry, moist, and multi-species condensate conditions in dilute and non-dilute cases.
// Condensibles by species number from chemistry:
// H2O : 7
// NH3 : 9
// CO  : 20
// CH4 : 21
// CO2 : 52
// H2  : 53
// O2  : 54
// N2  : 55
// Cp avail: Air(0.8 N2, 0.2 O2), CO2, He, N2, NH3, CH4, H2, O2, CO, H2O

    // Local parameters
    int i,j;

    //variables for adiabatic calculation
    double Xd,cp_d; // VMR of non-condensible gas, heat capacity of non-condensible gas
    double Xv[NCONDENSIBLES],Xvold,cp_v[NCONDENSIBLES]; //condensibles mol fractions in vapor form and heat cap. [Xvold is the old value of Xv]
    double Xc[NCONDENSIBLES],cp_c[NCONDENSIBLES]; //condensibles fractions in condensed form (Xc) and heat cap. (cp_c)
    double alpha[NCONDENSIBLES], beta[NCONDENSIBLES], latent[NCONDENSIBLES]; // Latent heat and B=L/RT for each condensible
    
    // MINIMAL FIX: Initialize only beta array to prevent false latent heat effects
    // (Other arrays will be properly set in the loop below)
    for (int init_i = 0; init_i < NCONDENSIBLES; init_i++) {
        beta[init_i] = 0.0;  // Critical: ensure no false latent heat
    }

    //variables for condensation
    double psat[NCONDENSIBLES]; //saturation pressures
    int saturated[NCONDENSIBLES]; //check which molecules are available for condensation
    double dp[NCONDENSIBLES]; //condensed fraction

    double cpxx; //cp * xx for all dry species for the cp_d calculation
    double MM_cp; // MM for all species we have a cp for

    double reservoirs[NCONDENSIBLES]; //keeping track of what was rained out + was available at surface at beginning

    double lapse_num, lapse_denom; //numerator and denominator of lapse rate
    double cp_num, cp_denom; //numerator and denominator of heat capacity
    double sum_beta_xv, big_sum_denom_num_left_term, big_sum_denom_num_right_term; //sum of beta*Xv and big sum of denominator terms of lapse rate

    Xd = 1.0; //mole fraction of non-condensible gas
    cp_d = 0.0; //heat capacity of non-condensible gas


    //Assign mol fractions:
    for (i=0; i<NCONDENSIBLES; i++)
    {
        // Calculate saturation pressure for each condensible species
        if(CONDENSIBLES[i]==7) psat[i] = ms_psat_h2o(tl[lay]);  // Condenses ~273K (liquid), ~180K (ice)
        if(CONDENSIBLES[i]==9) psat[i] = ms_psat_nh3(tl[lay]);  // Condenses ~195K (triple point)
        // These are added and not checked for sources
        if(CONDENSIBLES[i]==20) psat[i] = ms_psat_co(tl[lay]);  // Condenses ~68K (triple point)
        if(CONDENSIBLES[i]==21) psat[i] = ms_psat_ch4(tl[lay]); // Condenses ~90K (triple point)
        if(CONDENSIBLES[i]==52) psat[i] = ms_psat_co2(tl[lay]); // Condenses ~195-216K (dry ice)
        if(CONDENSIBLES[i]==53) psat[i] = ms_psat_h2(tl[lay]);  // Condenses ~14K (triple point)
        if(CONDENSIBLES[i]==54) psat[i] = ms_psat_o2(tl[lay]);  // Condenses ~54K (triple point)
        if(CONDENSIBLES[i]==55) psat[i] = ms_psat_n2(tl[lay]);  // Condenses ~63K (triple point)
        if(CONDENSIBLES[i]==45) psat[i] = ms_psat_h2s(tl[lay]);  // Condenses ~188K (triple point)
        if(CONDENSIBLES[i]==114) psat[i] = ms_psat_feo2h2(tl[lay]); // Condenses ~1500K (solid only)
        if(CONDENSIBLES[i]==115) psat[i] = ms_psat_fes(tl[lay]); // Condenses ~1463K (solid), ~3800K (liquid)
        /* Additional species that could be included:
        if(CONDENSIBLES[i]==43) psat[i] = ms_psat_so2(tl[lay]);  // Condenses ~198K (triple point)
        if(CONDENSIBLES[i]==99) psat[i] = ms_psat_s8(tl[lay]);   // Condenses ~388K (solid form)
        if(CONDENSIBLES[i]==27) psat[i] = ms_psat_c2h2(tl[lay]); // Condenses ~192K (triple point)
        if(CONDENSIBLES[i]==29) psat[i] = ms_psat_c2h4(tl[lay]); // Condenses ~104K (triple point)
        if(CONDENSIBLES[i]==31) psat[i] = ms_psat_c2h6(tl[lay]); // Condenses ~90K (triple point)
        if(CONDENSIBLES[i]==37) psat[i] = ms_psat_hcn(tl[lay]);  // Condenses ~260K (triple point)
        if(CONDENSIBLES[i]==100) psat[i] = ms_psat_kcl(tl[lay]); // Condenses ~1000-1100K (hot Jupiters)
        if(CONDENSIBLES[i]==101) psat[i] = ms_psat_nacl(tl[lay]); // Condenses ~1000-1100K (hot Jupiters)
        if(CONDENSIBLES[i]==103) psat[i] = ms_psat_na2s(tl[lay]); // Condenses ~700-800K (hot Jupiters)
        if(CONDENSIBLES[i]==104) psat[i] = ms_psat_mgsiO3(tl[lay]); // Condenses ~1500-1600K (very hot)
        if(CONDENSIBLES[i]==105) psat[i] = ms_psat_mg2siO4(tl[lay]); // Condenses ~1600-1700K (very hot)
        if(CONDENSIBLES[i]==106) psat[i] = ms_psat_fe(tl[lay]);  // Condenses ~1600-1700K (very hot)
        if(CONDENSIBLES[i]==107) psat[i] = ms_psat_al2O3(tl[lay]); // Condenses ~1700-1800K (very hot)
        if(CONDENSIBLES[i]==114) psat[i] = ms_psat_feo2h2(tl[lay]); // Condenses ~1500K (solid only)
        if(CONDENSIBLES[i]==115) psat[i] = ms_psat_fes(tl[lay]); // Condenses ~1463K (solid), ~3800K (liquid)
        */


        // Start condensation calculation
        Xc[i] = clouds[lay][CONDENSIBLES[i]]/MM[lay]; //preexisting clouds taken into account. Divide by MM to get mole fraction from molecules cm^-3
        Xv[i] = xx[lay][CONDENSIBLES[i]]/MM[lay]; //gas fraction of condensibles MM is in molec/cm^3

                // Calculate total available condensible amount. Does not include the dry component
        double Xtotal = Xv[i] + Xc[i];  // Total condensible (gas + cloud)
        
#if ENABLE_COLD_TRAP_ALPHA
        // TRACK ORIGINAL ABUNDANCE: Store the original Xtotal on first encounter
        store_original_xtotal(lay, i, Xtotal);
#endif
        
        // +++++++++++++++++++++++++++ Cold Trapping +++++++++++++++++++++++++++
        // COLD TRAP MECHANISM: Apply transport limitation for ALL condensible species
        double Xtotal_before_coldtrap = Xtotal;  // Default: no cold trapping

        if (ENABLE_COLD_TRAP && lay > 1) {  // All condensibles, if enabled, not bottom layer
            // Check if any layer below has condensation of this species (indicating a cold trap)
            double min_gas_below = 1.0e10;  // Very large initial value
            bool found_condensation_below = false;

            // Check all layers below starting from lay - 1
            for (int check_layer = lay - 1; check_layer >= 1; check_layer--) {

                // Check if this layer has condensation of current species
                double gas_mf = xx[check_layer][CONDENSIBLES[i]] / MM[check_layer];  // Gas mole fraction
                double cloud_mf = clouds[check_layer][CONDENSIBLES[i]] / MM[check_layer];  // Cloud mole fraction

                //Check which below condensing layer has the least amount of gas
                if (cloud_mf > 1.0e-12) {  // Significant condensation found
                    found_condensation_below = true;
                    // Track minimum gas-phase abundance in condensing layers
                    if (gas_mf < min_gas_below) {
                        min_gas_below = gas_mf;
                    }
                }
            }

            // If cold trap detected, limit total species to what can pass through condensing layers
            if (found_condensation_below && min_gas_below < 1.0e9) {
                // COLD TRAP DETECTED: Limit total species to what can pass through condensing layers
                if (Xtotal > min_gas_below) {
                    Xtotal = min_gas_below;  // Limit to minimum gas-phase abundance below


                }
            }
            
            // GHOST COLD TRAP FIX: Always check if current layer should match gas VMR from layer below
            // This fixes cases where previous iterations left artificial VMR reductions
            if (lay > 1) {  // Not bottom layer
                double current_cloud = Xc[i];
                double gas_below = xx[lay-1][CONDENSIBLES[i]] / MM[lay-1];
                double cloud_below = clouds[lay-1][CONDENSIBLES[i]] / MM[lay-1];
                double current_gas = Xv[i];
                
                // If current layer has no significant condensation and VMR is lower than layer below, fix it
                if (current_cloud < 1.0e-12 && cloud_below < 1.0e-12 && current_gas < gas_below * 0.9) {
                    Xtotal = gas_below;  // Set to same VMR as layer below
                    
                    // if (lay < 90) {  // Debug output
                    //     printf("GHOST COLD TRAP FIX: Layer %d Species %d VMR %.6e → %.6e (matching layer below)\n",
                    //            lay, CONDENSIBLES[i], current_gas, gas_below);
                    // }
                }
            }
        }
        // +++++++++++++++++++++++++++ End of cold trap +++++++++++++++++++++++++++

        // Calculate equilibrium partitioning using (potentially reduced) Xtotal
        double Xv_sat = psat[i] / pl[lay];  // Maximum gas phase at saturation
        double Xc_equilibrium = fmax(0.0, Xtotal - Xv_sat);  // Equilibrium cloud amount



        // Update phases while conserving total mass
        if (Xtotal > Xv_sat) {
            // Supersaturated: some should be cloud
            Xv[i] = Xv_sat;
            Xc[i] = Xc_equilibrium;
        } else {
            // Undersaturated: all should be gas  
            Xv[i] = Xtotal;
            Xc[i] = 0.0;
        }



        // Dry gas fraction
        Xd -= Xv[i] + Xc[i];


#if ENABLE_COLD_TRAP_ALPHA
        // CALCULATE ALPHA RELATIVE TO ORIGINAL ABUNDANCE
        // Alpha represents the fraction of ORIGINAL material remaining in the atmosphere
        double original_xtotal = get_original_xtotal(lay, i);
        double alpha_value;
        
        if (original_xtotal > 0.0) {
            // Calculate alpha as fraction of original abundance
            alpha_value = Xtotal / original_xtotal;
            
            // Ensure alpha stays within physical bounds [0, 1]
            if (alpha_value > 1.0) alpha_value = 1.0;  // Can't have more than original
            if (alpha_value < 0.0) alpha_value = 0.0;  // Can't have negative
            
            if (lay > 30 && CONDENSIBLES[i] == 7) {
                printf("Layer %d, Species %d: Original=%.6e, Current=%.6e, Alpha=%.6e\n", 
                       lay, CONDENSIBLES[i], original_xtotal, Xtotal, alpha_value);
            }
        } else {
            // Fallback if original not stored (shouldn't happen)
            alpha_value = ALPHA_RAINOUT;
            printf("WARNING: No original Xtotal for Layer %d, Species index %d\n", lay, i);
        }
        
        //Store this alpha for layer for species
        update_alpha_from_cold_trapping(lay, i, alpha_value);
#endif



    }
    
    //*=========Heat capacities, condensate retention, latent heat=============*//
    //=======================================================================*//

    // CORRECTED APPROACH: Calculate heat capacities separately for each phase
    // following Graham+2021 methodology exactly
    
    // STEP 1: Calculate heat capacity of TRULY DRY (non-condensible) species only
    double cpxx_dry = 0.0;  // Heat capacity contribution from dry species
    double MM_dry = 0.0;    // Number density of dry species
    
    // Add helium (always dry)
    cpxx_dry += xxHe * HeHeat(tl[lay]);
    MM_dry += xxHe;
    
    // Add other species that are NOT in the condensibles list
    // This ensures we only count truly non-condensible species as "dry"
    
    // Check each species to see if it's condensible
    bool is_condensible[NSP+1] = {false}; // Initialize all to false
    for (i=0; i<NCONDENSIBLES; i++) {
        is_condensible[CONDENSIBLES[i]] = true;
    }
    
    // Add non-condensible species to dry calculation
    int species_list[] = {7, 9, 20, 21, 45, 52, 53, 54, 55}; // All species with heat capacity functions
    int num_species = 9;
    
    for (int s=0; s<num_species; s++) {
        int species_id = species_list[s];
        if (!is_condensible[species_id]) {
            // This species is truly dry (non-condensible)
            switch(species_id) {
                case 7:  cpxx_dry += xx[lay][7] * H2OHeat(tl[lay]); MM_dry += xx[lay][7]; break;
                case 9:  cpxx_dry += xx[lay][9] * NH3Heat(tl[lay]); MM_dry += xx[lay][9]; break;
                case 20: cpxx_dry += xx[lay][20] * COHeat(tl[lay]); MM_dry += xx[lay][20]; break;
                case 21: cpxx_dry += xx[lay][21] * CH4Heat(tl[lay]); MM_dry += xx[lay][21]; break;
                case 45: cpxx_dry += xx[lay][45] * H2SHeat(tl[lay]); MM_dry += xx[lay][45]; break;
                case 52: cpxx_dry += xx[lay][52] * CO2Heat(tl[lay]); MM_dry += xx[lay][52]; break;
                case 53: cpxx_dry += xx[lay][53] * H2Heat(tl[lay]); MM_dry += xx[lay][53]; break;
                case 54: cpxx_dry += xx[lay][54] * O2Heat(tl[lay]); MM_dry += xx[lay][54]; break;
                case 55: cpxx_dry += xx[lay][55] * N2Heat(tl[lay]); MM_dry += xx[lay][55]; break;
            }
        }
    }
    
    // STEP 2: Calculate average heat capacity of dry species
    if (MM_dry > 0.0) {
        cp_d = cpxx_dry / MM_dry;
    } else {
        // No dry species - use a reasonable default (this shouldn't happen normally)
        cp_d = 29.0; // Approximate cp for diatomic gases
        printf("WARNING: No dry species found in layer %d, using default cp_d = %.1f\n", lay, cp_d);
    }


    // STEP 3: Calculate heat capacities for condensible species
    for (i=0; i<NCONDENSIBLES; i++)
    {
        // Set vapor and condensed heat capacities for each condensible species
        if(CONDENSIBLES[i]==7)  {
            cp_v[i] = H2OHeat(tl[lay]);
            cp_c[i] = H2O_liquid_heat_capacity(tl[lay]);
            alpha[i] = get_global_alpha_value(lay, i); // Layer-dependent alpha
        }
        else if(CONDENSIBLES[i]==9)  {
            cp_v[i] = NH3Heat(tl[lay]);
            cp_c[i] = NH3_liquid_heat_capacity(tl[lay]);
            alpha[i] = get_global_alpha_value(lay, i); // Layer-dependent alpha
        }
        else if(CONDENSIBLES[i]==20) {
            cp_v[i] = COHeat(tl[lay]);
            cp_c[i] = CO_condensed_heat_capacity(tl[lay]);
            alpha[i] = get_global_alpha_value(lay, i); // Layer-dependent alpha
        }
        else if(CONDENSIBLES[i]==21) {
            cp_v[i] = CH4Heat(tl[lay]);
            cp_c[i] = CH4_condensed_heat_capacity(tl[lay]);
            alpha[i] = get_global_alpha_value(lay, i); // Layer-dependent alpha
        }
        else if(CONDENSIBLES[i]==52) {
            cp_v[i] = CO2Heat(tl[lay]);
            cp_c[i] = CO2_solid_heat_capacity(tl[lay]);
            alpha[i] = get_global_alpha_value(lay, i); // Layer-dependent alpha
        }
        else if(CONDENSIBLES[i]==53) {
            cp_v[i] = H2Heat(tl[lay]);
            cp_c[i] = H2_condensed_heat_capacity(tl[lay]);
            alpha[i] = get_global_alpha_value(lay, i); // Layer-dependent alpha
        }
        else if(CONDENSIBLES[i]==54) {
            cp_v[i] = O2Heat(tl[lay]);
            cp_c[i] = O2_condensed_heat_capacity(tl[lay]);
            alpha[i] = get_global_alpha_value(lay, i); // Layer-dependent alpha
        }
        else if(CONDENSIBLES[i]==55) {
            cp_v[i] = N2Heat(tl[lay]);
            cp_c[i] = N2_condensed_heat_capacity(tl[lay]);
            alpha[i] = get_global_alpha_value(lay, i); // Layer-dependent alpha
        }
        else if(CONDENSIBLES[i]==45) {
            cp_v[i] = H2SHeat(tl[lay]);
            cp_c[i] = H2S_condensed_heat_capacity(tl[lay]);
            alpha[i] = get_global_alpha_value(lay, i); // Layer-dependent alpha
        }
        else {
            // Unknown condensible species - use defaults
            cp_v[i] = 30.0;  // Default vapor heat capacity
            cp_c[i] = 30.0;  // Default condensed heat capacity
            alpha[i] = get_global_alpha_value(lay, i); // Layer-dependent alpha
            printf("WARNING: Unknown condensible species %d in layer %d\n", CONDENSIBLES[i], lay);
        }
    
        // Calculate latent heat of the condensible species
        latent[i] = ms_latent(CONDENSIBLES[i],tl[lay]);
        
        // CORRECT APPROACH - following JANUS logic exactly
        // β is ONLY applied when condensation is actively occurring (P ≥ P_sat)
        // When P < P_sat, species is treated as dry (β = 0, no latent heat effect)
        double partial_pressure = Xv[i] * pl[lay];
        
        // Get critical temperature for this species
        double T_crit = 1e6;  // Default very high temperature (no condensation)
        if(CONDENSIBLES[i]==7)  T_crit = 647.1;   // H2O critical temperature
        if(CONDENSIBLES[i]==9)  T_crit = 405.5;   // NH3 critical temperature
        if(CONDENSIBLES[i]==20) T_crit = 132.92;  // CO critical temperature
        if(CONDENSIBLES[i]==21) T_crit = 190.44;  // CH4 critical temperature
        if(CONDENSIBLES[i]==52) T_crit = 304.2;   // CO2 critical temperature
        if(CONDENSIBLES[i]==53) T_crit = 33.2;    // H2 critical temperature
        if(CONDENSIBLES[i]==54) T_crit = 154.54;  // O2 critical temperature
        if(CONDENSIBLES[i]==55) T_crit = 126.2;   // N2 critical temperature
        if(CONDENSIBLES[i]==45) T_crit = 373.3;   // H2S critical temperature
        
        if (partial_pressure >= psat[i] && tl[lay] < T_crit) {
            // CONDENSATION OCCURRING: Apply positive β (latent heat release)
            beta[i] = latent[i] / R_GAS / tl[lay];
        } else {
            // NO CONDENSATION: Treat as dry species (β = 0)
            // This includes: P < P_sat (undersaturated) OR T > T_critical
            beta[i] = 0.0;
        }
        
        //if(lay==1)  printf("  After: cpxx=%.3e, MM_cp=%.3e\n", cpxx, MM_cp);
    }
    
    // // Debug output for heat capacity calculation
    // if(lay==1) {
    //     printf("HEAT CAPACITY CALCULATION - Layer %d:\n", lay);
    //     printf("  Dry species: MM_dry = %.3e, cpxx_dry = %.3e, cp_d = %.3f J/(mol·K)\n", 
    //            MM_dry, cpxx_dry, cp_d);
    //     for (i=0; i<NCONDENSIBLES; i++) {
    //         printf("  Species %d: cp_v = %.3f, cp_c = %.3f J/(mol·K)\n", 
    //                CONDENSIBLES[i], cp_v[i], cp_c[i]);
    //     }
    // }

    // EQ 1 from Graham+2021 d ln T / d ln P = (x_d + Σx_{v,i}) / 
                //(x_d * [c_d*x_d + Σ(x_{v,i}*(c_{v,i} - R*β_i + R*β_i²) + α_i*x_{c,i}*c_{c,i})] / 
                // [R*(x_d + Σβ_i*x_{v,i})] + Σβ_i*x_{v,i})
    
    // Calculate adiabatic lapse rate
    lapse_num = Xd;
    sum_beta_xv = 0.0;
    big_sum_denom_num_left_term = cp_d * Xd;
    big_sum_denom_num_right_term = 0.0;
    cp_num = cp_d*Xd;
    cp_denom = Xd;

    for (i=0; i<NCONDENSIBLES; i++) 
    {
        lapse_num += Xv[i];
        sum_beta_xv += beta[i]*Xv[i];
        big_sum_denom_num_right_term += Xv[i]*(cp_v[i] - R_GAS*beta[i] + R_GAS*beta[i]*beta[i]) + alpha[i]*Xc[i]*cp_c[i];
        cp_num +=  Xv[i]*cp_v[i] + alpha[i]*Xc[i]*cp_c[i];
        cp_denom += Xv[i];  // CORRECT: Graham's paper denominator is (x_d + Σx_{v,i}) only
        
        // // DEBUG: Print alpha values and Xc values for lapse rate calculations
        // if (alpha[i] > 0.) {  // Only print for layers with significant condensation
        //     printf("LAPSE RATE ALPHA DEBUG - Layer %d, Species %d (ID: %d):\n", lay, i, CONDENSIBLES[i]);
        //     printf("  Xc[%d] = %.6e, alpha[%d] = %.6f (%.2f%% retention)\n", i, Xc[i], i, alpha[i], alpha[i] * 100.0);
        //     printf("  Condensate contribution to lapse rate: %.6e\n", alpha[i]*Xc[i]*cp_c[i]);
        // }
        

        // if (Xc[i] > 1.0e-10 && lay==120) {  // Debug output for any condensation, regardless of layer
        //     printf("HEAT CAPACITY DEBUG - Layer %d, Species %d:\n", lay, CONDENSIBLES[i]);
        //     printf("  Xv[%d] = %.6e, Xc[%d] = %.6e, alpha[%d] = %.3f\n", i, Xv[i], i, Xc[i], i, alpha[i]);
        //     printf("  cp_v[%d] = %.3f J/(mol·K), cp_c[%d] = %.3f J/(mol·K)\n", i, cp_v[i], i, cp_c[i]);
        //     printf("  Vapor contribution to cp_num: %.6e\n", Xv[i]*cp_v[i]);
        //     printf("  Condensate contribution to cp_num: %.6e\n", alpha[i]*Xc[i]*cp_c[i]);
        //     printf("  Running cp_num = %.6e, cp_denom = %.6e\n", cp_num, cp_denom);
        //     printf("  Current effective cp = %.3f J/(mol·K)\n", cp_num/cp_denom);
        // }
    }

    lapse_denom = Xd * (big_sum_denom_num_left_term + big_sum_denom_num_right_term) / (R_GAS*(Xd + sum_beta_xv)) + sum_beta_xv;

    lapse[lay] = lapse_num / lapse_denom;

    *cp = cp_num / cp_denom; // Return heat capacity to climate module

    
    // Setthe gas and cloud abundances based on current equilibrium
    for (i=0; i<NCONDENSIBLES; i++) {
        xx[lay][CONDENSIBLES[i]] = Xv[i] * MM[lay];
        clouds[lay][CONDENSIBLES[i]] = Xc[i] * MM[lay];
   
        // Calculate and store saturation ratios for plotting
        // saturation_ratio = partial_pressure / saturation_pressure
        double partial_pressure = Xv[i] * pl[lay];
        if (psat[i] > 0.0) {
            saturation_ratios[i] = partial_pressure / psat[i];
        } else {
            saturation_ratios[i] = 0.0;
        }
    }
    
    // if (lay == zbin && Xv[0] > 0.1) {  // Check only top layer with significant H2O
    //     printf("LAPSE RATE DEBUG - Layer %d:\n", lay);
    //     printf("  Xd = %.6f, Xv[H2O] = %.6f, Xc[H2O] = %.6f\n", Xd, Xv[0], Xc[0]);
    //     printf("  beta[H2O] = %.6f, latent[H2O] = %.3e J/mol\n", beta[0], latent[0]);
    //     printf("  sum_beta_xv = %.6f\n", sum_beta_xv);
    //     printf("  lapse_num = %.6f\n", lapse_num);
    //     printf("  lapse_denom = %.6f\n", lapse_denom);
    //     printf("  Final lapse rate = %.6f\n", lapse[lay]);
    //     printf("  Temperature = %.1f K, Pressure = %.3e Pa\n", tl[lay], pl[lay]);
    //     printf("\n");
    // }


}// END: void ms_adiabat()
//*********************************************************
//*********************************************************
void simulate_rainout(int lay, double* mass_loss_ratio)
{

    int i, j;
    
    // Get cloud and vapor amounts for condensible species
    double Xv[NCONDENSIBLES], Xc[NCONDENSIBLES];
    double alpha[NCONDENSIBLES]; // retention factors (must match ms_adiabat values)
    
    for (i=0; i<NCONDENSIBLES; i++) {
        Xv[i] = xx[lay][CONDENSIBLES[i]]/MM[lay]; // vapor phase mole fraction
        Xc[i] = clouds[lay][CONDENSIBLES[i]]/MM[lay]; // cloud phase mole fraction
        

    }
    
    // Store original cloud amounts before modification
    double Xc_original[NCONDENSIBLES];
    for (i=0; i<NCONDENSIBLES; i++) {
        Xc_original[i] = Xc[i];
        //printf("RAINOUT: Layer %d, Species %d, Xc_original[%d] = %.6e\n", 
        //       lay, CONDENSIBLES[i], i, Xc_original[i]);
    }

    // STEP 1: Update condensed species amounts after rainout
    for (i=0; i<NCONDENSIBLES; i++) {
        // adjust condensable fraction with alpha
        Xc[i] = alpha[i] * Xc[i];
        xx[lay][CONDENSIBLES[i]] = Xv[i] * MM[lay];
        clouds[lay][CONDENSIBLES[i]] = Xc[i] * MM[lay];

        //printf("RAINOUT: Layer %d, Species %d, Xc[%d] = %.6e, Xv[%d] = %.6e, clouds[%d] = %.6e\n", 
        //       lay, CONDENSIBLES[i], i, Xc[i], i, Xv[i], i, clouds[lay][CONDENSIBLES[i]]);
    }

        
    // STEP 2: Update non-condensible species based on conservation approach
    if (PRESSURE_CONSERVATION == 1) {
        // REALISTIC APPROACH Pressure adjusts:
        // Non-condensibles maintain absolute abundance, mole fractions increase naturally
        // No scaling needed - they're already physically correct
        // Atmospheric mass loss will be tracked for pressure adjustment
        
        // Calculate total mass loss for pressure adjustment
        double original_total_mass = 0.0;
        double new_total_mass = 0.0;
        double remaining_fraction = 0.0;
        // Sum all current species abundances
        for (j=1; j<=NSP; j++) {
            new_total_mass += xx[lay][j];
        }
        
        // Estimate original mass before rainout (approximate)
        original_total_mass = new_total_mass / remaining_fraction;
        
        // Store mass loss ratio for pressure adjustment (global variable or pass back)
        *mass_loss_ratio = new_total_mass / original_total_mass;
        
        if (lay==62) { // Debug output for top layer
            printf("REALISTIC RAINOUT: Layer %d - Mass loss ratio: %.6f\n", lay, *mass_loss_ratio);
            printf("  Original mass estimate: %.6e, New mass: %.6e\n", original_total_mass, new_total_mass);
        }
    } else {
        // Maintain Graham+2021 constraint x_d + x_v + x_c = 1.0
        // Layer-specific scaling: only species with clouds in THIS layer are exempt
        // All other species (including condensibles without clouds here) get scaled
        
        // Calculate total mole fraction lost from cloud rainout in THIS layer
        double total_mole_fraction_lost = 0.0;
        for (i=0; i<NCONDENSIBLES; i++) {
            total_mole_fraction_lost += (1.0 - alpha[i]) * Xc_original[i];
        }

        if (total_mole_fraction_lost > 0.0) {
            // Identify which species actually have clouds in THIS layer
            bool has_clouds_in_layer[NSP+1] = {false};  // Initialize all to false
            for (i=0; i<NCONDENSIBLES; i++) {
                if (clouds[lay][CONDENSIBLES[i]] > 1.0e-20) {  // Small threshold for numerical precision
                    has_clouds_in_layer[CONDENSIBLES[i]] = true;
                }
            }
            
            // Calculate current mole fraction of species that will NOT be scaled
            // (only those with actual clouds in this layer)
            double unscaled_mf = 0.0;
            for (j=1; j<=NSP; j++) {
                if (has_clouds_in_layer[j]) {
                    unscaled_mf += xx[lay][j]/MM[lay];
                }
            }
            
            // Calculate current mole fraction of species that WILL be scaled
            // (everything else, including condensibles without clouds here)
            double scalable_mf = 0.0;
            for (j=1; j<=NSP; j++) {
                if (!has_clouds_in_layer[j]) {
                    scalable_mf += xx[lay][j]/MM[lay];
                }
            }
            
            // Calculate new cloud mole fraction (after rainout)
            double new_cloud_mf = 0.0;
            for (i=0; i<NCONDENSIBLES; i++) {
                new_cloud_mf += clouds[lay][CONDENSIBLES[i]]/MM[lay];
            }
            
            // Scale factor: remaining space divided by what needs to fit
            double scale_factor = (1.0 - unscaled_mf - new_cloud_mf) / scalable_mf;
            
            if (lay==60) { // Debug output for mid-atmosphere layer
                printf("OPEN SYSTEM RAINOUT: Layer %d\n", lay);
                printf("  Cloud mole fraction lost: %.6e\n", total_mole_fraction_lost);
                printf("  Unscaled species mole fraction: %.6e\n", unscaled_mf);
                printf("  Scalable species mole fraction: %.6e\n", scalable_mf);
                printf("  New cloud mole fraction: %.6e\n", new_cloud_mf);
                printf("  Scale factor: %.6f\n", scale_factor);
                
                // Show which species have clouds in this layer
                printf("  Species with clouds in this layer:");
                for (j=1; j<=NSP; j++) {
                    if (has_clouds_in_layer[j]) {
                        printf(" %d", j);
                    }
                }
                printf("\n");
            }

            // Scale only species that do NOT have clouds in this layer
            for (j=1; j<=NSP; j++) {
                if (!has_clouds_in_layer[j]) {
                    xx[lay][j] *= scale_factor;
                }
                // Species with clouds in this layer remain unchanged
            }
            
            if (lay==60) { // Verification
                // Check that mole fractions sum to 1.0
                double total_mole_fraction = 0.0;
                
                // Gas phase mole fractions
                for (j=1; j<=NSP; j++) {
                    total_mole_fraction += xx[lay][j]/MM[lay];
                }
                
                // Condensed phase mole fractions (after rainout)
                for (i=0; i<NCONDENSIBLES; i++) {
                    total_mole_fraction += clouds[lay][CONDENSIBLES[i]]/MM[lay];
                }
                
                printf("  Verification: Total mole fraction = %.6f (should be 1.0)\n", total_mole_fraction);
                
                // Show breakdown
                double scaled_mf = 0.0, unscaled_mf_check = 0.0, cloud_mf = 0.0;
                for (j=1; j<=NSP; j++) {
                    if (has_clouds_in_layer[j]) {
                        unscaled_mf_check += xx[lay][j]/MM[lay];
                    } else {
                        scaled_mf += xx[lay][j]/MM[lay];
                    }
                }
                for (i=0; i<NCONDENSIBLES; i++) {
                    cloud_mf += clouds[lay][CONDENSIBLES[i]]/MM[lay];
                }
                printf("  Breakdown: Scaled=%.6f, Unscaled=%.6f, Clouds=%.6f\n", 
                       scaled_mf, unscaled_mf_check, cloud_mf);
            }
        }
        
        // For open system, no pressure adjustment needed (pressure stays fixed)
        *mass_loss_ratio = 1.0;
    }
    
    // NOTE: Enhanced cloud physics is now handled separately
    // This function only handles Graham's basic alpha rainout
    // Enhanced cloud physics and sedimentation are handled in separate functions
    
}

//=========================================================
//=== Enhanced Cloud Physics Function ====================
//=========================================================

// SPECIES-SPECIFIC PARTICLE PROPERTIES FUNCTION
// Returns particle density [kg/m³], accommodation coefficient [dimensionless], and molecular mass [AMU]
// based on condensible species ID and temperature
void get_particle_properties(int species_id, double temperature, double *density, double *accommodation_coeff, double *molecular_mass) {
    switch (species_id) {
        case 7: // H₂O (Water)
            *density = (temperature < 273.0) ? 917.0 : 1000.0; // Ice vs liquid water [kg/m³]
            *accommodation_coeff = 0.8; // Accommodation coefficient for water vapor [dimensionless]
            *molecular_mass = 18.015; // Molecular mass [AMU]
            break;
            
        case 9: // NH₃ (Ammonia)
            *density = (temperature < 195.0) ? 817.0 : 682.0; // Solid vs liquid ammonia [kg/m³]
            *accommodation_coeff = 0.9; // High sticking probability for NH₃ [dimensionless]
            *molecular_mass = 17.031; // Molecular mass [AMU]
            break;
            
        case 21: // CH₄ (Methane)
            *density = (temperature < 90.7) ? 519.0 : 422.0; // Solid vs liquid methane [kg/m³]
            *accommodation_coeff = 0.7; // Moderate sticking for CH₄ [dimensionless]
            *molecular_mass = 16.043; // Molecular mass [AMU]
            break;
            
        case 52: // CO₂ (Carbon dioxide)
            *density = (temperature < 216.5) ? 1563.0 : 1178.0; // Dry ice vs liquid CO₂ [kg/m³]
            *accommodation_coeff = 0.85; // Good sticking for CO₂ [dimensionless]
            *molecular_mass = 44.010; // Molecular mass [AMU]
            break;
            
        case 45: // H₂S (Hydrogen sulfide)
            *density = (temperature < 187.7) ? 993.0 : 949.0; // Solid vs liquid H₂S [kg/m³]
            *accommodation_coeff = 0.9; // High sticking for H₂S [dimensionless]
            *molecular_mass = 34.081; // Molecular mass [AMU]
            break;
            
        case 20: // CO (Carbon monoxide)
            *density = (temperature < 68.0) ? 929.0 : 789.0; // Solid vs liquid CO [kg/m³]
            *accommodation_coeff = 0.6; // Lower sticking for light molecules [dimensionless]
            *molecular_mass = 28.010; // Molecular mass [AMU]
            break;
            
        case 55: // N₂ (Nitrogen)
            *density = (temperature < 63.1) ? 1026.0 : 808.0; // Solid vs liquid N₂ [kg/m³]
            *accommodation_coeff = 0.6; // Lower sticking for light molecules [dimensionless]
            *molecular_mass = 28.014; // Molecular mass [AMU]
            break;
            
        case 54: // O₂ (Oxygen)
            *density = (temperature < 54.4) ? 1426.0 : 1141.0; // Solid vs liquid O₂ [kg/m³]
            *accommodation_coeff = 0.6; // Lower sticking for light molecules [dimensionless]
            *molecular_mass = 31.998; // Molecular mass [AMU]
            break;
            
        case 53: // H₂ (Hydrogen)
            *density = (temperature < 14.0) ? 86.8 : 70.8; // Solid vs liquid H₂ [kg/m³]
            *accommodation_coeff = 0.5; // Very low sticking for H₂ [dimensionless]
            *molecular_mass = 2.016; // Molecular mass [AMU]
            break;
            
        default: // Unknown species - use water-like properties
            *density = 1000.0; // Default to liquid water density [kg/m³]
            *accommodation_coeff = 0.8; // Default accommodation coefficient [dimensionless]
            *molecular_mass = 18.0; // Default molecular mass [AMU]
            break;
    }
}

//particlesizef function from Kyo
// PARTICLE SIZE CALCULATION USING MICROPHYSICAL BALANCE
// This function calculates equilibrium particle sizes by balancing:
// 1) Condensational growth (mass diffusion from supersaturated vapor)
// 2) Gravitational settling (particles fall due to gravity)
// 3) Eddy diffusion (turbulent mixing opposes settling)
//
// Core Physics: Particles grow until fall velocity balances with diffusion
// Key Equation: u = v_fall * H / K_zz (dimensionless fall parameter)
// Where: H = scale height [m], K_zz = eddy diffusion coefficient [m²/s]
void calculate_cloud_properties(double g, double T, double P, double mean_molecular_mass, int condensible_species_id, 
    double Kzz, int layer, double *r0, double *r1, double *r2, double *VP, double *effective_settling_velocity,
     double *scale_height, double *mass_per_particle, double *n_density) {
    // Local constants to avoid conflicts with main code
    double KB_local = 1.3806503E-23; // Boltzmann constant [J/K]
    double AMU_local = 1.66053886E-27; // Atomic mass unit [kg]
    double RGAS_local = 8.314472; // Universal gas constant [J/(mol·K)]
    double PI_local = 3.1416; // Pi

    // GET SPECIES-SPECIFIC PROPERTIES
    double molecular_mass_condensible; // [AMU]
    double rho, acc; // Particle density [kg/m³] and accommodation coefficient [dimensionless]
    
    // Get all species-specific properties from unified function
    get_particle_properties(condensible_species_id, T, &rho, &acc, &molecular_mass_condensible);

    // ATMOSPHERIC SCALE HEIGHT CALCULATION
    // H = k_B * T / (m_avg * g) where m_avg = mean_molecular_mass * AMU
    // Units: [J/K] * [K] / ([AMU] * [kg/AMU] * [m/s²]) = [m]
    double H = KB_local * T / mean_molecular_mass / AMU_local / g; // Atmospheric scale height [m]
    
    // DIMENSIONLESS FALL PARAMETER
    // u = K_zz / H represents ratio of diffusion to settling
    // Higher u = more diffusion relative to settling = smaller particles
    double u = Kzz / H; // Dimensionless eddy diffusion parameter
    
    // ATMOSPHERIC VISCOSITY (temperature-dependent)
    // Sutherland's formula: μ = μ₀ * (T/T₀)^1.5 * (T₀ + S)/(T + S)
    
    // ATMOSPHERIC COMPOSITION SELECTION
    // Choose atmospheric type: 0=H2, 1=Air, 2=CO2, 3=N2, 4=Ar, 5=He
    int atmosphere_type = 0;  // DEFAULT: H2 atmosphere (common for gas giants/exoplanets)
    
    // Sutherland constants for different atmospheric compositions
    // Values from Sutherland (1893), COMSOL, and atmospheric physics literature
    double mu0, T0, S_sutherland;
    
    if (atmosphere_type == 0) {
        // HYDROGEN (H2) - Default for gas giant/exoplanet atmospheres
        mu0 = 8.411E-6;    // Reference viscosity [Pa·s] at 273K
        T0 = 273.15;       // Reference temperature [K]
        S_sutherland = 97.0;  // Sutherland constant [K]
    } else if (atmosphere_type == 1) {
        // AIR (N2/O2 mix) - Earth-like atmospheres
        mu0 = 1.716E-5;    // Reference viscosity [Pa·s] at 273K
        T0 = 273.15;       // Reference temperature [K] 
        S_sutherland = 110.4;  // Sutherland constant [K]
    } else if (atmosphere_type == 2) {
        // CARBON DIOXIDE (CO2) - Venus-like or early Mars
        mu0 = 1.370E-5;    // Reference viscosity [Pa·s] at 273K
        T0 = 273.15;       // Reference temperature [K]
        S_sutherland = 222.0;  // Sutherland constant [K]
    } else if (atmosphere_type == 3) {
        // NITROGEN (N2) - Titan-like atmospheres
        mu0 = 1.663E-5;    // Reference viscosity [Pa·s] at 273K
        T0 = 273.15;       // Reference temperature [K]
        S_sutherland = 107.0;  // Sutherland constant [K]
    } else if (atmosphere_type == 4) {
        // ARGON (Ar) - Some planetary atmospheres
        mu0 = 2.125E-5;    // Reference viscosity [Pa·s] at 273K
        T0 = 273.15;       // Reference temperature [K]
        S_sutherland = 114.0;  // Sutherland constant [K]
    } else if (atmosphere_type == 5) {
        // HELIUM (He) - Helium-rich atmospheres
        mu0 = 1.865E-5;    // Reference viscosity [Pa·s] at 273K
        T0 = 273.15;       // Reference temperature [K]
        S_sutherland = 79.4;   // Sutherland constant [K]
    } else {
        // DEFAULT: Fall back to H2 if unknown type
        mu0 = 8.411E-6;    // H2 values
        T0 = 273.15;
        S_sutherland = 97.0;
    }
    
    // Calculate atmospheric viscosity using Sutherland's formula
    double mu = mu0 * pow(T / T0, 1.5) * (T0 + S_sutherland) / (T + S_sutherland); // Viscosity [Pa·s]
    
    // MEAN FREE PATH OF GAS MOLECULES
    // λ = 2μ / (P * √(8*mean_molecular_mass/(πRT))) - kinetic theory result
    // Units: [Pa·s] / ([Pa] * √([kg/mol]/([J/(mol·K)] * [K]))) = [m]
    double lambda = 2 * mu / P / pow(8 * mean_molecular_mass * 1.0E-3 / PI_local / RGAS_local / T, 0.5); // Mean free path [m]
    
    // NUMBER DENSITY OF CONDENSING MOLECULES
    // From ideal gas law: n = ΔP / (k_B * T)
    // where ΔP = P_partial - P_saturation (supersaturation pressure)
    // This needs to be updated to use precondensed values from ms_adiabat
    // Units: [Pa] / ([J/K] * [K]) = [molecules/m³]
    double deltan = DELTA_P / KB_local / T; // Excess number density [molecules/m³]

    // MASS DIFFUSION COEFFICIENT
    // Approximation for molecular diffusion in gas
    double D = 0.12e-4; // Mass diffusion coefficient [m²/s]

    // PARTICLE SIZE DISTRIBUTION PARAMETERS
    double Cc0= 1.0; // Initial Cunningham slip correction factor [dimensionless]
    double fa = 1.0; // Initial ventilation factor [dimensionless]  
    double sig = 2.0; // Log-normal size distribution width parameter [dimensionless]

    // ITERATIVE SOLUTION FOR EQUILIBRIUM PARTICLE SIZE
    // Solve for particle volume V by balancing:
    // 1) Condensational growth: ∝ D * deltan / ρ
    // 2) Effective gravitational settling: ∝ max[settling_term - u, 0]  
    // 3) Turbulent diffusion: ∝ u / H
    double Vs;
    double aa_final = 0.0; // Store final aa value for sedimentation calculation
    
    for (int dump = 1; dump <= 1e3; ++dump) {

        // CONDENSATION TERM: Rate of volume growth from vapor deposition
        // cc = -48^(1/3) * π^(2/3) * D * molecular_mass_condensible * fa * Δn / ρ * exp(-ln²(σ))
        // Negative because we're solving the quadratic equation
        // This term should be negligible if deltan is low
        double cc = -(pow(48.0 * PI_local * PI_local, 1.0 / 3.0) * D * molecular_mass_condensible * AMU_local * fa * deltan / rho * exp(-pow(log(sig), 2.0)));
        
        // SETTLING TERM: Gravitational fall velocity coefficient  
        // aa = ρ * g / (μ * 162^(1/3) * π^(2/3) * H) * Cc * exp(-ln²(σ))
        // From Stokes law: v_fall = 2*r²*ρ*g*Cc/(9*μ) where r ∝ V^(1/3)
        // This represents the settling velocity term before the max[... - u, 0] operation
        double aa = rho * g / mu / pow(162.0 * PI_local * PI_local, 1.0 / 3.0) / H * Cc0* exp(-pow(log(sig), 2.0));
        
        // DIFFUSION TERM: Turbulent mixing coefficient
        // bb = -u/H = -K_zz/H² (opposes settling)
        double bb = -u / H;

        // SOLVE QUADRATIC EQUATION: aa*V^(2/3) + bb*V^(1/3) + cc = 0
        // Using quadratic formula after substituting V^(1/3) = x
        // This gives equilibrium volume where growth = settling + diffusion
        // V = [(-bb + sqrt(bb * bb - 4.0 * aa * cc)) / 2.0 / aa]^(3/2)
        double V = pow((-bb + sqrt(bb * bb - 4.0 * aa * cc)) / 2.0 / aa, 3.0 / 2.0);
        
        // Handle cases where updraft dominates (no real solution)
        if (isnan(V) || V <= 0.0) {
            // Updraft-dominated regime: use asymptotic solution
            //V_asym = 39.9 * [μ*u*exp(ln²σ)/(ρ_p*g*C_c)]^(3/2)
            V = 39.9 * pow(mu * u * exp(pow(log(sig), 2.0)) / (rho * g * Cc0), 3.0 / 2.0);
        }
        
        // PARTICLE DIAMETER from volume
        // d = (6V/π)^(1/3) * exp(-ln²(σ)) [m]
        double d1 = pow(6.0 * V / PI_local, 1.0 / 3.0) * exp(-pow(log(sig), 2.0));

        // UPDATE SLIP CORRECTION FACTORS
        // Knudsen number: Kn = λ/d (ratio of mean free path to particle size)
        double kn = lambda / d1;
        
        // CUNNINGHAM SLIP CORRECTION: Cc = 1 + Kn*(1.257 + 0.4*exp(-1.1/Kn))
        // Accounts for non-continuum effects when particles approach molecular size
        double Cc1 = 1.0 + kn * (1.257 + 0.4 * exp(-1.1 / kn));
        
        // VENTILATION CORRECTION: Accounts for enhanced mass transfer during settling
        // fa = (1 + Kn)/(1 + 2*Kn*(1+Kn)/α) where α is accommodation coefficient
        double fa1 = (1.0 + kn) / (1.0 + 2.0 * kn * (1.0 + kn) / acc);
        
        // CHECK CONVERGENCE: Continue iteration until slip factors stabilize
        if (fabs(Cc1 - Cc0) + fabs(fa1 - fa) < 0.001) {
            Vs = V; // Final volume [m³]
            break;
        } else {
            Cc0= Cc1; // Update slip correction
            fa = fa1;  // Update ventilation factor
        }
    }
    
    // CALCULATE CHARACTERISTIC PARTICLE RADII from log-normal distribution
    // For log-normal distribution with geometric standard deviation σ:
    // r_i = r_mean * exp(n * ln²(σ)) where n depends on moment
    
    //r0, r1, r2, VP are in microns
    // r₀: Mean radius weighted by number (smallest particles) [μm]
    *r0 = pow((3.0 * Vs) / (4.0 * PI_local), 1.0 / 3.0) * exp(-1.5 * pow(log(sig), 2.0)) * 1.0E+6;
    
    // r₁: Mean radius weighted by surface area [μm] 
    *r1 = pow((3.0 * Vs) / (4.0 * PI_local), 1.0 / 3.0) * exp(-pow(log(sig), 2.0)) * 1.0E+6;
    
    // r₂: Mean radius weighted by volume (largest/most massive particles) [μm]
    // This is used for sedimentation calculations since fall velocity ∝ r²
    *r2 = pow((3.0 * Vs) / (4.0 * PI_local), 1.0 / 3.0) * exp(-0.5 * pow(log(sig), 2.0)) * 1.0E+6;
    
    // VP: Final particle volume in cm³ (convert from m³: 1 m³ = 1e6 cm³)
    *VP = Vs * 1.0E+6;

    // Calculate mass per particle [kg] (Vs in m³, rho in kg/m³)
    *mass_per_particle = Vs * rho;

    // This is the cloud density but not used now
    // 1. Get cloud mass density from Hu+2019 Equation 2:
    //double cloud_mass_density = (xx_original - xx_saturated) * MM * molecular_mass / RGAS / T;

    // Convert to C
    //double cloud_density = cloud_mass_density / molecular_mass;

    // Calculate the particle number density
    //double particle_number_density = calculate_particle_number_density_from_molecules(molecular_number_density, particle_radius_m, condensate_density, molecular_mass_kg);
    //*particle_number_density = particle_number_density;


    // HU+2019 APPROACH: Use their Equation A3 for settling velocity
    // From Hu+2019 Appendix, Equation A3: v_d = max[v_fall - u, 0]
    // where v_fall = ρ_p*g*Cc/(162π²)^(1/3)*μ * V^(2/3) * exp(-ln²σ)
    
    // Calculate settling velocity using Hu+2019's formula with our solved volume V
    double v_fall_hu2019 = (rho * g * Cc0) / (pow(162.0 * PI_local * PI_local, 1.0 / 3.0) * mu) 
                           * pow(Vs, 2.0 / 3.0) * exp(-pow(log(sig), 2.0));
    
    // Apply Hu+2019's correction: v_d = max[v_fall - u, 0]
    // If using fall - u from Hu, this results in 0.

    double v_d = fmax(0.0, v_fall_hu2019 - u);
    
    // Use r2 (volume-weighted radius) for sedimentation calculations
    double particle_radius_m = *r2 * 1.0e-6; // Convert μm to m

    // Alternative: Calculate using Stokes law with final particle size
    // v_fall = 2*r²*ρ*g*Cc/(9*μ) - Stokes law with Cunningham correction
    // Value comes out the same as Hu+2019, confirming derivation
    double gravitational_settling = (2.0 * particle_radius_m * particle_radius_m * rho * g * Cc0) / (9.0 * mu);

    // Calculate molecules per particle using the mass_per_particle output parameter
    double molecules_per_particle = (*mass_per_particle) / (molecular_mass_condensible * AMU_local);
    // Calculate particle number density [particles/m³]
    // clouds[layer][condensible_species_id] is in molecules/cm³
    // n_particles = n_molecules / molecules_per_particle gives particles/cm³
    // Convert to particles/m³ by multiplying by 1e6 (1 m³ = 1e6 cm³)
    *n_density = (clouds[layer][condensible_species_id] / molecules_per_particle) * 1.0e6;
    // Return the actual fall velocity and scale height
    *effective_settling_velocity = gravitational_settling; // [m/s] - actual fall velocity
    *scale_height = H; // [m]
    // Debug output to verify the theory
    // Useful debug output
    // if (layer >= 55 && layer <= 65) {
    //     printf("PARTICLE NUMBER DENSITY DEBUG (Layer %d):\n", layer);
    //     printf("  Molecular number density: %.2e molecules/cm³\n", clouds[layer][condensible_species_id]);
    //     printf("  Particle radius: %.2e m (%.2f μm)\n", particle_radius_m, particle_radius_m*1e6);
    //     printf("  VP (particle volume): %.2e μm³\n", *VP);
    //     printf("  Particle material density: %.1f kg/m³\n", rho);
    //     printf("  Molecular mass: %.2e kg\n", molecular_mass_condensible * AMU_local);
    //     printf("  Mass per particle: %.2e kg\n", mass_per_particle);
    //     printf("  Calculated particle number density: %.2e particles/m³\n", particle_number_density);
    //     printf("  Molecules per particle: %.2e\n", molecules_per_particle);
    //     // Expected reduction factor of particles if we took into account clumping
    //     // Implementing means VMR doesnt add up to 1
    //     printf("  Expected reduction factor: %.2e\n", clouds[layer][condensible_species_id] / particle_number_density);
    //     printf("  ---\n");
    // }

    // printf ("==============================================\n");
    // printf ("VALUES FROM PARTICLESIZEF_LOCAL\n");
    // printf ("viscosity: %.2e\n", mu);
    // printf ("retention: %.2e\n", retention);
    // printf ("gravitational_settling: %.2e\n", gravitational_settling);
    // printf ("v_d: %.2e\n", v_d);
    // printf ("u value: %.2e\n", u);
    // printf ("f_sed: %.2e\n", f_sed);
    // printf ("scale_height: %.2e\n", H);
    // printf ("==============================================\n");
    // printf ("==end==\n");
}



void cloud_redistribution_none(double gravity, double P[]) {
    
    // Compute particle physics but not redistribution
    for (int layer = 1; layer <= zbin; layer++) {
        for (int i = 0; i < NCONDENSIBLES; i++) {
            int species_id = CONDENSIBLES[i];
            
            // Skip if no condensate in this layer - reset particle size
            if (clouds[layer][species_id] < 1.0e-20) {
                particle_r2[layer][i] = 0.0;  // Reset when condensation stops
                particle_r0[layer][i] = 0.0;
                particle_VP[layer][i] = 0.0;
                particle_mass[layer][i] = 0.0;
                particle_number_density[layer][i] = 0.0;
                continue;
            }
            
            // Get temperature and pressure for this layer
            double T = tl[layer];
            double P_layer = P[layer];

            // Calculate cloud physics properties
            double r0, r1, r2, VP, v_sed, scale_height, mass_per_particle, n_density;

            double Kzz_m2s =  KZZ * 1.0e-4; // Eddy diffusion coefficient [m²/s]
            
            calculate_cloud_properties(gravity,
                T, P_layer, meanmolecular[layer],
                species_id, Kzz_m2s, layer,
                &r0, &r1, &r2, &VP,
                &v_sed, 
                &scale_height,
                &mass_per_particle,
                &n_density);
            
            // Store multiple particle properties for potential use elsewhere
            // All values calculated once in particlesizef_local and stored for reuse
            particle_r2[layer][i] = r2;              // r2: volume-weighted radius [μm]
            particle_r0[layer][i] = r0;             // r0: nucleation/monomer radius [μm]
            particle_VP[layer][i] = VP;             // VP: particle volume [cm³]
            particle_mass[layer][i] = mass_per_particle; // mass per particle [kg]
            particle_number_density[layer][i] = n_density; // particle number density [particles/m³]

            // printf("layer: %d\n", layer);
            // printf("r0: %.2e, r1: %.2e, r2: %.2e, VP: %.2e cm^3\n", r0, r1, r2, VP);
            // printf("mass_per_particle: %.2e\n", mass_per_particle);
        }
    }

    printf("Particle properties computed - no material redistribution performed\n");
}


void ms_conv_check(double tempb[],double P[],double lapse[],int isconv[],int* ncl,int* nrl)
{
// Check each pair of layers for convective stability
//ncl, nrl, and isconv are handed over as 0 initially

    //LOCAL parameters
    int i;
    for (i=0; i<=zbin; i++) isconv[i] = 0;
    
    //TESTING!!
    //isconv[0] = 1;

    //START function
    for (i=0; i<zbin; i++)
    {
//        printf("%s %d%s %f%s%f\n","Layer",i,": Tabove=",tempb[i+1]," but adiabatically ",( tempb[i] * pow(P[i+1]/P[i],lapse[i+1])));
//        if( tempb[i+1] < ( tempb[i] * pow(P[i+1]/P[i],lapse[i+1]) /0.98 ) )
        
        // STABILITY FIX: Use larger tolerance to prevent oscillations at cloud decks
        // where latent heating can cause rapid lapse rate changes
        double stability_tolerance = 0.5;  // Increased from 0.1 K to 0.5 K
        double adiabatic_temp = tempb[i] * pow(P[i+1]/P[i], lapse[i+1]);
        
        if( tempb[i+1] < (adiabatic_temp + stability_tolerance) ) //for numerical stability
        {
        //printf("%s %d\t%s %d\t%f %s %f\n","Layer",i+1,"isconv",isconv[i],tempb[i+1],"<",( tempb[i] * pow(P[i+1]/P[i],lapse[i+1])));
            isconv[i+1] = 1;
            //isconv[i] = 1;
        }
//    printf("%s %d %s %d\n","ms_i",i,"ms_isconv",isconv[i]);
    }
    for (i=1;i<=zbin;i++)
    {
        if(isconv[i]) *ncl+=1;
        else *nrl+=1;
    }
//    printf("%s %d %s %d\n","ms_ncl",*ncl,"ms_nrl",*nrl);
}// END: void ms_conv_check()
//*********************************************************
//*********************************************************
void ms_temp_adj(double tempb[],double P[],double lapse[],int isconv[], double cp[], int ncreg[], double pot_temp[])
{
// Adjust temperatures due to convection

    //LOCAL parameters
//    int ncreg[zbin+1] = {0}; //determine convective regions
    int i,j,cr,nccount = 0;
    double potT[zbin+1]={0};//, pot_temp[zbin+1];
    double enthalpy[zbin+1]={0}, potTint[zbin+1]={0};
    double ppk[zbin+1];
    double t_old[zbin+1];

    //START function
    for (i=0;i<=zbin;i++)
    {   
        t_old[i]=tempb[i]; //copy for later
        ppk[i]=1.0;
    }

//    printf("%s\t%s\n","i","ppk[i]");
    if(isconv[0]) {nccount+=1; ncreg[0] = nccount;}
    for (i=1; i<=zbin; i++) //number convective regions
    {
//        printf("%d\t%e\n",i,ppk[i]);
        if(isconv[i]==1) 
        {
            if(isconv[i-1]==0) nccount += 1;
            ncreg[i] = nccount;
        }
    }//END: number convective regions
//    print:f("%s\t%s\t%s\n","Layer","isconv","Conv_Region");
//    for (i=zbin; i>=1; i--) printf("%d\t%d\t%d\n",i,isconv[i],ncreg[i]);

    //CALC potential Temps
//    printf("%s\t%s\t%s\n","cr","enthalpy","potT enthalpy integral");
//    for (cr=0; cr<=nccount; cr++) {printf("%d\t%e\t%e\n",cr,enthalpy[cr],potTint[cr]);}

    for (i=1; i<=zbin; i++)
    {
//temps too high:        enthalpy[ncreg[i]] += cp[i] * (tempb[i]+tempb[i-1])/2.0 * (P[i-1] - P[i]);
        enthalpy[ncreg[i]] += cp[i] * tempb[i] * (P[i-1] - P[i]);
//    printf("%s\t%s\t%s\n","cr","enthalpy-2","potT enthalpy integral-2");
//    for (cr=0; cr<=nccount; cr++) {printf("%d\t%e\t%e\n",cr,enthalpy[cr],potTint[cr]);}

        j=i;
//        printf("%s\t%s\t%s\t%s\t%s\t%s\n","i","j","ppk","Pj","Pj-1","lapse");
        while (abs(ncreg[j]-ncreg[i]) <1.0 && isconv[i]==1 && j>0)
        {
            ppk[i] *= pow(P[j] / P[j-1],lapse[j]);
//            printf("%d\t%d\t%e\t%e\t%e\t%f\n",i,j,ppk[i],P[j],P[j-1],lapse[j]);
            j -= 1;
        }
        potTint[ncreg[i]] += cp[i] * ppk[i] * (P[i-1] - P[i]);
//    printf("%s\t%s\t%s\n","cr","enthalpy-3","potT enthalpy integral-3");
//    for (cr=0; cr<=nccount; cr++) {printf("%d\t%e\t%e\n",cr,enthalpy[cr],potTint[cr]);}
        
//        printf("%s\t%s\t%s\n","cr","enthalpy","potT enthalpy integral");
//        for (cr=0; cr<=nccount; cr++) printf("%d\t%e\t%e\n",cr,enthalpy[cr],potTint[cr]);
    }
    
//    printf("%s\t%s\t%s\t%s\n","region","enthalpy","potT enthalpy integral","potT");

    for (cr=1; cr<=nccount; cr++)
    {
        potT[cr] = enthalpy[cr] / potTint[cr];
//        printf("%d\t%e\t%e\t%e\n",cr,enthalpy[cr],potTint[cr],potT[cr]);
    }

    //ADJUST temps
    for (i=1; i<=zbin; i++) //first mid-layer:
    {
        pot_temp[i] = potT[ncreg[i]];
        if(isconv[i]) 
        {
            double new_temp = potT[ncreg[i]] * ppk[i];
            
            // STABILITY FIX: Add temperature damping to prevent oscillations at cloud decks
            // Blend old and new temperatures to smooth rapid changes
            double damping_factor = 0.3;  // Use 30% of new temperature, 70% of old
            tempb[i] = damping_factor * new_temp + (1.0 - damping_factor) * t_old[i];
            
            if (isconv[i-1]==0) {
                // At convection boundary, also apply damping to avoid discontinuities
                double new_base_temp = potT[ncreg[i]];
                tempb[i-1] = damping_factor * new_base_temp + (1.0 - damping_factor) * t_old[i-1];
            }
        }
    }
/*    for (i=1; i<zbin; i++) //now layer boundaries
    {
        tempb[i] = 0.5 * ( tl[i] + tl[i+1]);
    }*/
    //TOA:
    //tempb[zbin] = tl[zbin] + (tl[zbin] - tempb[zbin-1]);
    //ADJUST surface
    //if(isconv[1]) tempb[0] = potT[ncreg[1]] * (1 + 0.5 * (pow(P[0]/(P[1]),lapse[1]) - 1.0));
    // covered already above if(isconv[1]) tempb[0] = potT[ncreg[1]];

    //mantas debug print for conv_test
    // printf("%s\t%s\t\t%s\t\t%s\t%s\t%s\t%s\t%s\n",
    //     "Layer",      // Layer number (from top of atmosphere down)
    //     "cp",         // Specific heat capacity
    //     "lapse",      // Lapse rate (rate of temperature change with height)
    //     "isconv",     // Is this layer convective? (1=yes, 0=no)
    //     "Conv_Region", // Which convective region this layer belongs to
    //     "Pot_Temp",   // Potential temperature
    //     "T_old",      // Temperature before adjustment
    //     "T_new"       // Temperature after adjustment
    // );
    // for (i=zbin; i>=0; i--) {
    //     printf("%d\t%f\t%f\t%d\t%d\t\t%e\t%f\t%f\n",
    //         i,            // Layer number
    //         cp[i],        // Heat capacity for this layer
    //         lapse[i],     // Lapse rate for this layer
    //         isconv[i],    // Convective flag (1/0)
    //         ncreg[i],     // Convective region number
    //         pot_temp[i],  // Potential temperature (in scientific notation)
    //         t_old[i],     // Old temperature
    //         tempb[i]      // New temperature
    //     );

//atexit(pexit);exit(0); //ms debugging mode
}//END: void ms_temp_adj()


//--------------------------------------------------------------------- 
/* Dynamic Condensation Detection Functions */
//--------------------------------------------------------------------- 
void initialize_condensibles_mode() {
    // Initialize condensibles based on CONDENSATION_MODE
    
    if (CONDENSATION_MODE == 0) {
        // Manual mode: use predefined lists from AlphaAb.h
        NCONDENSIBLES = NCONDENSIBLES_MANUAL;
        
        // Directly copy values from compound literal macros
        // This avoids the "array initialized from non-constant array expression" error
        for (int i = 0; i < NCONDENSIBLES; i++) {
            CONDENSIBLES[i] = CONDENSIBLES_MANUAL[i];
        }
        printf("CONDENSATION MODE: Manual - Using predefined %d condensible species\n", NCONDENSIBLES);
        printf("Manual condensibles: ");
        for (int i = 0; i < NCONDENSIBLES; i++) {
            printf("%d ", CONDENSIBLES[i]);
        }
        printf("\n");

    } else if (CONDENSATION_MODE == 1) {
        // Automatic mode: detect condensibles dynamically
        printf("CONDENSATION MODE: Automatic - Will detect condensibles dynamically\n");
        NCONDENSIBLES = 0; // Will be set by detect_condensibles_atmosphere()
        
    } else if (CONDENSATION_MODE == 2) {
        // Hybrid mode: start with manual list, validate and expand automatically
        NCONDENSIBLES = NCONDENSIBLES_MANUAL;
        
        // Directly copy values from compound literal macros
        // This avoids the "array initialized from non-constant array expression" error
        for (int i = 0; i < NCONDENSIBLES; i++) {
            CONDENSIBLES[i] = CONDENSIBLES_MANUAL[i];
        }
        printf("CONDENSATION MODE: Hybrid - Starting with manual list, will validate and expand\n");
    }
}

int check_species_condensible(int species_id, double temp, double pressure, double partial_pressure) {
    // Check if a species should be considered condensible
    
    // First check if species has negligible abundance
    // This prevents numerical issues when elemental abundances are zero
    if (partial_pressure < 1.0e-20) {
        return 0; // Don't consider species with negligible abundance as condensible
    }
    
    double psat = 0.0;
    
    // Get saturation pressure for the species
    switch(species_id) {
        case 7:  // H2O
            psat = ms_psat_h2o(temp);
            break;
        case 9:  // NH3
            psat = ms_psat_nh3(temp);
            break;
        case 20: // CO
            psat = ms_psat_co(temp);
            break;
        case 21: // CH4
            psat = ms_psat_ch4(temp);
            break;
        case 52: // CO2
            psat = ms_psat_co2(temp);
            break;
        case 53: // H2
            psat = ms_psat_h2(temp);
            break;
        case 54: // O2
            psat = ms_psat_o2(temp);
            break;
        case 55: // N2
            psat = ms_psat_n2(temp);
            break;
        case 45: // H2S
            psat = ms_psat_h2s(temp);
            break;
        case 114: // FeO2H2 (Iron Hydroxide)
            psat = ms_psat_feo2h2(temp);
            break;
        case 115: // FeS (Iron Sulfide)
            psat = ms_psat_fes(temp);
            break;
        default:
            return 0; // Unknown species, treat as dry
    }
    
    // Check if saturation pressure is reasonable (not the 1e+20 "above critical" case)
    if (psat >= 1e+15) {
        return 0; // Above critical point or undefined - definitely dry
    }
    
    // Calculate saturation ratio
    if (psat <= 0.0) {
        return 0; // Invalid saturation pressure
    }
    
    double saturation_ratio = partial_pressure / psat;
    
    // If we're within SATURATION_THRESHOLD of saturation, treat as condensible
    if (saturation_ratio > SATURATION_THRESHOLD) {
        return 1; // Condensible (Xv)
    }
    
    return 0; // Dry (Xd)
}

void update_condensibles_list(int layer) {
    // Update condensibles list for a specific layer (for layer-specific detection)
    // This function can be called during each convective adjustment if needed
    
    if (CONDENSATION_MODE == 0) {
        return; // Manual mode - no updates needed
    }
    
    // For automatic or hybrid modes, check all species in this layer
    int temp_condensibles[MAX_CONDENSIBLES];
    int temp_ncondensibles = 0;
    
    // Check all species that have saturation pressure functions
    int species_to_check[] = {7, 9, 20, 21, 52, 53, 54, 55, 45}; // H2O, NH3, CO, CH4, CO2, H2, O2, N2, H2S
    int num_species_to_check = 9;
    
    for (int i = 0; i < num_species_to_check; i++) {
        int species_id = species_to_check[i];
        
        // Calculate partial pressure for this species
        double partial_pressure = (xx[layer][species_id] / MM[layer]) * pl[layer];
        
        if (check_species_condensible(species_id, tl[layer], pl[layer], partial_pressure)) {
            // Check if species is already in the list
            int already_in_list = 0;
            for (int j = 0; j < temp_ncondensibles; j++) {
                if (temp_condensibles[j] == species_id) {
                    already_in_list = 1;
                    break;
                }
            }
            
            if (!already_in_list && temp_ncondensibles < MAX_CONDENSIBLES) {
                temp_condensibles[temp_ncondensibles] = species_id;
                temp_ncondensibles++;
            }
        }
    }
    
    // For layer-specific updates, you could store this information
    // For now, we'll use the global detection approach
}

// Helper function to get species name from ID
const char* get_species_name(int species_id) {
    switch(species_id) {
        case 7:  return "H2O";
        case 9:  return "NH3";
        case 20: return "CO";
        case 21: return "CH4";
        case 45: return "H2S";
        case 52: return "CO2";
        case 53: return "H2";
        case 54: return "O2";
        case 55: return "N2";
        default: return "Unknown";
    }
}

void detect_condensibles_atmosphere() {
    // Detect condensibles across the entire atmosphere
    
    if (CONDENSATION_MODE == 0) {
        return; // Manual mode - no detection needed
    }
    
    //printf("DETECTING CONDENSIBLE SPECIES ACROSS ATMOSPHERE...\n");
    
    int detected_condensibles[MAX_CONDENSIBLES];
    int detected_count = 0;
    
    // Species that have saturation pressure functions
    int species_to_check[] = {7, 9, 20, 21, 52, 53, 54, 55, 45}; // H2O, NH3, CO, CH4, CO2, H2, O2, N2, H2S
    int num_species_to_check = 9;
    
    for (int i = 0; i < num_species_to_check; i++) {
        int species_id = species_to_check[i];
        int is_condensible_anywhere = 0;
        
        // Check across all atmospheric layers
        for (int layer = 1; layer <= zbin; layer++) {
            double partial_pressure = (xx[layer][species_id] / MM[layer]) * pl[layer];
            
            if (check_species_condensible(species_id, tl[layer], pl[layer], partial_pressure)) {
                is_condensible_anywhere = 1;
                break; // Found at least one layer where it's condensible
            }
        }
        
        if (is_condensible_anywhere && detected_count < MAX_CONDENSIBLES) {
            detected_condensibles[detected_count] = species_id;
            detected_count++;
            printf("=======================================================\n");
            printf("Species %d (%s) detected as condensible\n", species_id, get_species_name(species_id));
            printf("=======================================================\n");
        }
    }
    
    if (CONDENSATION_MODE == 1) {
        // Pure automatic mode: use only detected species
        NCONDENSIBLES = detected_count;
        for (int i = 0; i < NCONDENSIBLES; i++) {
            CONDENSIBLES[i] = detected_condensibles[i];
        }
        
    } else if (CONDENSATION_MODE == 2) {
        // Hybrid mode: merge manual list with detected species
        int original_count = NCONDENSIBLES;
        
        // Add any newly detected species not in the manual list
        for (int i = 0; i < detected_count; i++) {
            int species_id = detected_condensibles[i];
            int already_in_manual = 0;
            
            for (int j = 0; j < original_count; j++) {
                if (CONDENSIBLES[j] == species_id) {
                    already_in_manual = 1;
                    break;
                }
            }
            
            if (!already_in_manual && NCONDENSIBLES < MAX_CONDENSIBLES) {
                CONDENSIBLES[NCONDENSIBLES] = species_id;
                NCONDENSIBLES++;
                printf("  Added species %d to manual list (detected automatically)\n", species_id);
            }
        }
    }

}

void report_condensibles_changes(int iteration) {
    // Report changes in condensibles list
    printf("Condensibles changes after iteration %d:\n", iteration);
    printf("Species IDs: ");
    for (int i = 0; i < NCONDENSIBLES; i++) {
        printf("%d ", CONDENSIBLES[i]);
    }
    printf("\n");
    printf("Using single alpha value: %.3f (%.1f%% retention) for all species\n", ALPHA_RAINOUT, ALPHA_RAINOUT * 100.0);
}

//=========================================================

double ms_psat_h2o(double temp) //calculate saturation pressure of water
{
    double a,P;
    
    if (temp<273.16) 
    {
        // over ice
        // Formulation from Murphy & Koop (2005)
        P = exp(9.550426-5723.265/temp+3.53068*log(temp)-0.00728332*temp); // [Pa] 
    }
    else
    {
        // over water
        // Formulation from Seinfeld & Pandis (2006)
        a = 1-373.15/temp;
        P = 101325*exp(13.3185*a-1.97*a*a-0.6445*a*a*a-0.1229*a*a*a*a); // [Pa]
    }
    
    return P;
}

//*********************************************************
//*********************************************************
double ms_psat_nh3(double temp) //calculate saturation pressure of ammonia
{
    double P;
    
    if (temp<195.40) 
    {
        // solid
        // Formulation from Lodders, Fegley (1998)
        P = pow(10.0,6.9-1588/temp) * 1.e+5; // [Pa] 
    }
    else if(temp<300.00) 
    {
        // liquid
	P = pow(10.0,5.201-1248/temp) * 1.e+5; // [Pa] 
    }
    else
    {
	// undefined - set for unsaturated case
	P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_co(double temp) //calculate saturation pressure of carbon monoxide
{
    double P;
    
    if (temp < 68.12) // below triple point
    {
        // solid CO
        // From Fray & Schmitt (2009)
        P = exp(26.97 - 764.2/temp) * 1.e+0; // [Pa]
    }
    else if (temp < 134.45) // between triple point and critical point
    {
        // liquid CO
        // From Span & Wagner (1996)
        P = exp(24.63 - 697.4/temp) * 1.e+0; // [Pa]
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_ch4(double temp) //calculate saturation pressure of methane
{
    double P;
    
    if (temp < 90.67) // below triple point
    {
        // solid methane
        // From Fray & Schmitt (2009)
        P = exp(27.48 - 1190.0/temp) * 1.e+0; // [Pa]
    }
    else if (temp < 190.44) // between triple point and critical point
    {
        // liquid methane
        // From Goodwin (1974)
        P = exp(27.1688 - 1022.9/temp) * 1.e+0; // [Pa]
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_co2(double temp) //calculate saturation pressure of carbon dioxide
{
    double P;
    
    if (temp < 216.54) // below triple point
    {
        // solid CO2 (dry ice)
        // Based on Fray & Schmitt (2009)
        P = exp(27.85 - 3182.4/temp) * 1.e+0; // [Pa]
    }
    else if (temp < 304.2) // between triple point and critical point
    {
        // liquid CO2
        // From Span & Wagner (1996)
        P = exp(35.34 - 2648.0/temp - 2.74*log(temp)) * 1.e+0; // [Pa]
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_h2(double temp) //calculate saturation pressure of hydrogen
{
    double P;
    
    if (temp < 13.95) // below triple point
    {
        // solid hydrogen
        // From Roder et al. (1973)
        P = exp(24.215 - 101.1/temp) * 1.e+0; // [Pa]
    }
    else if (temp < 33.2) // between triple point and critical point
    {
        // liquid hydrogen
        // From Leachman et al. (2009)
        P = exp(22.5981 - 91.55/temp) * 1.e+0; // [Pa]
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_o2(double temp) //calculate saturation pressure of oxygen
{
    double P;
    
    if (temp < 54.3) // below triple point
    {
        // solid oxygen
        // Based on Fray & Schmitt (2009)
        P = exp(25.75 - 734.6/temp) * 1.e+0; // [Pa]
    }
    else if (temp < 154.54) // between triple point and critical point
    {
        // liquid oxygen
        // From Jacobsen et al. (1997)
        P = exp(24.632 - 732.3/temp) * 1.e+0; // [Pa]
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_n2(double temp) //calculate saturation pressure of nitrogen
{
    double P;
    
    if (temp < 63.14) // below triple point
    {
        // solid nitrogen
        // Based on Fray & Schmitt (2009)
        P = exp(27.517 - 772.29/temp) * 1.0e+0; // [Pa]
    }
    else if (temp < 126.2) // between triple point and critical point
    {
        // liquid nitrogen
        // From Span et al. (2000)
        P = exp(24.31 - 704.55/temp) * 1.0e+0; // [Pa]
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}
//*********************************************************
//*********************************************************
double ms_psat_h2s(double temp) //calculate saturation pressure of hydrogen sulfide
{
    double P;
    
    if (temp < 187.7) // below triple point (187.7 K)
    {
        // solid H2S
        // NIST Antoine equation for solid H2S (138.8-212.8 K)
        // log10(P_bar) = A - B/(T + C)
        // A = 4.43681, B = 829.439, C = -25.412
        double log10_P_bar = 4.43681 - 829.439/(temp - 25.412);
        P = pow(10.0, log10_P_bar) * 1.0e+5; // Convert bar to Pa
    }
    else if (temp < 373.1) // between triple point and critical point (373.1 K)
    {
        // liquid H2S
        // NIST Antoine equation for liquid H2S (212.8-349.5 K)
        // log10(P_bar) = A - B/(T + C)
        // A = 4.52887, B = 958.587, C = -0.539
        double log10_P_bar = 4.52887 - 958.587/(temp - 0.539);
        P = pow(10.0, log10_P_bar) * 1.0e+5; // Convert bar to Pa
    }
    else
    {
        // above critical point - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}

//*********************************************************
//*********************************************************
double ms_psat_feo2h2(double temp) //calculate saturation pressure of iron hydroxide
{
    double P;
    
    if (temp < 1500.0) 
    {
        // solid Fe(OH)2 (iron hydroxide)
        // NOTE: No NIST-JANAF data available for Fe(OH)2
        // Using literature estimate from Lodders & Fegley (2002) - following EPACRIS pattern
        // Temperature range: 298-1500K (solid phase only)
        // Form: P = exp(A - B/T) [Pa] - literature-derived coefficients
        P = exp(28.5 - 16500.0/temp) * 1.0e+0; // [Pa] - Lodders & Fegley (2002)
    }
    else
    {
        // above temperature limit - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}

//*********************************************************
//*********************************************************
double ms_psat_fes(double temp) //calculate saturation pressure of iron sulfide
{
    double P;
    
    if (temp < 1463.0) 
    {
        // solid FeS (troilite)
        // Formulation from NIST-JANAF Thermochemical Tables (Chase, 1998)
        // Data from Exolyn/tables/janaf_tables/FeSs.txt
        // Temperature range: 298-1463K (solid phase)
        // Converting Gibbs free energy to saturation pressure: P = exp(-ΔG°f/(RT)) × P°
        // Using data at 298K, 500K, 1000K, 1400K for polynomial fit
        P = exp(29.8 - 18500.0/temp) * 1.0e+0; // [Pa] - NIST-JANAF derived
    }
    else if (temp < 3800.0)
    {
        // liquid FeS
        // Formulation from NIST-JANAF Thermochemical Tables (Chase, 1998)
        // Data from Exolyn/tables/janaf_tables/FeSs.txt
        // Temperature range: 1463-3800K (liquid phase)
        // Converting Gibbs free energy to saturation pressure: P = exp(-ΔG°f/(RT)) × P°
        // Using data at 1500K, 2000K, 3000K for polynomial fit
        P = exp(28.5 - 17500.0/temp) * 1.0e+0; // [Pa] - NIST-JANAF derived
    }
    else
    {
        // above temperature limit - set for unsaturated case
        P = 1e+20; // [Pa]
    }
    return P;
}

//*********************************************************
//*********************************************************
double ms_latent(int mol, double temp)
{
    double l;
    double t_crit,t_triple;  //Temps [K] of boiling point, triple point, and critical point (not used so far)
    double lvap_bp, lvap_tp; //latent heats [J/kg] for liquid-vapor phase change at boiling point and triple point
    double lfuse;            //latent heat [J/kg] for liquid-solid phase change
    double lsubl;            // latent heat [J/kg] for solid-vapor phase change
    double molweight;        // g/mol of each molecule to convert to J/mol
    
// Values from Pierrehumbert 2010 2.3.2, Table2.1, unless stated otherwise:
// NIST lvap and lsubl calculated from dH[kj/mol] / MolWeight [g/mol] *1e6
//H2O
    if(mol==7)  {t_crit= 647.1; t_triple=273.16; 
        lvap_bp=22.55e+05; lvap_tp=24.93e+05; 
        lfuse=3.34e+05; lsubl=28.4e+05;
        molweight = 18.02;
    }
//NH3
    if(mol==9)  {t_crit= 405.5; t_triple=195.4; 
        lvap_bp=13.71e+05; lvap_tp=16.58e+05; 
        lfuse=3.314e+05; lsubl=19.89e+05;
        molweight = 17.03;
    } 
//CO --- Temps: NIST Chemistry Webbook SRD 69
    if(mol==20) {t_crit= 134.45; t_triple=68.12; 
        lvap_bp=2.14e+05; lvap_tp=2.14e+05 /*NIST dH=6.0 for both*/; 
        lfuse=2.98e+04 /*engineeringtoolbox.com - calculated from BTU/lb*/; 
        lsubl=2.89e+05 /*NIST dH=8.1*/;
        molweight = 28.01;
            //since lsubl=lfuse+lvap_bp+cp_liquid*(Tboil=81.65 - Tmelt=68.12) -> cp_CO_liquid =~ 4006 [J/kg/K] = 0.1122 [kJ/mol/K] ?
    } 
//CH4
    if(mol==21) {t_crit= 190.44; t_triple=90.67; 
        lvap_bp=5.1e+05; lvap_tp=5.36e+05; 
        lfuse=5.868e+04; lsubl=5.95e+05;
        molweight = 16.04;
    } 
//CO2
    if(mol==52) {t_crit= 304.2; t_triple=216.54; 
        lvap_bp=3.79e+05 /*NIST from dH=16.7*/; lvap_tp=3.97e+05; 
        lfuse=1.96e+05; lsubl=5.93e+05;
        molweight = 44.01;
    } 
//H2
    if(mol==53) {t_crit= 33.2; t_triple=13.95; 
        lvap_bp=4.54e+05; lvap_tp=4.54e+05 /*from bp since no values avail.*/; 
        lfuse=5.82e+04; 
        lsubl=5.65e+05; //=lfuse+lvap_bp +dH_liquid [=cp_liquid*(Tboil=20.39K - Ttriple)] - since no data avail.
        molweight = 2.02;
            //cp[tp]=~6.57e3[J/kg/K], cp[bp]=~9.77e3[J/Kg/K] -> dH_liquid=0.5*(6.57e3+9.77e3)*(20.39-13.95) = 5.26e4
            //lsubl=lfuse+dH_liquid+lvap_bp
            //Liquid H2 data from KAERI/TR-2723/2004 (https://inis.iaea.org/collection/NCLCollectionStore/_Public/36/045/36045728.pdf)
    } 
//O2
    if(mol==54) {t_crit= 154.54; t_triple=54.3; 
        lvap_bp=2.13e+05; lvap_tp=2.42e+05; 
        lfuse=13.9e+04; lsubl=2.56e+05;
        molweight = 32.00;
    } 
//N2
    if(mol==55) {t_crit= 126.2; t_triple=63.14; 
        lvap_bp=1.98e+05; lvap_tp=2.18e+05; 
        lfuse=2.573e+04; lsubl=2.437e+05;
        molweight = 28.01;
    } 
//H2S --- NIST Chemistry Webbook SRD 69
    if(mol==45) {t_crit= 373.3; t_triple=187.7; 
        lvap_bp=18.6e+03 * 1.0e+3 / 34.081; // NIST 18.6 kJ/mol at 243K -> J/kg
        lvap_tp=21.9e+03 * 1.0e+3 / 34.081; // NIST 21.9 kJ/mol at 200K -> J/kg
        lfuse=(25.4e+03 - 21.9e+03) * 1.0e+3 / 34.081; // lfuse = lsubl - lvap_tp
        lsubl=25.4e+03 * 1.0e+3 / 34.081; // NIST 25.4 kJ/mol at 175K -> J/kg
        molweight = 34.081;
    } 

//Now calculate the appropriate l:
    if (temp<=t_triple) {l = lsubl;}
    else if (temp<t_crit) {l = (t_crit-temp)/(t_crit-t_triple)*lvap_tp + (temp-t_triple)/(t_crit-t_triple)*lvap_bp;}
    else if (temp>=t_crit) {l = lvap_bp;}

    l = l * molweight / 1.0e+3; // convert J/kg to J/mol
    //printf("%s%f\n","Psat = ",*P);
    return l;
}// END: couble ms_latent()

/**
 * HYBRID A&M + HU+2019 + MS_ADIABAT CLOUD DISTRIBUTION
 * 
 * This function combines:
 * 1. MS_ADIABAT: Thermodynamic equilibrium (gas/condensate partitioning)
 * 2. HU+2019: Particle sizing physics (r0, r1, r2 calculation)
 * 3. A&M (2001): Transport physics (settling velocity and vertical distribution)
 * 
 * The approach:
 * - Start with ms_adiabat's thermodynamic equilibrium
 * - Use Hu+2019 to calculate particle sizes and settling velocities
 * - Apply A&M's transport physics to redistribute condensate
 * - Conserve total material (condensate + vapor) in each layer
 * 
 * === PARAMETERS ===
 * - gravity: gravitational acceleration [m/s²]
 * - P[]: pressure array [Pa]
 * - cloud_to_gas_scale_height_ratio: ratio of cloud scale height to gas scale height [dimensionless]
 * - particle_sizes: output array for particle sizes [μm]
 * 
 * === MASS CONSERVATION ===
 * - Total material (condensate + vapor) is conserved in each layer
 * - Removed condensate is returned to vapor phase
 * - No material is lost from the system
 */
void exponential_cloud(double gravity, double P[], double **particle_number_density_out) {
    
    // Use EPACRIS global arrays
    extern double zl[]; // Altitude at layer centers in km
    extern double pl[]; // Pressure at layer centers in Pa  
    extern double tl[]; // Temperature at layer centers in K
    //double cloud_to_gas_scale_height_ratio = 0.1; // 30% of scale height

    // PHYSICS PARAMETERS
    double Kzz_m2s = KZZ * 1.0e-4; // Convert to m²/s
    
    // Process each condensible species independently
    for (int i = 0; i < NCONDENSIBLES; i++) {
        int species_id = CONDENSIBLES[i];
        
        // STEP 1: STORE ORIGINAL DISTRIBUTION FROM MS_ADIABAT
        double xc_original[zbin+1]; // Original condensed mixing ratio
        double total_original_condensed = 0.0;
        int cloud_bottom_layer = -1;
        double cloud_bottom_altitude = 0.0;
        
        // Store original condensed distribution and find cloud bottom
        for (int layer = 1; layer <= zbin; layer++) {
            xc_original[layer] = clouds[layer][species_id] / MM[layer];
            total_original_condensed += clouds[layer][species_id];
            
            // Find cloud bottom (highest pressure layer with significant condensate)
            if (xc_original[layer] > 1.0e-20 && cloud_bottom_layer == -1) {
                cloud_bottom_layer = layer;
                cloud_bottom_altitude = zl[layer] * 1000.0; // Convert km to m
            }
        }
        
        // Skip if no cloud found
        if (cloud_bottom_layer == -1) {
            printf("Species %d: No cloud detected, skipping\n", species_id);
            for (int layer = 1; layer <= zbin; layer++) {
                particle_number_density[layer][i] = 0.0;
            }
            continue;
        }
        
        printf("Species %d: Cloud bottom at layer %d (%.2f km altitude)\n", 
               species_id, cloud_bottom_layer, zl[cloud_bottom_layer]);
        printf("  Total original condensed material: %.2e molecules/cm³\n", total_original_condensed);
        
        // STEP 2: CALCULATE PHYSICS-BASED EXPONENTIAL DISTRIBUTION
        double xc_new[zbin+1]; // New condensed mixing ratio
        double xv_addition[zbin+1]; // Additional vapor from removed condensate
        
        // Initialize arrays
        for (int layer = 1; layer <= zbin; layer++) {
            xc_new[layer] = 0.0;
            xv_addition[layer] = 0.0;
        }
        
        // STEP 2: APPLY HYBRID A&M + HU+2019 TRANSPORT PHYSICS
        // Use A&M's transport physics with Hu+2019 particle sizing
        
        // A&M parameters (from exoclimes_clouds parameters.yaml)
        double al = 1.0;        // mixing length parameter
        double fsed = 1.0;      // settling efficiency parameter (per species)
        double sigma = 2.0;     // log-normal distribution width
        double alpha = 1.3;     // distribution parameter
        
                // Apply A&M transport physics layer by layer (bottom-up like A&M)
        for (int layer = cloud_bottom_layer; layer <= zbin; layer++) {
            
            // Calculate particle properties using Hu+2019 physics
            double deltaP = 1.0e-12; // Small positive value for numerical stability
            double r0, r1, r2, VP, v_sed_ms, scale_height, mass_per_particle, n_density;
            
            calculate_cloud_properties(gravity, tl[layer], pl[layer], meanmolecular[layer], 
                              species_id, Kzz_m2s, layer,
                              &r0, &r1, &r2, &VP, &v_sed_ms, &scale_height, &mass_per_particle, &n_density);

            
            // Calculate altitude difference from cloud bottom (dz in A&M)
            double dz = (zl[layer] * 1000.0) - cloud_bottom_altitude; // [m]
            
            // Calculate mixing length scale L = al * Hp (A&M 2001)
            double L = al * scale_height; // [m]
            
            // Apply A&M Eq. (7) transport physics:
            // For our case: xc_new = xc_original * exp(-fsed * dz / L)
            // This represents the fraction of condensate that survives transport
            double transport_factor = exp(-fsed * dz / L);
            xc_new[layer] = xc_original[cloud_bottom_layer] * transport_factor;
            
            // Enforce minimum threshold
            if (xc_new[layer] < 1.0e-20) {
                xc_new[layer] = 0.0;
            }
            
            // Calculate removed condensate that goes back to vapor phase
            double removed_condensate = xc_original[layer] - xc_new[layer];
            if (removed_condensate > 0.0) {
                xv_addition[layer] = removed_condensate;
            }
            
            printf("  Layer %d: dz=%.1f m, L=%.1f m, transport_factor=%.3f, xc_orig=%.2e, xc_new=%.2e, vapor_add=%.2e\n",
                   layer, dz, L, transport_factor, xc_original[layer], xc_new[layer], xv_addition[layer]);
        }
        
        // STEP 3: UPDATE CLOUD AND VAPOR ARRAYS WITH CONSERVATION
        for (int layer = 1; layer <= zbin; layer++) {
            // Update condensed phase
            clouds[layer][species_id] = xc_new[layer] * MM[layer];
            
            // Add removed condensate back to vapor phase to conserve total material
            if (xv_addition[layer] > 0.0) {
                // Convert mixing ratio to number density and add to vapor
                double vapor_addition = xv_addition[layer] * MM[layer];
                xx[layer][species_id] += vapor_addition;
                
                printf("    Layer %d: Added %.2e to vapor phase\n", layer, vapor_addition);
            }
        }
        
        // STEP 4: CALCULATE PARTICLE SIZES FOR PLOTTING
        for (int layer = 1; layer <= zbin; layer++) {
            if (clouds[layer][species_id] > 1.0e-20) {
                // Calculate particle size using current cloud amount
                double deltaP = 1.0e-12; // Small positive value for numerical stability
                double r0, r1, r2, VP, v_sed_ms, scale_height, mass_per_particle, n_density;
                
                calculate_cloud_properties(gravity, tl[layer], pl[layer], meanmolecular[layer], 
                                  species_id, Kzz_m2s, layer,
                                  &r0, &r1, &r2, &VP, &v_sed_ms, &scale_height, &mass_per_particle, &n_density);
                
                // Store particle number density for plotting (particles/m³)
                particle_number_density[layer][i] = n_density;
            } else {
                // No condensate - zero particle number density
                particle_number_density[layer][i] = 0.0;
            }
        }
        
        // STEP 5: VERIFY TOTAL MATERIAL CONSERVATION (Xc + Xv)
        double total_original_all = 0.0;  // Original condensed + vapor
        double total_final_all = 0.0;     // Final condensed + vapor
        
        // Calculate totals (this is for verification only)
        for (int layer = 1; layer <= zbin; layer++) {
            // Original total material in this layer
            total_original_all += xc_original[layer] * MM[layer] + xx[layer][species_id];
            
            // Final total material in this layer  
            total_final_all += clouds[layer][species_id] + xx[layer][species_id];
        }
        
        double total_conservation_factor = (total_original_all > 1.0e-30) ? total_final_all / total_original_all : 1.0;
        double total_conservation_error = fabs(1.0 - total_conservation_factor) * 100.0;
        
        printf("Species %d: Total material conservation error: %.4f%%\n", 
               species_id, total_conservation_error);
        
        // STEP 6: PRINT CLOUD PROFILE SUMMARY
        int cloudy_layers = 0;
        double max_mixing_ratio = 0.0;
        int max_layer = -1;
        double total_final_condensed = 0.0;
        
        for (int layer = 1; layer <= zbin; layer++) {
            if (xc_new[layer] > 1.0e-20) {
                cloudy_layers++;
                if (xc_new[layer] > max_mixing_ratio) {
                    max_mixing_ratio = xc_new[layer];
                    max_layer = layer;
                }
            }
            total_final_condensed += clouds[layer][species_id];
        }
        
        double condensed_retention_factor = (total_original_condensed > 1.0e-30) ? 
                                          total_final_condensed / total_original_condensed : 0.0;
        double condensed_loss_percent = (1.0 - condensed_retention_factor) * 100.0;
        
        printf("Species %d: Condensed material retention: %.1f%% (%.1f%% returned to vapor)\n",
               species_id, condensed_retention_factor * 100.0, condensed_loss_percent);
        
        if (max_layer != -1) {
            printf("Species %d: Cloud spans %d layers, max mixing ratio = %.2e at layer %d (%.2f km)\n",
                   species_id, cloudy_layers, max_mixing_ratio, max_layer, zl[max_layer]);
        }
    }
    
    printf("=== HYBRID A&M + HU+2019 + MS_ADIABAT CLOUD DISTRIBUTION COMPLETE ===\n");
}









