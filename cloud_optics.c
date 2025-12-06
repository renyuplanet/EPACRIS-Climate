/**
 * cloud_optics.c
 * 
 * Functions to read and store cloud optical property lookup tables from MIE tables
 * Cloud species is configured via CLOUD_SPECIES in config.h (species ID, e.g., 7 for H2O)
**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include "cloud_optics.h"
#include "global_temp.h"  // For zbin, NSP, NLAMBDA, NCONDENSIBLES, CONDENSIBLES
#include "conv_cond_funcs.h"  // For get_species_name()
#include "config.h"  // For CLOUD_MIE_DIRECTORY and CLOUD_SPECIES
#include "ms_functions.h"

// Include config file for debug switches (if not already included)
#ifndef CLOUD_DEBUG_MIE
#define CLOUD_DEBUG_MIE 0  // Default to disabled if not defined
#endif

// Global optical property tables - array for multiple cloud species
CloudOpticalTableMie cloud_mie_optics[MAX_CLOUD_SPECIES];
int cloud_species_ids[MAX_CLOUD_SPECIES];
int n_cloud_species_loaded = 0;

//--------------------------------------------------------------------- 
//  For LX-Mie format only
//  Extract particle radius from filename like "r0.010000.dat"
//  Returns radius in micrometers, or -1.0 on error
//--------------------------------------------------------------------- 
static double extract_radius_from_filename(const char *filename) {
    // Skip "r" prefix
    if (filename[0] != 'r') {
        return -1.0;
    }
    
    // Find ".dat" suffix
    char *dot_dat = strstr(filename, ".dat");
    if (dot_dat == NULL) {
        return -1.0;
    }
    
    // Extract number between "r" and ".dat"
    char num_str[64];
    int len = dot_dat - (filename + 1);
    if (len <= 0 || len >= 64) {
        return -1.0;
    }
    
    strncpy(num_str, filename + 1, len);
    num_str[len] = '\0';
    
    return atof(num_str);
}

/**
 * Compare function for sorting filenames by particle radius
 */
static int compare_particle_files(const void *a, const void *b) {
    const char *file_a = *(const char **)a;
    const char *file_b = *(const char **)b;
    
    double radius_a = extract_radius_from_filename(file_a);
    double radius_b = extract_radius_from_filename(file_b);
    
    if (radius_a < radius_b) return -1;
    if (radius_a > radius_b) return 1;
    return 0;
}

/**
 * Read a single LX-Mie output file
 * 
 * File format:
 * # Header line (skip)
 * wavelength(μm)  size_param  extinction(cm²)  scattering(cm²)  absorption(cm²)  albedo  asymmetry_g
 */
static int read_mie_file(const char *filepath, double *wavelengths, 
                         double *extinction, double *albedo, double *asymmetry,
                         int *n_wavelengths) {
    FILE *fp = fopen(filepath, "r");
    if (fp == NULL) {
        printf("ERROR: Cannot open LX-Mie file: %s\n", filepath);
        return -1;
    }
    
    // Skip header line (starts with #)
    char header[1024];
    if (fgets(header, sizeof(header), fp) == NULL) {
        printf("ERROR: Failed to read header from %s\n", filepath);
        fclose(fp);
        return -1;
    }
    
    // Read data lines
    // Format: wavelength(μm)  size_param  extinction(cm²)  scattering(cm²)  absorption(cm²)  albedo  asymmetry_g
    int count = 0;
    double wl, size_param, ext, scat, abs, alb, asym_val;
    
    while (fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", 
                 &wl, &size_param, &ext, &scat, &abs, &alb, &asym_val) == 7) {
        if (count >= MAX_MIE_WAVELENGTHS) {
            printf("WARNING: Too many wavelengths in %s, truncating at %d\n", 
                   filepath, MAX_MIE_WAVELENGTHS);
            break;
        }
        
        wavelengths[count] = wl;           // Column 0: wavelength in microns
        extinction[count] = ext;           // Column 2: extinction cross section (cm²) - total = scattering + absorption
        albedo[count] = alb;               // Column 5: single scattering albedo
        asymmetry[count] = asym_val;       // Column 6: asymmetry parameter g
        
        count++;
    }
    
    fclose(fp);
    
    if (count == 0) {
        printf("ERROR: No data read from %s\n", filepath);
        return -1;
    }
    
    *n_wavelengths = count;
    return 0;
}

/**
 * Generate EPACRIS wavelength grid
 * w(0) = 0.1 μm
 * w(i) = w(i-1) * (1+1/200) for i=1 to 1386
 * Total: 1387 wavelengths (approximately 0.1 to 100 μm)
 * 
 * Note: This formula gives w(1386) ≈ 100.499 μm, not exactly 100 μm.
 * This matches the EPACRIS file format specification.
 */
static void generate_epacris_wavelength_grid(double *wavelength, int n_wavelengths) {
    wavelength[0] = 0.1;  // w(0) = 0.1 μm
    for (int i = 1; i < n_wavelengths; i++) {

        wavelength[i] = wavelength[i-1] * (1.0 + 1.0/200.0);  // Constant resolution: 1+1/200
    }
}

/**
 * Read EPACRIS format MIE table file (Albedo, Cross, or Geo)
 * 
 * File format:
 * - First column: particle radius (log-spaced, starting at 0.01)
 * - Remaining columns: wavelength-dependent values (1387 columns)
 * - No header line
 * 
 * @param filepath Path to the file
 * @param data Output array [n_particle_sizes][n_wavelengths] to store values
 * @param radius_out Output array [n_particle_sizes] to store radii from first column
 * @param n_particle_sizes Number of particle sizes (rows) - will be determined from file
 * @param n_wavelengths Number of wavelengths (columns - 1, since first column is radius)
 * @return Number of particle sizes read on success, -1 on error
 */
static int read_epacris_mie_file(const char *filepath, double **data, double *radius_out, 
                                 int max_particle_sizes, int n_wavelengths) {
    FILE *fp = fopen(filepath, "r");
    if (fp == NULL) {
        printf("ERROR: Cannot open EPACRIS MIE file: %s\n", filepath);
        return -1;
    }
    
    // Buffer needs to be large enough for 1388 columns (1 radius + 1387 wavelengths)
    // Each value can be up to ~15 characters (e.g., "1.234567e-10"), plus spaces
    // 1388 * 20 = ~27760 bytes, use 32KB to be safe
    char line[32768];
    int row = 0;
    
    while (fgets(line, sizeof(line), fp) != NULL && row < max_particle_sizes) {
        // Check if line was truncated (fgets doesn't include newline if buffer is full)
        size_t len = strlen(line);
        if (len > 0 && line[len-1] != '\n' && !feof(fp)) {
            printf("WARNING: Line %d in %s was truncated (buffer too small)\n", row+1, filepath);
        }
        
        // Parse the line: first value is radius, rest are wavelength values
        char *token = strtok(line, " \t\n");
        if (token == NULL) continue;
        
        // Read radius from first column
        radius_out[row] = atof(token);
        
        // Read wavelength values (columns 1 to n_wavelengths)
        int col = 0;
        token = strtok(NULL, " \t\n");  // Start reading from second column
        
        while (token != NULL && col < n_wavelengths) {
            data[row][col] = atof(token);
            col++;
            token = strtok(NULL, " \t\n");
        }
        
        if (col != n_wavelengths) {
            printf("WARNING: Row %d in %s has %d wavelength columns (expected %d)\n", 
                   row+1, filepath, col, n_wavelengths);
        }
        
        row++;
    }
    
    fclose(fp);
    
    if (row == 0) {
        printf("ERROR: No data read from %s\n", filepath);
        return -1;
    }
    
    return row;  // Return actual number of particle sizes read
}

/**
 * Read EPACRIS format cloud optical tables for a single species
 * 
 * Reads 3 files: Albedo.dat, Cross.dat, Geo.dat
 * File naming convention: Just Albedo.dat, Cross.dat, Geo.dat (no species prefix)
 * Directory structure: CLOUD_MIE_DIRECTORY_EPACRIS/SPECIES_NAME/ (e.g., EPACRIS_MIE/H2O/)
 * 
 * @param species_id Species ID (e.g., 7 for H2O)
 * @param species_name Species name (e.g., "H2O" or "NH3")
 * @param directory Directory containing the files (e.g., "EPACRIS_MIE/H2O")
 * @param optics Output structure to populate
 * @return 0 on success, -1 on error
 */
static int read_epacris_species_tables(int species_id, const char *species_name, 
                                       const char *directory, CloudOpticalTableMie *optics) {
    // EPACRIS format constants
    const int EPACRIS_N_WAVELENGTHS = 1387;
    
    if (EPACRIS_N_WAVELENGTHS > MAX_MIE_WAVELENGTHS) {
        printf("ERROR: EPACRIS format requires %d wavelengths, but MAX_MIE_WAVELENGTHS = %d\n",
               EPACRIS_N_WAVELENGTHS, MAX_MIE_WAVELENGTHS);
        return -1;
    }
    
    // Build file paths - files are just Albedo.dat, Cross.dat, Geo.dat
    char actual_albedo[1024], actual_cross[1024], actual_geo[1024];
    snprintf(actual_albedo, sizeof(actual_albedo), "%s/Albedo.dat", directory);
    snprintf(actual_cross, sizeof(actual_cross), "%s/Cross.dat", directory);
    snprintf(actual_geo, sizeof(actual_geo), "%s/Geo.dat", directory);
    
    // Check if files exist
    FILE *test_file = fopen(actual_albedo, "r");
    if (test_file == NULL) {
        printf("ERROR: Cannot open EPACRIS Albedo file: %s\n", actual_albedo);
        return -1;
    }
    fclose(test_file);
    
    test_file = fopen(actual_cross, "r");
    if (test_file == NULL) {
        printf("ERROR: Cannot open EPACRIS Cross file: %s\n", actual_cross);
        return -1;
    }
    fclose(test_file);
    
    test_file = fopen(actual_geo, "r");
    if (test_file == NULL) {
        printf("ERROR: Cannot open EPACRIS Geo file: %s\n", actual_geo);
        return -1;
    }
    fclose(test_file);
    
    printf("Found EPACRIS files:\n");
    printf("  Albedo: %s\n", actual_albedo);
    printf("  Cross: %s\n", actual_cross);
    printf("  Geo: %s\n", actual_geo);
    
    // Initialize wavelength grid first (needed for reading)
    optics->n_wavelengths = EPACRIS_N_WAVELENGTHS;
    generate_epacris_wavelength_grid(optics->wavelength, optics->n_wavelengths);
    
    // Allocate temporary arrays for reading (use MAX_PARTICLE_SIZES as upper bound)
    double **albedo_data = NULL;
    double **cross_data = NULL;
    double **geo_data = NULL;
    double *radius_temp = NULL;
    
    albedo_data = (double **)malloc(MAX_PARTICLE_SIZES * sizeof(double *));
    cross_data = (double **)malloc(MAX_PARTICLE_SIZES * sizeof(double *));
    geo_data = (double **)malloc(MAX_PARTICLE_SIZES * sizeof(double *));
    radius_temp = (double *)malloc(MAX_PARTICLE_SIZES * sizeof(double));
    
    if (albedo_data == NULL || cross_data == NULL || geo_data == NULL || radius_temp == NULL) {
        printf("ERROR: Memory allocation failed\n");
        goto cleanup_error;
    }
    
    for (int i = 0; i < MAX_PARTICLE_SIZES; i++) {
        albedo_data[i] = (double *)malloc(optics->n_wavelengths * sizeof(double));
        cross_data[i] = (double *)malloc(optics->n_wavelengths * sizeof(double));
        geo_data[i] = (double *)malloc(optics->n_wavelengths * sizeof(double));
    }
    
    // Read the three files - they will determine the actual number of particle sizes
    int n_particle_sizes_albedo = read_epacris_mie_file(actual_albedo, albedo_data, radius_temp, 
                                                        MAX_PARTICLE_SIZES, optics->n_wavelengths);
    if (n_particle_sizes_albedo < 0) {
        printf("ERROR: Failed to read Albedo file\n");
        goto cleanup_error;
    }
    
    int n_particle_sizes_cross = read_epacris_mie_file(actual_cross, cross_data, radius_temp, 
                                                        MAX_PARTICLE_SIZES, optics->n_wavelengths);
    if (n_particle_sizes_cross < 0) {
        printf("ERROR: Failed to read Cross file\n");
        goto cleanup_error;
    }
    
    int n_particle_sizes_geo = read_epacris_mie_file(actual_geo, geo_data, radius_temp, 
                                                     MAX_PARTICLE_SIZES, optics->n_wavelengths);
    if (n_particle_sizes_geo < 0) {
        printf("ERROR: Failed to read Geo file\n");
        goto cleanup_error;
    }
    
    // Verify all files have the same number of particle sizes
    if (n_particle_sizes_albedo != n_particle_sizes_cross || 
        n_particle_sizes_albedo != n_particle_sizes_geo) {
        printf("ERROR: Mismatch in particle size counts: Albedo=%d, Cross=%d, Geo=%d\n",
               n_particle_sizes_albedo, n_particle_sizes_cross, n_particle_sizes_geo);
        goto cleanup_error;
    }
    
    // Set the actual number of particle sizes
    optics->n_particle_sizes = n_particle_sizes_albedo;
    
    // Copy radii from first file (they should be the same in all files)
    for (int i = 0; i < optics->n_particle_sizes; i++) {
        optics->radius[i] = radius_temp[i];
    }
    
    free(radius_temp);
    
    // Convert EPACRIS format to CloudOpticalTableMie structure
    // Cross file contains extinction cross section (cm²)
    // Albedo file contains single scattering albedo
    // Geo file contains asymmetry parameter g
    for (int i = 0; i < optics->n_particle_sizes; i++) {
        for (int j = 0; j < optics->n_wavelengths; j++) {
            optics->extinction_cross[i][j] = cross_data[i][j];
            optics->albedo[i][j] = albedo_data[i][j];
            optics->asymmetry_g[i][j] = geo_data[i][j];
        }
    }
    
    // Free temporary arrays
    for (int i = 0; i < optics->n_particle_sizes; i++) {
        free(albedo_data[i]);
        free(cross_data[i]);
        free(geo_data[i]);
    }
    free(albedo_data);
    free(cross_data);
    free(geo_data);
    
    printf("Successfully loaded EPACRIS format tables:\n");
    printf("  Particle sizes: %d (%.6f - %.6f μm)\n", 
           optics->n_particle_sizes, 
           optics->radius[0], 
           optics->radius[optics->n_particle_sizes-1]);
    printf("  Wavelengths: %d (%.6f - %.6f μm)\n",
           optics->n_wavelengths,
           optics->wavelength[0],
           optics->wavelength[optics->n_wavelengths-1]);
    
    return 0;
    
cleanup_error:
    // Free temporary arrays on error
    if (albedo_data != NULL) {
        for (int i = 0; i < MAX_PARTICLE_SIZES; i++) {
            if (albedo_data[i] != NULL) free(albedo_data[i]);
            if (cross_data[i] != NULL) free(cross_data[i]);
            if (geo_data[i] != NULL) free(geo_data[i]);
        }
        free(albedo_data);
        free(cross_data);
        free(geo_data);
    }
    if (radius_temp != NULL) {
        free(radius_temp);
    }
    return -1;
}

/**
 * Read all EPACRIS format cloud optical tables for cloud species (defined in config.h)
 * Reads Albedo.dat, Cross.dat, Geo.dat files for each species
 */
static void read_cloud_optical_tables_epacris(void) {
    printf("%s\n",fillmi); 
    printf("READING CLOUD OPTICAL PROPERTY TABLES (EPACRIS FORMAT)\n");
    printf("%s\n",fillmi); 
    
    // Initialize counter
    n_cloud_species_loaded = 0;
    
    // Build the cloud species list from CLOUD_SPECIES_LIST macro
    static const int cloud_species_list[] = {CLOUD_SPECIES_LIST};
    const int n_cloud_species = sizeof(cloud_species_list) / sizeof(cloud_species_list[0]);
    
    // Loop through all cloud species
    for (int species_idx = 0; species_idx < n_cloud_species && species_idx < MAX_CLOUD_SPECIES; species_idx++) {
        int species_id = cloud_species_list[species_idx];
        
        // Get species name from ID
        const char *species_name = get_species_name(species_id);
        if (strcmp(species_name, "Unknown") == 0) {
            printf("WARNING: CLOUD_SPECIES[%d] (%d) is not a recognized species ID. Skipping.\n", species_idx, species_id);
            continue;
        }
        
        printf("--- Processing cloud species %d/%d: %s (ID: %d) ---\n", 
               species_idx + 1, n_cloud_species, species_name, species_id);
        
        // Build directory path: CLOUD_MIE_DIRECTORY_EPACRIS/SPECIES_NAME
        // Files are in subdirectories like EPACRIS_MIE/H2O/, EPACRIS_MIE/NH3/
        char directory[1024];
        snprintf(directory, sizeof(directory), "%s/%s", CLOUD_MIE_DIRECTORY_EPACRIS, species_name);
        
        printf("Reading cloud optical tables for species: %s (ID: %d)\n", species_name, species_id);
        printf("Directory: %s\n", directory);
        
        // Initialize structure for this species
        CloudOpticalTableMie *optics = &cloud_mie_optics[n_cloud_species_loaded];
        cloud_species_ids[n_cloud_species_loaded] = species_id;
        
        // Files are named Albedo.dat, Cross.dat, Geo.dat (no species prefix)
        // Pass species_name as-is (e.g., "H2O", "NH3") - it's only used for error messages
        if (read_epacris_species_tables(species_id, species_name, directory, optics) != 0) {
            printf("ERROR: Failed to load EPACRIS tables for %s. Skipping.\n", species_name);
            continue;
        }
        
        // Increment counter for successfully loaded species
        n_cloud_species_loaded++;
        printf("Successfully loaded optical tables for %s\n", species_name);
    }
    
    printf("%s\n",fillmi); 
    printf("Loaded optical tables for %d cloud species\n", n_cloud_species_loaded);
    printf("%s\n",fillmi); 
    
    // Print summary for each loaded species
    for (int i = 0; i < n_cloud_species_loaded; i++) {
        const char *species_name = get_species_name(cloud_species_ids[i]);
        CloudOpticalTableMie *optics = &cloud_mie_optics[i];
        printf("Species %d: %s (ID: %d)\n", i+1, species_name, cloud_species_ids[i]);
        printf("  Particle size range: %.6f - %.6f micrometers\n", 
               optics->radius[0], 
               optics->radius[optics->n_particle_sizes-1]);
        printf("  Number of particle sizes: %d\n", optics->n_particle_sizes);
        printf("  Number of wavelengths: %d\n", optics->n_wavelengths);
        printf("  Wavelength range: %.6f - %.6f micrometers\n",
               optics->wavelength[0],
               optics->wavelength[optics->n_wavelengths-1]);
    }
    printf("%s\n",fillmi); 
}

//--------------------------------------------------------------------- 
//  Read all MIE output files for cloud species (defined in config.h)
//  Main function called from epacris_main.c
//--------------------------------------------------------------------- 
void read_cloud_optical_tables_mie(void) {
    // Check which format to use
    // Currently can use LX-MIE (like in HELIOS) or EPACRIS format for MIE tables

    #if USE_EPACRIS_FORMAT
        read_cloud_optical_tables_epacris();
        return;
    #endif
    
    printf("%s\n",fillmi); 
    printf("READING CLOUD OPTICAL PROPERTY TABLES (LX-MIE FORMAT)\n");
    printf("%s\n",fillmi); 
    
    // Initialize counter
    n_cloud_species_loaded = 0;
    
    // Build the cloud species list from CLOUD_SPECIES_LIST macro
    // This automatically handles any number of species defined in config.h
    static const int cloud_species_list[] = {CLOUD_SPECIES_LIST};
    const int n_cloud_species = sizeof(cloud_species_list) / sizeof(cloud_species_list[0]);
    
    // Loop through all cloud species
    for (int species_idx = 0; species_idx < n_cloud_species && species_idx < MAX_CLOUD_SPECIES; species_idx++) {
        int species_id = cloud_species_list[species_idx];
        
        // Get species name from ID
        const char *species_name = get_species_name(species_id);
        if (strcmp(species_name, "Unknown") == 0) {
            printf("WARNING: CLOUD_SPECIES[%d] (%d) is not a recognized species ID. Skipping.\n", species_idx, species_id);
            continue;
        }
        
        printf("\n--- Processing cloud species %d/%d: %s (ID: %d) ---\n", 
               species_idx + 1, n_cloud_species, species_name, species_id);
        
        // Build directory path: CLOUD_MIE_DIRECTORY_LXMIE/SPECIES_NAME
        char directory[512];
        snprintf(directory, sizeof(directory), "%s/%s", CLOUD_MIE_DIRECTORY_LXMIE, species_name);
        
        
        printf("Reading cloud optical tables for species: %s (ID: %d)\n", species_name, species_id);
        printf("Directory: %s\n", directory);
        
        // Initialize structure for this species
        CloudOpticalTableMie *optics = &cloud_mie_optics[n_cloud_species_loaded];
        optics->n_particle_sizes = 0;
        optics->n_wavelengths = 0;
        cloud_species_ids[n_cloud_species_loaded] = species_id;
    
    // Open directory and collect all r*.dat filenames
    DIR *dir = opendir(directory);
    if (dir == NULL) {
        printf("FATAL ERROR: Cannot open directory: %s\n", directory);
        exit(1);
    }
    
    // Collect all r*.dat filenames
    char **filenames = (char **)malloc(MAX_PARTICLE_SIZES * sizeof(char *));
    int file_count = 0;
    
    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL && file_count < MAX_PARTICLE_SIZES) {
        const char *name = entry->d_name;
        
        // Check if filename matches pattern r*.dat
        if (name[0] == 'r' && strstr(name, ".dat") != NULL) {
            // Check if radius extraction works
            double radius = extract_radius_from_filename(name);
            if (radius > 0.0) {
                filenames[file_count] = (char *)malloc(strlen(name) + 1);
                strcpy(filenames[file_count], name);
                file_count++;
            }
        }
    }
    closedir(dir);
    
    if (file_count == 0) {
        printf("FATAL ERROR: No valid r*.dat files found in %s\n", directory);
        exit(1);
    }
    
    // Sort filenames by particle radius (ascending)
    qsort(filenames, file_count, sizeof(char *), compare_particle_files);
    
    printf("Found %d particle size files\n", file_count);
    
    // Read first file to determine number of wavelengths
    char first_filepath[1024];
    snprintf(first_filepath, sizeof(first_filepath), "%s/%s", directory, filenames[0]);
    
    //Store temp var for reading out
    double temp_wavelengths[MAX_MIE_WAVELENGTHS];
    double temp_extinction[MAX_MIE_WAVELENGTHS];
    double temp_albedo[MAX_MIE_WAVELENGTHS];
    double temp_asymmetry[MAX_MIE_WAVELENGTHS];
    
        if (read_mie_file(first_filepath, temp_wavelengths, temp_extinction, 
                         temp_albedo, temp_asymmetry, &optics->n_wavelengths) != 0) {
            printf("FATAL ERROR: Failed to read first file: %s\n", first_filepath);
            exit(1);
        }
        
        // Copy wavelengths from first file (should be same for all particle sizes)
        for (int j = 0; j < optics->n_wavelengths; j++) {
            optics->wavelength[j] = temp_wavelengths[j];
        }
        
        printf("Number of wavelengths per file: %d\n", optics->n_wavelengths);
        
        // Read all particle size files
        for (int i = 0; i < file_count; i++) {
            char filepath[1024];
            snprintf(filepath, sizeof(filepath), "%s/%s", directory, filenames[i]);
            
            double radius = extract_radius_from_filename(filenames[i]);
            optics->radius[i] = radius;
            
            if (read_mie_file(filepath, temp_wavelengths, temp_extinction,
                             temp_albedo, temp_asymmetry, &optics->n_wavelengths) != 0) {
                printf("WARNING: Failed to read %s, skipping\n", filepath);
                continue;
            }
            
            // Verify wavelength grid matches
            int wavelengths_match = 1;
            for (int j = 0; j < optics->n_wavelengths; j++) {
                if (fabs(temp_wavelengths[j] - optics->wavelength[j]) > 1e-10) {
                    wavelengths_match = 0;
                    break;
                }
            }
            
            if (!wavelengths_match) {
                printf("WARNING: Wavelength grid mismatch in %s, skipping\n", filepath);
                continue;
            }
            
            // Copy optical properties, i is the particle size index, j is the wavelength index
            for (int j = 0; j < optics->n_wavelengths; j++) {
                optics->extinction_cross[i][j] = temp_extinction[j];
                optics->albedo[i][j] = temp_albedo[j];
                optics->asymmetry_g[i][j] = temp_asymmetry[j];
            }
            
            optics->n_particle_sizes++;
        }
        
        // Print debug in scientific notation
        printf("n_particle_sizes = %d\n", optics->n_particle_sizes);
        printf("n_wavelengths = %d\n", optics->n_wavelengths);
        printf("radius[0] = %e\n", optics->radius[0]);
        printf("wavelength[0] = %e\n", optics->wavelength[0]);
        printf("extinction_cross[0][0] = %e\n", optics->extinction_cross[0][0]);
        printf("albedo[0][0] = %e\n", optics->albedo[0][0]);
        printf("asymmetry_g[0][0] = %e\n", optics->asymmetry_g[0][0]);
        printf("extinction_cross[0][1] = %e\n", optics->extinction_cross[0][1]);
        printf("albedo[0][1] = %e\n", optics->albedo[0][1]);
        printf("asymmetry_g[0][1] = %e\n", optics->asymmetry_g[0][1]);
        printf("extinction_cross[0][2] = %e\n", optics->extinction_cross[0][2]);
        
        // Free filenames memory
        for (int i = 0; i < file_count; i++) {
            free(filenames[i]);
        }
        free(filenames);
        
        // Increment counter for successfully loaded species
        n_cloud_species_loaded++;
        printf("Successfully loaded optical tables for %s\n", species_name);
    }
    
    printf("%s\n",fillmi); 
    printf("Loaded optical tables for %d cloud species\n", n_cloud_species_loaded);
    printf("%s\n",fillmi); 
    
    // Print summary for each loaded species
    for (int i = 0; i < n_cloud_species_loaded; i++) {
        const char *species_name = get_species_name(cloud_species_ids[i]);
        CloudOpticalTableMie *optics = &cloud_mie_optics[i];
        printf("\nSpecies %d: %s (ID: %d)\n", i+1, species_name, cloud_species_ids[i]);
        printf("  Particle size range: %.6f - %.6f micrometers\n", 
               optics->radius[0], 
               optics->radius[optics->n_particle_sizes-1]);
        printf("  Number of particle sizes: %d\n", optics->n_particle_sizes);
        printf("  Number of wavelengths: %d\n", optics->n_wavelengths);
        printf("  Wavelength range: %.6f - %.6f micrometers\n",
               optics->wavelength[0],
               optics->wavelength[optics->n_wavelengths-1]);
    }
    printf("%s\n",fillmi); 
}


/**
 * Linear interpolation in log-space (for particle size interpolation)
 * Interpolates y (linear) vs log10(x) - used for albedo and asymmetry
 * Similar to MATLAB interp1(log10(x), y, log10(xi))
 */
static double interp1_log(double *x, double *y, int n, double xi) {
    if (n < 2) return y[0];
    
    double log_xi = log10(fmax(x[0], fmin(xi, x[n-1])));
    
    // Handle boundaries
    if (log_xi <= log10(x[0])) return y[0];
    if (log_xi >= log10(x[n-1])) return y[n-1];
    
    // Find interpolation interval
    for (int i = 0; i < n-1; i++) {
        double log_x0 = log10(x[i]);
        double log_x1 = log10(x[i+1]);
        
        if (log_xi >= log_x0 && log_xi <= log_x1) {
            // Linear interpolation in log space
            double t = (log_xi - log_x0) / (log_x1 - log_x0);
            return y[i] + t * (y[i+1] - y[i]);
        }
    }
    
    return y[n-1]; // fallback
}

/**
 * Log-log interpolation (for cross-section interpolation)
 * Interpolates log10(y) vs log10(x), then converts back
 * Similar to MATLAB pow(10, interp1(log10(x), log10(y), log10(xi)))
 * Matches cloud_interp.c line 167: pow(10, interp1(log10(H2OI_r), log10(H2OI_c), log10(...)))
 */
static double interp1_loglog(double *x, double *y, int n, double xi) {
    if (n < 2) return y[0]; // check if fewer than 2 points, no interpolation needed
    
    // convert r0 to log10 space while maintining boundaries
    double log_xi = log10(fmax(x[0], fmin(xi, x[n-1])));
    
    // If boundary conditions are met, return the value at the boundary
    if (log_xi <= log10(x[0])) return y[0];
    if (log_xi >= log10(x[n-1])) return y[n-1];
    
    // Find interpolation interval
    // Loop through all particle sizes
    for (int i = 0; i < n-1; i++) {
        double log_x0 = log10(x[i]);
        double log_x1 = log10(x[i+1]);
        
        // If between boundary conditions, perform linear interpolation in log-log space
        if (log_xi >= log_x0 && log_xi <= log_x1) {
            // Linear interpolation in log-log space
            double log_y0 = log10(fmax(1e-30, y[i]));      // Protect against zero/negative
            double log_y1 = log10(fmax(1e-30, y[i+1]));
            double t = (log_xi - log_x0) / (log_x1 - log_x0); // calculate the weight factor
            double log_y_interp = log_y0 + t * (log_y1 - log_y0);
            return pow(10.0, log_y_interp);
        }
    }
    
    return y[n-1]; // fallback
}

/**
 * Linear interpolation in wavelength (for wavelength grid mapping)
 */
static double interp1_wavelength(double *wavelength_mie, double *value_mie, int n_mie,
                                double wavelength_rt) {
    if (n_mie < 2) return value_mie[0];
    
    // Handle boundaries
    if (wavelength_rt <= wavelength_mie[0]) return value_mie[0];
    if (wavelength_rt >= wavelength_mie[n_mie-1]) return value_mie[n_mie-1];
    
    // Find interpolation interval
    for (int i = 0; i < n_mie-1; i++) {
        if (wavelength_rt >= wavelength_mie[i] && wavelength_rt <= wavelength_mie[i+1]) {
            double t = (wavelength_rt - wavelength_mie[i]) / (wavelength_mie[i+1] - wavelength_mie[i]);
            return value_mie[i] + t * (value_mie[i+1] - value_mie[i]);
        }
    }
    
    return value_mie[n_mie-1]; // fallback
}

/**
 * Helper function to get output array pointers for a given species ID
 * Returns 1 on success, 0 on failure
 */
static int get_cloud_output_arrays(int species_id, double ***c_out, double ***a_out, double ***g_out) {
    if (species_id == 7) {  // H2O
        *c_out = cH2O;
        *a_out = aH2O;
        *g_out = gH2O;
        return 1;
    } else if (species_id == 9) {  // NH3
        *c_out = cNH3;
        *a_out = aNH3;
        *g_out = gNH3;
        return 1;
    }
    return 0;  // Unsupported species
}

/**
 * Calculate cloud optical properties for all layers and wavelengths
 * 
 * This function processes ALL loaded cloud species:
 * 1. Interpolates Mie scattering data in particle size dimension (log-space)
 * 2. Interpolates in wavelength dimension to match EPACRIS RT grid
 * 3. Converts cloud mass density to opacity (cm^-1)
 * 
 * Inputs:
 * - particle_r2: particle radii from cloud physics (micrometers, volume-weighted) [zbin+1][MAX_CONDENSIBLES]
 * - clouds: cloud mass densities (molecules/cm³) [zbin+1][NSP+1]
 * - wavelength: EPACRIS wavelength grid (nm, from wavelength[] array)
 * - cH2O, aH2O, gH2O, cNH3, aNH3, gNH3: output arrays [zbin+1][NLAMBDA]
 * 
 * Based on cloud_interp.c logic with HELIOS Mie data format
 */
void calculate_cloud_opacity_arrays(void) {
    // All arrays accessed as globals:
    // Inputs: clouds, particle_r2, particle_number_density, wavelength
    // Outputs: cH2O, aH2O, gH2O, cNH3, aNH3, gNH3
        
    // Loop through all loaded cloud species
    for (int species_table_idx = 0; species_table_idx < n_cloud_species_loaded; species_table_idx++) {
        int cloud_species_id = cloud_species_ids[species_table_idx];
        CloudOpticalTableMie *optics = &cloud_mie_optics[species_table_idx];
        
        // Get output array pointers for this species
        double **c_out, **a_out, **g_out;
        if (!get_cloud_output_arrays(cloud_species_id, &c_out, &a_out, &g_out)) {
            const char *species_name = get_species_name(cloud_species_id);
            printf("WARNING: Cloud species %s (ID: %d) not supported for output arrays. Skipping.\n", 
                   species_name, cloud_species_id);
            continue;
        }
        
        // Temporary arrays for interpolation
        double *ext_cross_interp = (double *)malloc(optics->n_wavelengths * sizeof(double));
        double *albedo_interp = (double *)malloc(optics->n_wavelengths * sizeof(double));
        double *asym_interp = (double *)malloc(optics->n_wavelengths * sizeof(double));
        
        // Find cloud species index in CONDENSIBLES array
        int cloud_species_idx = -1;
        for (int k = 0; k < NCONDENSIBLES; k++) {
            if (CONDENSIBLES[k] == cloud_species_id) {
                cloud_species_idx = k;
                break;
            }
        }
        
        // Check if cloud species is in condensibles list
        if (cloud_species_idx < 0) {
            // Cloud species not in condensibles list - skip this species
            const char *species_name = get_species_name(cloud_species_id);
            printf("WARNING: Cloud species %s (ID: %d) not in condensibles list. Skipping cloud opacity calculation.\n", 
                   species_name, cloud_species_id);
            free(ext_cross_interp);
            free(albedo_interp);
            free(asym_interp);
            continue;
        }
        
        // Process each layer for this species
        for (int j = 1; j <= zbin; j++) {
            // if cloud density very small, set to zero opacity and continue to next layer
            if (clouds[j][cloud_species_id] < 1e-12) {
                // No clouds or very small particles - set to zero opacity
                for (int i = 0; i < NLAMBDA; i++) {
                    c_out[j][i] = 0.0;
                    a_out[j][i] = 1.0;
                    g_out[j][i] = 0.0;
                }
                continue;
            }
            
            // Use volume-weighted radius (r2, 2nd moment) for cloud optics interpolation
            // This matches the EPACRIS format specification and is used for both formats
            // Note: Interpolation functions handle out-of-range values automatically
            // by clamping to the actual lookup table range
            double particle_radius = particle_r2[j][cloud_species_idx];
            
            // Particle number density: stored in particles/m³, convert to particles/cm³
            double particle_number_density_cm3 = particle_number_density[j][cloud_species_idx] * 1.0e-6;  // particles/m³ -> particles/cm³
            
            // Debug: Print interpolation parameters only for layers with non-zero particle size
            // Check both particle_radius and particle number density to ensure we only print for actual clouds
#if CLOUD_DEBUG_MIE
            if (particle_radius > 1e-10 && particle_number_density[j][cloud_species_idx] > 1e-10) {
                printf("\n=== CLOUD OPTICS INTERPOLATION DEBUG (Layer %d) ===\n", j);
                printf("particle_r2[%d][%d] = %.6e μm (volume-weighted, 2nd moment)\n", j, cloud_species_idx, particle_radius);
                printf("Lookup table radius range: %.6e - %.6e μm\n", 
                       optics->radius[0], 
                       optics->radius[optics->n_particle_sizes-1]);
                printf("Number of particle sizes in lookup table: %d\n", optics->n_particle_sizes);
                printf("particle_number_density[%d][%d] = %.6e particles/m³\n", j, cloud_species_idx, particle_number_density[j][cloud_species_idx]);
                printf("particle_number_density_cm3 = %.6e particles/cm³\n", particle_number_density_cm3);
            }
#endif


            //////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Interpolate optical properties in particle size dimension
            // Following cloud_interp.c rules (matching MATLAB cloud_interp.m):
            // - Extinction cross-section: log-log interpolation (log10 vs log10) - matches H2OI_c interpolation
            // - Albedo: linear in log-radius space - matches H2OI_a interpolation  
            // - Asymmetry: linear in log-radius space - matches H2OI_g interpolation
            
            // Interpolate for each Mie wavelength
            for (int w = 0; w < optics->n_wavelengths; w++) {
                
                // Extract values for this wavelength across all particle sizes
                double *ext_per_r = (double *)malloc(optics->n_particle_sizes * sizeof(double));
                double *alb_per_r = (double *)malloc(optics->n_particle_sizes * sizeof(double));
                double *asym_per_r = (double *)malloc(optics->n_particle_sizes * sizeof(double));
                
                for (int r = 0; r < optics->n_particle_sizes; r++) {
                    ext_per_r[r] = optics->extinction_cross[r][w];
                    alb_per_r[r] = optics->albedo[r][w];
                    asym_per_r[r] = optics->asymmetry_g[r][w];
                }
                
                // Extinction cross-section: log-log interpolation (matches cloud_interp.c line 167: H2OI_c)
                ext_cross_interp[w] = interp1_loglog(optics->radius, ext_per_r, 
                                                       optics->n_particle_sizes, particle_radius);
                
                // Albedo: linear in log-radius space (matches cloud_interp.c line 168: H2OI_a)
                albedo_interp[w] = interp1_log(optics->radius, alb_per_r,
                                               optics->n_particle_sizes, particle_radius);
                
                // Asymmetry: linear in log-radius space (matches cloud_interp.c line 169: H2OI_g)
                asym_interp[w] = interp1_log(optics->radius, asym_per_r,
                                            optics->n_particle_sizes, particle_radius);
            
                // Debug: Print interpolation results for layers with non-zero particle size (first 3 wavelengths)
#if CLOUD_DEBUG_MIE
                if (particle_radius > 1e-10 && particle_number_density[j][cloud_species_idx] > 1e-10 && w < 3) {
                    // Find bracketing indices for particle_radius
                    int idx_low = -1, idx_high = -1;
                    for (int r = 0; r < optics->n_particle_sizes - 1; r++) {
                        if (particle_radius >= optics->radius[r] && particle_radius <= optics->radius[r+1]) {
                            idx_low = r;
                            idx_high = r + 1;
                            break;
                        }
                    }
                    
                    printf("\n  Wavelength %d (%.6e μm):\n", w, optics->wavelength[w]);
                    if (idx_low >= 0) {
                        printf("    Lookup table radii: %.6e - %.6e μm (bracketing radius=%.6e)\n",
                               optics->radius[idx_low], optics->radius[idx_high], particle_radius);
                        printf("    Extinction at radii: %.6e, %.6e cm²\n",
                               ext_per_r[idx_low], ext_per_r[idx_high]);
                        printf("    Albedo at radii: %.6f, %.6f\n",
                               alb_per_r[idx_low], alb_per_r[idx_high]);
                        printf("    Asymmetry at radii: %.6f, %.6f\n",
                               asym_per_r[idx_low], asym_per_r[idx_high]);
                    } else if (particle_radius <= optics->radius[0]) {
                        printf("    particle_radius below lookup range - using first table value\n");
                        printf("    Extinction: %.6e cm² (from radius %.6e μm)\n",
                               ext_per_r[0], optics->radius[0]);
                    } else {
                        printf("    particle_radius above lookup range - using last table value\n");
                        printf("    Extinction: %.6e cm² (from radius %.6e μm)\n",
                               ext_per_r[optics->n_particle_sizes-1],
                               optics->radius[optics->n_particle_sizes-1]);
                    }
                    printf("    Interpolated extinction: %.6e cm²\n", ext_cross_interp[w]);
                    printf("    Interpolated albedo: %.6f\n", albedo_interp[w]);
                    printf("    Interpolated asymmetry: %.6f\n", asym_interp[w]);
                }
#endif
                
                free(ext_per_r);
                free(alb_per_r);
                free(asym_per_r);
            }  // End Mie wavelength loop
            
            // Debug: Print summary for layers with non-zero particle size
#if CLOUD_DEBUG_MIE
            if (particle_radius > 1e-10 && particle_number_density[j][cloud_species_idx] > 1e-10) {
                printf("\n  Interpolation complete for layer %d\n", j);
                printf("  Sample interpolated values (first 3 wavelengths):\n");
                for (int w = 0; w < 3 && w < optics->n_wavelengths; w++) {
                    printf("    λ=%.6e μm: ext=%.6e cm², alb=%.6f, asym=%.6f\n",
                           optics->wavelength[w],
                           ext_cross_interp[w], albedo_interp[w], asym_interp[w]);
                }
                printf("%s\n",fillmi); 
            }
#endif
            
            
            // Now interpolate in wavelength dimension to EPACRIS RT grid
            for (int i = 0; i < NLAMBDA; i++) {
            double wavelength_rt_micron = wavelength[i] * 1.0e-3;  // Convert nm to microns
            
                // Interpolate optical properties to RT wavelength grid
                double ext_cross = interp1_wavelength(optics->wavelength, ext_cross_interp,
                                                      optics->n_wavelengths, wavelength_rt_micron);
                double alb_val = interp1_wavelength(optics->wavelength, albedo_interp,
                                                   optics->n_wavelengths, wavelength_rt_micron);
                double asym_val = interp1_wavelength(optics->wavelength, asym_interp,
                                                    optics->n_wavelengths, wavelength_rt_micron);
                
                // Debug: Print wavelength interpolation results for layers with non-zero particle size (first 5 wavelengths)
#if CLOUD_DEBUG_MIE
                if (particle_radius > 1e-10 && particle_number_density[j][cloud_species_idx] > 1e-10 && i < 5) {
                // Find bracketing indices for wavelength
                int idx_low = -1, idx_high = -1;
                for (int w = 0; w < optics->n_wavelengths - 1; w++) {
                    if (wavelength_rt_micron >= optics->wavelength[w] && 
                        wavelength_rt_micron <= optics->wavelength[w+1]) {
                        idx_low = w;
                        idx_high = w + 1;
                        break;
                    }
                }
                
                printf("\n  RT Wavelength %d (%.6e μm = %.6e nm):\n", i, wavelength_rt_micron, wavelength[i]);
                if (idx_low >= 0) {
                    printf("    MIE lookup wavelengths: %.6e - %.6e μm (bracketing λ_rt=%.6e)\n",
                           optics->wavelength[idx_low], optics->wavelength[idx_high], wavelength_rt_micron);
                    printf("    Extinction at wavelengths: %.6e, %.6e cm²\n",
                           ext_cross_interp[idx_low], ext_cross_interp[idx_high]);
                    printf("    Albedo at wavelengths: %.6f, %.6f\n",
                           albedo_interp[idx_low], albedo_interp[idx_high]);
                    printf("    Asymmetry at wavelengths: %.6f, %.6f\n",
                           asym_interp[idx_low], asym_interp[idx_high]);
                } else if (wavelength_rt_micron <= optics->wavelength[0]) {
                    printf("    λ_rt below lookup range - using first table value\n");
                    printf("    Extinction: %.6e cm² (from wavelength %.6e μm)\n",
                           ext_cross_interp[0], optics->wavelength[0]);
                } else {
                    printf("    λ_rt above lookup range - using last table value\n");
                    printf("    Extinction: %.6e cm² (from wavelength %.6e μm)\n",
                           ext_cross_interp[optics->n_wavelengths-1],
                           optics->wavelength[optics->n_wavelengths-1]);
                }
                    printf("    Interpolated extinction: %.6e cm²\n", ext_cross);
                    printf("    Interpolated albedo: %.6f\n", alb_val);
                    printf("    Interpolated asymmetry: %.6f\n", asym_val);
                    const char *species_name = get_species_name(cloud_species_id);
                    printf("    Final opacity c%s[%d][%d] = %.6e cm⁻¹ (n_density=%.6e * ext_cross=%.6e)\n",
                           species_name, j, i, particle_number_density_cm3 * ext_cross, particle_number_density_cm3, ext_cross);
                }
#endif
                
                // Calculate opacity (cm^-1) - matches MATLAB: crow = cloudden/VP * 1.0E-3 * pow(10, interp1(...))
                // Write to the correct arrays for this species
                c_out[j][i] = particle_number_density_cm3 * ext_cross;  // cm^-1
                a_out[j][i] = alb_val;
                g_out[j][i] = asym_val;
            }  // End RT wavelength loop
            
            // Debug: Print summary for wavelength interpolation
#if CLOUD_DEBUG_MIE
            if (particle_radius > 1e-10 && particle_number_density[j][cloud_species_idx] > 1e-10) {
                printf("\n  Wavelength interpolation complete for layer %d\n", j);
                printf("  Sample final RT values (first 5 wavelengths):\n");
                for (int i = 0; i < 5 && i < NLAMBDA; i++) {
                    const char *species_name = get_species_name(cloud_species_id);
                    printf("    λ=%.6e nm: c%s=%.6e cm⁻¹, a%s=%.6f, g%s=%.6f\n",
                           wavelength[i], species_name, c_out[j][i], species_name, a_out[j][i], species_name, g_out[j][i]);
                }
                printf("%s\n",fillmi); 
            }
#endif
        }  // End layer loop
        
        free(ext_cross_interp);
        free(albedo_interp);
        free(asym_interp);
    }  // End species loop
}

