/*Function to read in cross section data */
/* input cross section in m^2 */
/* input wavelength in microns */
/* output mean cross section in cm^2 */
/* Interpolation to specified T and P */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "constant.h"
#include "config.h"
#include "nrutil.h"  // For f3tensor() and free_f3tensor()
#include <stdio.h>

#define MAX_SPECIES_CACHE 50  // Maximum number of species to cache

// External references to variables defined in the main file
extern const char* species[];  // External reference to species array
extern double **opacCO2, **opacO2, **opacSO2, **opacH2O, **opacOH, **opacH2CO;
extern double **opacH2O2, **opacHO2, **opacH2S, **opacCO, **opacO3, **opacCH4;
extern double **opacNH3;
extern double **opacC2H2, **opacC2H4, **opacC2H6, **opacHCN, **opacCH2O2, **opacHNO3;
extern double **opacN2O, **opacN2, **opacNO, **opacNO2, **opacOCS;
extern double tl[], pl[];  // External reference to temperature and pressure arrays
// zbin is defined in a header as a macro, not a variable
// Use CROSSHEADING_STR which is the string version of the CROSSHEADING macro
extern char CROSSHEADING_STR[];  // External reference to opacity file path

// NUM_SPECIES is defined in epacris_main.c where species[] is actually defined
// Cannot use sizeof() on extern array declarations
extern const int NUM_SPECIES;  // Number of opacity species

// Function prototypes
void readcross(char Fname[], double **xsc);
void reinterpolate_opacities_by_file(double **xsc, char Fname[]);
void read_all_opacities(void);
void cleanup_opacity_cache(void);
void reinterpolate_all_opacities(void);
// Forward declaration for Interpolation2D (defined in Interpolation.c)
double Interpolation2D(double x, double y, double xx[], int nx, double yy[], int ny, double **data);

// Simple structure to hold opacity data for one species
typedef struct {
	char filename[256];          // File path
	char species_name[32];       // Species name (e.g., "H2O")
	double ***opac;              // 3D opacity data array [wavelength][temperature][pressure]
	double wave[NLAMBDA];        // Wavelength array
	double temp[NTEMP];          // Temperature array
	double pres[NPRESSURE];      // Pressure array
	int is_loaded;               // Whether data is loaded
} OpacityFile;

// Static array of all opacity files we support
static OpacityFile opacity_files[MAX_SPECIES_CACHE] = {0};
static int num_cached_files = 0;

// Helper function to get an opacity file by species name or filename
static OpacityFile* get_opacity_file(const char* filename) {
	// Extract species name from filename
	const char* last_slash = strrchr(filename, '/');
	const char* file_base = last_slash ? last_slash + 1 : filename;
	const char* opac_prefix = strstr(file_base, "opac");
	
	char species_name[32] = {0};
	if (opac_prefix) {
		strncpy(species_name, opac_prefix + 4, sizeof(species_name) - 1);
		// Remove .dat extension if present
		char* dot = strchr(species_name, '.');
		if (dot) *dot = '\0';
	} else {
		strcpy(species_name, "unknown");
	}
	
	// Check if we already have this file cached
	for (int i = 0; i < num_cached_files; i++) {
		if (strcmp(opacity_files[i].species_name, species_name) == 0) {
			return &opacity_files[i];
		}
	}
	
	// If not found and we have room, create a new entry
	if (num_cached_files < MAX_SPECIES_CACHE) {
		OpacityFile* new_file = &opacity_files[num_cached_files++];
		strcpy(new_file->filename, filename);
		strcpy(new_file->species_name, species_name);
		new_file->opac = NULL;
		new_file->is_loaded = 0;
		return new_file;
	}
	
	// If we're out of cache slots, return the first one (least recently used)
	printf("Warning: Opacity cache full, reusing first slot\n");
	return &opacity_files[0];
}

// Helper function to get opacity array pointer for a species name
// Returns NULL if species not found
static double** get_opacity_array(const char *species_name) {
	// Ordered by frequency of use - most common species checked first
	if (strcmp(species_name, "H2O") == 0) return opacH2O;
	if (strcmp(species_name, "CO2") == 0) return opacCO2;
	if (strcmp(species_name, "O2") == 0) return opacO2;
	if (strcmp(species_name, "CH4") == 0) return opacCH4;
	if (strcmp(species_name, "CO") == 0) return opacCO;
	if (strcmp(species_name, "NH3") == 0) return opacNH3;
	if (strcmp(species_name, "N2") == 0) return opacN2;
	if (strcmp(species_name, "SO2") == 0) return opacSO2;
	if (strcmp(species_name, "H2S") == 0) return opacH2S;
	if (strcmp(species_name, "OH") == 0) return opacOH;
	if (strcmp(species_name, "H2CO") == 0) return opacH2CO;
	if (strcmp(species_name, "H2O2") == 0) return opacH2O2;
	if (strcmp(species_name, "HO2") == 0) return opacHO2;
	if (strcmp(species_name, "O3") == 0) return opacO3;
	if (strcmp(species_name, "OCS") == 0) return opacOCS;
	if (strcmp(species_name, "C2H2") == 0) return opacC2H2;
	if (strcmp(species_name, "C2H4") == 0) return opacC2H4;
	if (strcmp(species_name, "C2H6") == 0) return opacC2H6;
	if (strcmp(species_name, "CH2O2") == 0) return opacCH2O2;
	if (strcmp(species_name, "HCN") == 0) return opacHCN;
	if (strcmp(species_name, "HNO3") == 0) return opacHNO3;
	if (strcmp(species_name, "N2O") == 0) return opacN2O;
	if (strcmp(species_name, "NO") == 0) return opacNO;
	if (strcmp(species_name, "NO2") == 0) return opacNO2;
	return NULL;  // Species not found
}

// Read all opacities from the species list
void read_all_opacities(void) {
	char crossfile[1024];
	double **opac_ptr;

	printf("Reading opacities for %d gas species\n", NUM_SPECIES); 
	for (int i = 0; i < NUM_SPECIES; i++) {
		strcpy(crossfile, CROSSHEADING_STR);
		strcat(crossfile, "opac");
		strcat(crossfile, species[i]);
		strcat(crossfile, ".dat");
		
		// Get opacity array pointer for this species
		opac_ptr = get_opacity_array(species[i]);
		if (opac_ptr != NULL) {
			readcross(crossfile, opac_ptr);
		} else {
			printf("Warning: No opacity array found for species %s\n", species[i]);
		}
	}
}

// Function to clean up opacity cache
void cleanup_opacity_cache(void) {
	printf("Cleaning up opacity cache...\n");
	
	for (int i = 0; i < num_cached_files; i++) {
		if (opacity_files[i].is_loaded && opacity_files[i].opac != NULL) {
			free_f3tensor(opacity_files[i].opac, 0, NLAMBDA-1, 0, NTEMP-1, 0, NPRESSURE-1);
			opacity_files[i].opac = NULL;
			opacity_files[i].is_loaded = 0;
			printf("Freed opacity cache for %s\n", opacity_files[i].species_name);
		}
	}
	
	num_cached_files = 0;
	printf("Opacity cache cleanup complete\n");
}

// Function to reinterpolate all opacities
// Called every 10 RT steps to reinterpolate opacities as T/P change
void reinterpolate_all_opacities(void) {
	char crossfile[1024];
	double **opac_ptr;

	for (int i = 0; i < NUM_SPECIES; i++) {
		strcpy(crossfile, CROSSHEADING_STR);
		strcat(crossfile, "opac");
		strcat(crossfile, species[i]);
		strcat(crossfile, ".dat");
		
		// Get opacity array pointer for this species (reuses optimized lookup)
		opac_ptr = get_opacity_array(species[i]);
		if (opac_ptr != NULL) {
			reinterpolate_opacities_by_file(opac_ptr, crossfile);
		}
		// Silently skip if species not found (already warned in read_all_opacities)
	}
}

// Main function to read opacity data for a specific file
// Now opacities are cached in the OpacityFile structure
// Does not need to be read in every time opacities are reinterpolated
void readcross(char Fname[], double **xsc)
{
	int i, j, k;
	FILE *fim;
	double wave1[NLAMBDA]; /* in nm */
	double presdummy, cross[NLAMBDA];
	double pfitting, tfitting;
	
	// Get the cache entry for this species
	OpacityFile* opac_file = get_opacity_file(Fname);

	// Check if we need to load the data
	if (!opac_file->is_loaded) {
		printf("Reading opacity file: %s\n", Fname);
		
		// Allocate opacity tensor if first time
		if (opac_file->opac == NULL) {
			opac_file->opac = f3tensor(0, NLAMBDA-1, 0, NTEMP-1, 0, NPRESSURE-1);
		}
		
		fim = fopen(Fname, "r");
		if (!fim) {
			printf("Error: Could not open opacity file %s\n", Fname);
			// Initialize with zeros if file not found
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					xsc[i][j] = 0.0;
				}
			}
			return;
		}
		
		/* Header Lines */
		for (i=0; i<NTEMP-1; i++) {
			fscanf(fim, "%lf", opac_file->temp+i);
		}
		i=NTEMP-1;
		fscanf(fim, "%lf\n", opac_file->temp+i);
		
		for (i=0; i<NPRESSURE-1; i++) {
			fscanf(fim, "%le", opac_file->pres+i);
		}
		i=NPRESSURE-1;
		fscanf(fim, "%le\n", opac_file->pres+i);
		
		/* Read in data */
		for (i=0; i<NLAMBDA; i++) {
			fscanf(fim, "%le\n", opac_file->wave+i);
			for (j=0; j<NPRESSURE; j++) {
				fscanf(fim, "%le", &presdummy);
				for (k=0; k<NTEMP-1; k++) {
					fscanf(fim, "%le", &opac_file->opac[i][k][j]);
				}
				k=NTEMP-1;
				fscanf(fim, "%le\n", &opac_file->opac[i][k][j]);
			}
		}
		
		fclose(fim);
		opac_file->is_loaded = 1;
	} else {
		printf("Using cached opacity data for: %s\n", opac_file->species_name);
	}
	
	/* get wave1 from wave */
	for (i=0; i<NLAMBDA; i++) {
		wave1[i] = opac_file->wave[i]*1.0E+9; /* convert to nm */
	}
	
	/* Calculate the cross sections */
	for (i=1; i<=zbin; i++) {
		for (j=0; j<NLAMBDA; j++) {
			if (pl[i] > opac_file->pres[NPRESSURE-1]) {
				pfitting = opac_file->pres[NPRESSURE-1];
			} else if (pl[i] < opac_file->pres[0]) {
				pfitting = opac_file->pres[0];
			} else {
				pfitting = pl[i];
			}
			
			if (tl[i] > opac_file->temp[NTEMP-1]) {
				tfitting = opac_file->temp[NTEMP-1];
			} else if (tl[i] < opac_file->temp[0]) {
				tfitting = opac_file->temp[0];
			} else {
				tfitting = tl[i];
			}
			
			cross[j] = fmax(Interpolation2D(tfitting, pfitting, 
										  opac_file->temp, NTEMP, 
										  opac_file->pres, NPRESSURE, 
										  *(opac_file->opac+j))*1.0E+4, 0.0);
			xsc[i][j] = cross[j];
		}
	}
}

// Function to reinterpolate opacity data without re-reading the file
void reinterpolate_opacities_by_file(double **xsc, char Fname[]) {
	// Find the cache entry for this species
	OpacityFile* opac_file = get_opacity_file(Fname);
	
	if (!opac_file->is_loaded || !opac_file->opac) {
		// If data not loaded, load it
		readcross(Fname, xsc); // Should never happen
		return;
	}

	int i, j;
	double wave1[NLAMBDA]; /* in nm */
	double pfitting, tfitting;
	double cross[NLAMBDA];

	/* get wave1 from wave */
	for (i=0; i<NLAMBDA; i++) {
		wave1[i] = opac_file->wave[i]*1.0E+9; /* convert to nm */
	}
	
	/* Calculate the cross sections */
	for (i=1; i<=zbin; i++) {
		for (j=0; j<NLAMBDA; j++) {
			if (pl[i] > opac_file->pres[NPRESSURE-1]) {
				pfitting = opac_file->pres[NPRESSURE-1];
			} else if (pl[i] < opac_file->pres[0]) {
				pfitting = opac_file->pres[0];
			} else {
				pfitting = pl[i];
			}
			
			if (tl[i] > opac_file->temp[NTEMP-1]) {
				tfitting = opac_file->temp[NTEMP-1];
			} else if (tl[i] < opac_file->temp[0]) {
				tfitting = opac_file->temp[0];
			} else {
				tfitting = tl[i];
			}
			
			cross[j] = fmax(Interpolation2D(tfitting, pfitting, 
										  opac_file->temp, NTEMP, 
										  opac_file->pres, NPRESSURE, 
										  *(opac_file->opac+j))*1.0E+4, 0.0);
			xsc[i][j] = cross[j];
		}
	}
}
