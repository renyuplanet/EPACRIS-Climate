/*Function to read in collision-induce absorption cross section data */
/* output opacity in cm+5 */
/* Interpolation to specified temperature in each layer */
// Function is setup to hold opacities in cache for reinterpolation during the run

#include <string.h>
#include "config.h"
#include <stdio.h>
#include "global_temp.h"
#include "nrutil.h"
#include "routine.h"

// Forward declaration for Interpolation2D (defined in Interpolation.c)
double Interpolation2D(double x, double y, double xx[], int nx, double yy[], int ny, double **data);

// Simple CIA file structure
typedef struct {
	char filename[256];        // Filename
	char pair_name[32];        // Gas pair name (e.g., "H2-H2")
	double **data;             // CIA data array
	double wavelength[20000];  // Wavelength array (large enough for all files)
	int num_lines;             // Number of data lines
	int loaded;                // Whether data is loaded
} CIAFile;

// Static array of all CIA files we support
#define MAX_CIA_FILES 6
static CIAFile cia_files[MAX_CIA_FILES] = {
	{ "", "H2-H2", NULL, {0.0}, 0, 0 },
	{ "", "H2-He", NULL, {0.0}, 0, 0 },
	{ "", "H2-H", NULL, {0.0}, 0, 0 },
	{ "", "N2-N2", NULL, {0.0}, 0, 0 },
	{ "", "N2-H2", NULL, {0.0}, 0, 0 },
	{ "", "CO2-CO2", NULL, {0.0}, 0, 0 }
};

// Function prototypes
void readcia(void);
void reinterpolate_cia(void);
void cleanup_cia_cache(void);
void reinterpolate_all_cia_opacities(void);

// Helper function to get a CIA file by pair name
static CIAFile* get_cia_file(const char* pair_name) {
	for (int i = 0; i < MAX_CIA_FILES; i++) {
		if (strcmp(cia_files[i].pair_name, pair_name) == 0) {
			return &cia_files[i];
		}
	}
	return NULL;
}



// Main function to read CIA data
void readcia(void) {
	int i, j, s, nl;
	FILE *fim, *fout;
	char *temp1;
	double **cia;
	char crossfile[1024];
	double tfitting[zbin+1];
	char dataline[10000];
	
	// Temperature fitting for each layer
	for (i=1; i<=zbin; i++) {
		tfitting[i] = tl[i];
		if (tfitting[i] > 2000.0) tfitting[i] = 2000.0;
		if (tfitting[i] < 100.0) tfitting[i] = 100.0;
	}
	
	// Temperature bins for interpolation
	double temp[20] = {100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
					1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0,1800.0,1900.0,2000.0};
	
	printf("Loading CIA opacity data...\n");
	
	// Process each CIA pair
	// H2-H2 CIA data
	CIAFile* h2h2_file = get_cia_file("H2-H2");
	if (h2h2_file == NULL) {
		printf("Error: H2-H2 CIA file structure not found - skipping and initializing to zeros\n");
		// Initialize with zeros and continue processing other species
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2H2CIA[i][j] = 0.0;
			}
		}
	} else {
		if (!h2h2_file->loaded) {
			strcpy(crossfile, CROSSHEADING);
			strcat(crossfile, "H2-H2_CIA.dat");
			strcpy(h2h2_file->filename, crossfile);
			printf("Reading CIA file: %s\n", crossfile);
			
			fim = fopen(crossfile, "r");
			if (fim == NULL) {
				printf("Warning: Could not open CIA file: %s\n", crossfile);
				// Initialize with zeros if file not found
				for (i=1; i<=zbin; i++) {
					for (j=0; j<NLAMBDA; j++) {
						H2H2CIA[i][j] = 0.0;
					}
				}
			} else {
				s = LineNumber(fim, 20000);  // Increased buffer for larger files
				fclose(fim);
				nl = s-1;
				h2h2_file->num_lines = nl;
				
				// Allocate matrix and read data
				cia = dmatrix(0, nl-1, 0, 19);
				fim = fopen(crossfile, "r");
				temp1 = fgets(dataline, 10000, fim); /* Read in the header line */
				
				for (i=0; i<nl; i++) {
					fscanf(fim, "%lf", &h2h2_file->wavelength[i]);
					for (j=0; j<19; j++) {
						fscanf(fim, "%le", &cia[i][j]);
					}
					fscanf(fim, "%le\n", &cia[i][19]);
				}
				fclose(fim);
				
				// Save to cache and calculate cross sections
				h2h2_file->data = cia;
				h2h2_file->loaded = 1;
				
				for (i=1; i<=zbin; i++) {
					for (j=0; j<NLAMBDA; j++) {
						H2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i], 
										 h2h2_file->wavelength, nl, temp, 20, cia);
					}
				}
				
				// Write debug output
				fout = fopen("AuxillaryOut/CheckCIA_H2H2.dat", "w");
				for (j=0; j<NLAMBDA; j++) {
					fprintf(fout, "%2.6f\t", wavelength[j]);
					for (i=1; i<=zbin; i++) {
						fprintf(fout, "%2.6e\t", H2H2CIA[i][j]);
					}
					fprintf(fout, "\n");
				}
				fclose(fout);
			}
		} else {
			printf("Using cached data for H2-H2 CIA\n");
			// Recalculate using cached data
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					H2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i], 
									h2h2_file->wavelength, h2h2_file->num_lines, 
									temp, 20, h2h2_file->data);
				}
			}
		}
	}
	
	// H2-He CIA data
	CIAFile* h2he_file = get_cia_file("H2-He");
	if (h2he_file == NULL) {
		printf("Error: H2-He CIA file structure not found - skipping and initializing to zeros\n");
		// Initialize with zeros and continue processing other species
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HeCIA[i][j] = 0.0;
			}
		}
	} else if (!h2he_file->loaded) {
		strcpy(crossfile, CROSSHEADING);
		strcat(crossfile, "H2-He_CIA.dat");
		strcpy(h2he_file->filename, crossfile);
		printf("Reading CIA file: %s\n", crossfile);
		
		fim = fopen(crossfile, "r");
		if (fim == NULL) {
			printf("Warning: Could not open CIA file: %s\n", crossfile);
			// Initialize with zeros
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					H2HeCIA[i][j] = 0.0;
				}
			}
		} else {
			s = LineNumber(fim, 10000);
			fclose(fim);
			nl = s-1;
			h2he_file->num_lines = nl;
			// printf("  Found %d lines of data\n", nl);
			
			cia = dmatrix(0, nl-1, 0, 19);
			fim = fopen(crossfile, "r");
			temp1 = fgets(dataline, 10000, fim);
			
			for (i=0; i<nl; i++) {
				fscanf(fim, "%lf", &h2he_file->wavelength[i]);
				for (j=0; j<19; j++) {
					fscanf(fim, "%le", &cia[i][j]);
				}
				fscanf(fim, "%le\n", &cia[i][19]);
			}
			fclose(fim);
			
			h2he_file->data = cia;
			h2he_file->loaded = 1;
			
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					H2HeCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i], 
									h2he_file->wavelength, nl, temp, 20, cia);
				}
			}
			
			fout = fopen("AuxillaryOut/CheckCIA_H2He.dat", "w");
			for (j=0; j<NLAMBDA; j++) {
				fprintf(fout, "%2.6f\t", wavelength[j]);
				for (i=1; i<=zbin; i++) {
					fprintf(fout, "%2.6e\t", H2HeCIA[i][j]);
				}
				fprintf(fout, "\n");
			}
			fclose(fout);
		}
	} else {
		printf("Using cached data for H2-He CIA\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HeCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								h2he_file->wavelength, h2he_file->num_lines,
								temp, 20, h2he_file->data);
			}
		}
	}
	
	// H2-H CIA data
	CIAFile* h2h_file = get_cia_file("H2-H");
	if (h2h_file == NULL) {
		printf("Error: H2-H CIA file structure not found - skipping and initializing to zeros\n");
		// Initialize with zeros and continue processing other species
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HCIA[i][j] = 0.0;
			}
		}
	} else if (!h2h_file->loaded) {
		strcpy(crossfile, CROSSHEADING);
		strcat(crossfile, "H2-H_CIA.dat");
		strcpy(h2h_file->filename, crossfile);
		printf("Reading CIA file: %s\n", crossfile);
		
		fim = fopen(crossfile, "r");
		if (fim == NULL) {
			printf("Warning: Could not open CIA file: %s\n", crossfile);
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					H2HCIA[i][j] = 0.0;
				}
			}
		} else {
			s = LineNumber(fim, 10000);
			fclose(fim);
			nl = s-1;
			h2h_file->num_lines = nl;
			// printf("  Found %d lines of data\n", nl);
			
			cia = dmatrix(0, nl-1, 0, 19);
			fim = fopen(crossfile, "r");
			temp1 = fgets(dataline, 10000, fim);
			
			for (i=0; i<nl; i++) {
				fscanf(fim, "%lf", &h2h_file->wavelength[i]);
				for (j=0; j<19; j++) {
					fscanf(fim, "%le", &cia[i][j]);
				}
				fscanf(fim, "%le\n", &cia[i][19]);
			}
			fclose(fim);
			
			h2h_file->data = cia;
			h2h_file->loaded = 1;
			
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					H2HCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
									h2h_file->wavelength, nl, temp, 20, cia);
				}
			}
			
			fout = fopen("AuxillaryOut/CheckCIA_H2H.dat", "w");
			for (j=0; j<NLAMBDA; j++) {
				fprintf(fout, "%2.6f\t", wavelength[j]);
				for (i=1; i<=zbin; i++) {
					fprintf(fout, "%2.6e\t", H2HCIA[i][j]);
				}
				fprintf(fout, "\n");
			}
			fclose(fout);
		}
	} else {
		printf("Using cached data for H2-H CIA\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								h2h_file->wavelength, h2h_file->num_lines,
								temp, 20, h2h_file->data);
			}
		}
	}
	
	// N2-N2 CIA data
	CIAFile* n2n2_file = get_cia_file("N2-N2");
	if (n2n2_file == NULL) {
		printf("Error: N2-N2 CIA file structure not found - skipping and initializing to zeros\n");
		// Initialize with zeros and continue processing other species
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2N2CIA[i][j] = 0.0;
			}
		}
	} else if (!n2n2_file->loaded) {
		strcpy(crossfile, CROSSHEADING);
		strcat(crossfile, "N2-N2_CIA.dat");
		strcpy(n2n2_file->filename, crossfile);
		printf("Reading CIA file: %s\n", crossfile);
		
		fim = fopen(crossfile, "r");
		if (fim == NULL) {
			printf("Warning: Could not open CIA file: %s\n", crossfile);
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					N2N2CIA[i][j] = 0.0;
				}
			}
		} else {
			s = LineNumber(fim, 10000);
			fclose(fim);
			nl = s-1;
			n2n2_file->num_lines = nl;
			// printf("  Found %d lines of data\n", nl);
			
			cia = dmatrix(0, nl-1, 0, 19);
			fim = fopen(crossfile, "r");
			temp1 = fgets(dataline, 10000, fim);
			
			for (i=0; i<nl; i++) {
				fscanf(fim, "%lf", &n2n2_file->wavelength[i]);
				for (j=0; j<19; j++) {
					fscanf(fim, "%le", &cia[i][j]);
				}
				fscanf(fim, "%le\n", &cia[i][19]);
			}
			fclose(fim);
			
			n2n2_file->data = cia;
			n2n2_file->loaded = 1;
			
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					N2N2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
									n2n2_file->wavelength, nl, temp, 20, cia);
				}
			}
			
			fout = fopen("AuxillaryOut/CheckCIA_N2N2.dat", "w");
			for (j=0; j<NLAMBDA; j++) {
				fprintf(fout, "%2.6f\t", wavelength[j]);
				for (i=1; i<=zbin; i++) {
					fprintf(fout, "%2.6e\t", N2N2CIA[i][j]);
				}
				fprintf(fout, "\n");
			}
			fclose(fout);
		}
	} else {
		printf("Using cached data for N2-N2 CIA\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2N2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								n2n2_file->wavelength, n2n2_file->num_lines,
								temp, 20, n2n2_file->data);
			}
		}
	}
	
	// N2-H2 CIA data
	CIAFile* n2h2_file = get_cia_file("N2-H2");
	if (n2h2_file == NULL) {
		printf("Error: N2-H2 CIA file structure not found - skipping and initializing to zeros\n");
		// Initialize with zeros and continue processing other species
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2H2CIA[i][j] = 0.0;
			}
		}
	} else if (!n2h2_file->loaded) {
		strcpy(crossfile, CROSSHEADING);
		strcat(crossfile, "N2-H2_CIA.dat");
		strcpy(n2h2_file->filename, crossfile);
		printf("Reading CIA file: %s\n", crossfile);
		
		fim = fopen(crossfile, "r");
		if (fim == NULL) {
			printf("Warning: Could not open CIA file: %s\n", crossfile);
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					N2H2CIA[i][j] = 0.0;
				}
			}
		} else {
			s = LineNumber(fim, 10000);
			fclose(fim);
			nl = s-1;
			n2h2_file->num_lines = nl;
			// printf("  Found %d lines of data\n", nl);
			
			cia = dmatrix(0, nl-1, 0, 19);
			fim = fopen(crossfile, "r");
			temp1 = fgets(dataline, 10000, fim);
			
			for (i=0; i<nl; i++) {
				fscanf(fim, "%lf", &n2h2_file->wavelength[i]);
				for (j=0; j<19; j++) {
					fscanf(fim, "%le", &cia[i][j]);
				}
				fscanf(fim, "%le\n", &cia[i][19]);
			}
			fclose(fim);
			
			n2h2_file->data = cia;
			n2h2_file->loaded = 1;
			
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					N2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
									n2h2_file->wavelength, nl, temp, 20, cia);
				}
			}
			
			fout = fopen("AuxillaryOut/CheckCIA_N2H2.dat", "w");
			for (j=0; j<NLAMBDA; j++) {
				fprintf(fout, "%2.6f\t", wavelength[j]);
				for (i=1; i<=zbin; i++) {
					fprintf(fout, "%2.6e\t", N2H2CIA[i][j]);
				}
				fprintf(fout, "\n");
			}
			fclose(fout);
		}
	} else {
		printf("Using cached data for N2-H2 CIA\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								n2h2_file->wavelength, n2h2_file->num_lines,
								temp, 20, n2h2_file->data);
			}
		}
	}
	
	// CO2-CO2 CIA data
	CIAFile* co2co2_file = get_cia_file("CO2-CO2");
	if (co2co2_file == NULL) {
		printf("Error: CO2-CO2 CIA file structure not found - skipping and initializing to zeros\n");
		// Initialize with zeros and continue processing other species
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				CO2CO2CIA[i][j] = 0.0;
			}
		}
	} else if (!co2co2_file->loaded) {
		strcpy(crossfile, CROSSHEADING);
		strcat(crossfile, "CO2-CO2_CIA.dat");
		strcpy(co2co2_file->filename, crossfile);
		printf("Reading CIA file: %s\n", crossfile);
		
		fim = fopen(crossfile, "r");
		if (fim == NULL) {
			printf("Warning: Could not open CIA file: %s\n", crossfile);
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					CO2CO2CIA[i][j] = 0.0;
				}
			}
		} else {
			s = LineNumber(fim, 10000);
			fclose(fim);
			nl = s-1;
			co2co2_file->num_lines = nl;
			// printf("  Found %d lines of data\n", nl);
			
			cia = dmatrix(0, nl-1, 0, 19);
			fim = fopen(crossfile, "r");
			temp1 = fgets(dataline, 10000, fim);
			
			for (i=0; i<nl; i++) {
				fscanf(fim, "%lf", &co2co2_file->wavelength[i]);
				for (j=0; j<19; j++) {
					fscanf(fim, "%le", &cia[i][j]);
				}
				fscanf(fim, "%le\n", &cia[i][19]);
			}
			fclose(fim);
			
			co2co2_file->data = cia;
			co2co2_file->loaded = 1;
			
			for (i=1; i<=zbin; i++) {
				for (j=0; j<NLAMBDA; j++) {
					CO2CO2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
									 co2co2_file->wavelength, nl, temp, 20, cia);
				}
			}
			
			fout = fopen("AuxillaryOut/CheckCIA_CO2CO2.dat", "w");
			for (j=0; j<NLAMBDA; j++) {
				fprintf(fout, "%2.6f\t", wavelength[j]);
				for (i=1; i<=zbin; i++) {
					fprintf(fout, "%2.6e\t", CO2CO2CIA[i][j]);
				}
				fprintf(fout, "\n");
			}
			fclose(fout);
		}
	} else {
		printf("Using cached data for CO2-CO2 CIA\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				CO2CO2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								 co2co2_file->wavelength, co2co2_file->num_lines,
								 temp, 20, co2co2_file->data);
			}
		}
	}
	
}

// Reinterpolate CIA opacities for new temperature profile
void reinterpolate_cia(void) {
	int i, j;
	double temp[20] = {100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
					1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0,1800.0,1900.0,2000.0};
	
	double tfitting[zbin+1];
	for (i=1; i<=zbin; i++) {
		tfitting[i] = tl[i];
		if (tfitting[i] > 2000.0) tfitting[i] = 2000.0;
		if (tfitting[i] < 100.0) tfitting[i] = 100.0;
	}
	
	
	// Get each CIA file and reinterpolate if loaded
	CIAFile* h2h2_file = get_cia_file("H2-H2");
	if (h2h2_file != NULL && h2h2_file->loaded && h2h2_file->data) {
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								h2h2_file->wavelength, h2h2_file->num_lines,
								temp, 20, h2h2_file->data);
			}
		}
	} else {
		printf("Warning: H2-H2 CIA data not loaded, using zeros\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2H2CIA[i][j] = 0.0;
			}
		}
	}
	
	CIAFile* h2he_file = get_cia_file("H2-He");
	if (h2he_file != NULL && h2he_file->loaded && h2he_file->data) {
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HeCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								h2he_file->wavelength, h2he_file->num_lines,
								temp, 20, h2he_file->data);
			}
		}
	} else {
		printf("Warning: H2-He CIA data not loaded, using zeros\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HeCIA[i][j] = 0.0;
			}
		}
	}
	
	CIAFile* h2h_file = get_cia_file("H2-H");
	if (h2h_file != NULL && h2h_file->loaded && h2h_file->data) {
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HCIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								h2h_file->wavelength, h2h_file->num_lines,
								temp, 20, h2h_file->data);
			}
		}
	} else {
		printf("Warning: H2-H CIA data not loaded, using zeros\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				H2HCIA[i][j] = 0.0;
			}
		}
	}
	
	CIAFile* n2n2_file = get_cia_file("N2-N2");
	if (n2n2_file != NULL && n2n2_file->loaded && n2n2_file->data) {
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2N2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								n2n2_file->wavelength, n2n2_file->num_lines,
								temp, 20, n2n2_file->data);
			}
		}
	} else {
		printf("Warning: N2-N2 CIA data not loaded, using zeros\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2N2CIA[i][j] = 0.0;
			}
		}
	}
	
	CIAFile* n2h2_file = get_cia_file("N2-H2");
	if (n2h2_file != NULL && n2h2_file->loaded && n2h2_file->data) {
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2H2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								n2h2_file->wavelength, n2h2_file->num_lines,
								temp, 20, n2h2_file->data);
			}
		}
	} else {
		printf("Warning: N2-H2 CIA data not loaded, using zeros\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				N2H2CIA[i][j] = 0.0;
			}
		}
	}
	
	CIAFile* co2co2_file = get_cia_file("CO2-CO2");
	if (co2co2_file != NULL && co2co2_file->loaded && co2co2_file->data) {
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				CO2CO2CIA[i][j] = Interpolation2D(wavelength[j], tfitting[i],
								 co2co2_file->wavelength, co2co2_file->num_lines,
								 temp, 20, co2co2_file->data);
			}
		}
	} else {
		printf("Warning: CO2-CO2 CIA data not loaded, using zeros\n");
		for (i=1; i<=zbin; i++) {
			for (j=0; j<NLAMBDA; j++) {
				CO2CO2CIA[i][j] = 0.0;
			}
		}
	}
	
	//printf("CIA opacity reinterpolation complete\n");
}

// Clean up allocated memory
void cleanup_cia_cache(void) {
	printf("Cleaning up CIA opacity cache...\n");
	
	CIAFile* h2h2_file = get_cia_file("H2-H2");
	if (h2h2_file != NULL && h2h2_file->loaded && h2h2_file->data) {
		free_dmatrix(h2h2_file->data, 0, h2h2_file->num_lines-1, 0, 19);
		h2h2_file->data = NULL;
		h2h2_file->loaded = 0;
		h2h2_file->num_lines = 0;
		printf("Freed H2-H2 CIA cache\n");
	}
	
	CIAFile* h2he_file = get_cia_file("H2-He");
	if (h2he_file != NULL && h2he_file->loaded && h2he_file->data) {
		free_dmatrix(h2he_file->data, 0, h2he_file->num_lines-1, 0, 19);
		h2he_file->data = NULL;
		h2he_file->loaded = 0;
		h2he_file->num_lines = 0;
		printf("Freed H2-He CIA cache\n");
	}
	
	CIAFile* h2h_file = get_cia_file("H2-H");
	if (h2h_file != NULL && h2h_file->loaded && h2h_file->data) {
		free_dmatrix(h2h_file->data, 0, h2h_file->num_lines-1, 0, 19);
		h2h_file->data = NULL;
		h2h_file->loaded = 0;
		h2h_file->num_lines = 0;
		printf("Freed H2-H CIA cache\n");
	}
	
	CIAFile* n2n2_file = get_cia_file("N2-N2");
	if (n2n2_file != NULL && n2n2_file->loaded && n2n2_file->data) {
		free_dmatrix(n2n2_file->data, 0, n2n2_file->num_lines-1, 0, 19);
		n2n2_file->data = NULL;
		n2n2_file->loaded = 0;
		n2n2_file->num_lines = 0;
		printf("Freed N2-N2 CIA cache\n");
	}
	
	CIAFile* n2h2_file = get_cia_file("N2-H2");
	if (n2h2_file != NULL && n2h2_file->loaded && n2h2_file->data) {
		free_dmatrix(n2h2_file->data, 0, n2h2_file->num_lines-1, 0, 19);
		n2h2_file->data = NULL;
		n2h2_file->loaded = 0;
		n2h2_file->num_lines = 0;
		printf("Freed N2-H2 CIA cache\n");
	}
	
	CIAFile* co2co2_file = get_cia_file("CO2-CO2");
	if (co2co2_file != NULL && co2co2_file->loaded && co2co2_file->data) {
		free_dmatrix(co2co2_file->data, 0, co2co2_file->num_lines-1, 0, 19);
		co2co2_file->data = NULL;
		co2co2_file->loaded = 0;
		co2co2_file->num_lines = 0;
		printf("Freed CO2-CO2 CIA cache\n");
	}
	
	printf("CIA cache cleanup complete\n");
}

// Implementation for reinterpolate_all_cia_opacities
void reinterpolate_all_cia_opacities(void) {
	reinterpolate_cia();
}




